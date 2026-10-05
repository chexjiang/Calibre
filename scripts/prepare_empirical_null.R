#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
if (length(setdiff(args, "--overwrite"))) {
  stop("Usage: Rscript scripts/prepare_empirical_null.R [--overwrite]",
       call. = FALSE)
}
file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)
if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml is required", call. = FALSE)
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(source_file)
settings <- yaml::read_yaml("config/empirical_null.yml")
dataset_configuration <- yaml::read_yaml("config/datasets.yml")
dataset_settings <- utils::modifyList(
  dataset_configuration$defaults,
  dataset_configuration[[settings$dataset]]
)
if (!requireNamespace("Seurat", quietly = TRUE) ||
    !requireNamespace("Signac", quietly = TRUE) ||
    !requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE) ||
    !requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  stop("Seurat, Signac, EnsDb.Hsapiens.v86, and hg38 BSgenome are required",
       call. = FALSE)
}

output_dir <- settings$output_dir
setup_dir <- file.path(output_dir, "setup")
dir.create(setup_dir, recursive = TRUE, showWarnings = FALSE)
for (folder in seq_len(ceiling(settings$n_replicates /
                               settings$replicates_per_folder))) {
  dir.create(file.path(output_dir, sprintf("%02d", folder)),
             recursive = TRUE, showWarnings = FALSE)
}
message("Loading PBMC/Tnk data")
rna <- readRDS(dataset_settings$rna)
atac <- readRDS(dataset_settings$atac)
metadata <- readRDS(dataset_settings$metadata)
if (!identical(colnames(rna), colnames(atac))) {
  stop("RNA and ATAC cell order differs", call. = FALSE)
}
metadata <- metadata[colnames(rna), , drop = FALSE]
keep_cells <- rownames(metadata)[
  metadata[[settings$celltype_column]] == settings$celltype
]
rna <- rna[, keep_cells, drop = FALSE]
atac <- atac[, keep_cells, drop = FALSE]
metadata <- metadata[keep_cells, , drop = FALSE]

annotation <- Signac::GetGRangesFromEnsDb(
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
)
annotation_levels <- GenomeInfoDb::seqlevels(annotation)
GenomeInfoDb::seqlevels(annotation) <- ifelse(
  annotation_levels == "MT", "chrM", paste0("chr", annotation_levels)
)
object <- Seurat::CreateSeuratObject(counts = rna, meta.data = metadata)
object[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = atac, sep = InferPeakSeparator(rownames(atac)), annotation = annotation
)

pairs_file <- file.path(setup_dir, "candidate_pairs.rds")
legacy_pairs <- file.path(output_dir, "candidate_pairs.rds")
if (file.exists(pairs_file) && !overwrite) {
  pairs <- readRDS(pairs_file)
} else if (file.exists(legacy_pairs) && !overwrite) {
  message("Reusing the fixed candidate pairs created by the earlier pilot")
  pairs <- readRDS(legacy_pairs)
  saveRDS(pairs, pairs_file)
} else {
  pairs <- sample_empirical_pairs(
    object, readRDS(dataset_settings$candidate_pairs),
    n_pairs = settings$n_pairs, selection = settings$selection,
    seed = settings$seed
  )
  saveRDS(pairs, pairs_file)
}
if (nrow(pairs) != settings$n_pairs) {
  stop("Prepared candidate-pair count does not match config", call. = FALSE)
}
utils::write.table(
  pairs, file.path(setup_dir, "candidate_pairs.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)

peak_metadata_file <- file.path(setup_dir, "signac_peak_metadata.rds")
legacy_peak_metadata <- file.path(output_dir, "signac_peak_metadata.rds")
if (file.exists(peak_metadata_file) && !overwrite) {
  peak_metadata <- readRDS(peak_metadata_file)
} else if (file.exists(legacy_peak_metadata) && !overwrite) {
  message("Reusing hg38 peak statistics created by the earlier pilot")
  peak_metadata <- readRDS(legacy_peak_metadata)
  saveRDS(peak_metadata, peak_metadata_file)
} else {
  message("Computing hg38 peak sequence statistics for Signac")
  object[["ATAC"]] <- Signac::RegionStats(
    object[["ATAC"]],
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    sep = c("-", "-")
  )
  peak_metadata <- object[["ATAC"]]@meta.features
  saveRDS(peak_metadata, peak_metadata_file)
}

permutation_file <- file.path(setup_dir, "permutation_indices.rds")
manifest_file <- file.path(setup_dir, "permutation_manifest.tsv")
replicate_seeds <- epilink_pair_seeds(
  settings$n_replicates, settings$seed + 1000L
)
if (file.exists(permutation_file) && !overwrite) {
  permutation_indices <- readRDS(permutation_file)
  if (!identical(dim(permutation_indices),
                 c(length(keep_cells), as.integer(settings$n_replicates)))) {
    stop("Stored permutation dimensions do not match config; use --overwrite",
         call. = FALSE)
  }
} else {
  message("Generating and storing ", settings$n_replicates,
          " deterministic ATAC cell permutations")
  permutation_indices <- vapply(
    replicate_seeds,
    function(seed) with_epilink_seed(seed, sample.int(length(keep_cells))),
    integer(length(keep_cells))
  )
  colnames(permutation_indices) <- sprintf(
    "replicate_%04d", seq_len(settings$n_replicates)
  )
  saveRDS(permutation_indices, permutation_file)
}
saveRDS(keep_cells, file.path(setup_dir, "cell_names.rds"))
manifest <- data.frame(
  replicate = seq_len(settings$n_replicates),
  folder = vapply(
    seq_len(settings$n_replicates), empirical_null_folder, character(1),
    replicates_per_folder = settings$replicates_per_folder
  ),
  permutation_seed = replicate_seeds,
  stringsAsFactors = FALSE
)
utils::write.table(manifest, manifest_file, sep = "\t", quote = FALSE,
                   row.names = FALSE)
metadata_out <- data.frame(
  dataset = settings$dataset, celltype = settings$celltype,
  n_cells = length(keep_cells), n_pairs = nrow(pairs),
  n_replicates = settings$n_replicates,
  replicates_per_folder = settings$replicates_per_folder,
  permutation = "one shared unrestricted within-Tnk cell permutation",
  rna_fixed = TRUE, seed = settings$seed,
  prepared = as.character(Sys.time()), R = R.version.string,
  stringsAsFactors = FALSE
)
utils::write.table(metadata_out, file.path(setup_dir, "preparation_metadata.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("Empirical-null preparation complete: ", setup_dir)
