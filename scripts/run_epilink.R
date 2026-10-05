#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
args <- setdiff(args, "--overwrite")
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/run_epilink.R <pbmc|bmmc> <celltype> <poisson|nb> [--overwrite]", call. = FALSE)
}
dataset_name <- args[[1]]
celltype <- args[[2]]
family <- args[[3]]
if (!family %in% c("poisson", "nb")) stop("family must be poisson or nb", call. = FALSE)

project_root <- normalizePath(file.path(dirname(commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
  |> sub("^--file=", "", x = _)), ".."), mustWork = TRUE)
setwd(project_root)

if (!requireNamespace("yaml", quietly = TRUE)) stop("install the yaml package", call. = FALSE)
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(source_file)
configuration <- yaml::read_yaml("config/datasets.yml")
if (!dataset_name %in% names(configuration)) {
  stop("unknown dataset: ", dataset_name, call. = FALSE)
}
settings <- utils::modifyList(configuration$defaults, configuration[[dataset_name]])
inputs <- unlist(settings[c("rna", "atac", "metadata", "candidate_pairs")], use.names = TRUE)
missing_inputs <- names(inputs)[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("missing configured inputs: ", paste(paste0(missing_inputs, "=", inputs[missing_inputs]), collapse = ", "), call. = FALSE)
}

output_dir <- file.path(settings$output_dir, paste0("EpiLink_", family, "_s"), celltype)
output_file <- file.path(output_dir, "epilink_result.rds")
if (file.exists(output_file) && !overwrite) {
  stop("output already exists: ", output_file, "; use --overwrite to replace it", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading configured inputs for ", dataset_name, " ...")
rna <- readRDS(settings$rna)
atac <- readRDS(settings$atac)
metadata <- readRDS(settings$metadata)
pairs <- readRDS(settings$candidate_pairs)
pairs$peak <- CanonicalPeakNames(pairs$peak)
if (!is.null(settings$metadata_barcode_column)) {
  rownames(metadata) <- metadata[[settings$metadata_barcode_column]]
}
if (!identical(colnames(rna), colnames(atac))) stop("RNA and ATAC cell order differs", call. = FALSE)
if (!all(colnames(rna) %in% rownames(metadata))) stop("metadata is missing matrix cells", call. = FALSE)
metadata <- metadata[colnames(rna), , drop = FALSE]

if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Signac", quietly = TRUE)) {
  stop("Seurat and Signac are required for the real-data runner", call. = FALSE)
}
object <- Seurat::CreateSeuratObject(counts = rna, meta.data = metadata, assay = "RNA")
object[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = atac, sep = InferPeakSeparator(rownames(atac))
)

started <- Sys.time()
result <- EpiLink(
  object = object,
  gene_peak_df = pairs,
  celltype_col = settings$celltype_column,
  celltype = celltype,
  covariates = unlist(settings$covariates, use.names = FALSE),
  family = family,
  null_generation = "s",
  stat = settings$statistic,
  fdrcutoff = settings$fdr,
  bin = settings$binarize_atac,
  selection = settings$selection,
  correction = settings$correction,
  conservative = settings$conservative,
  ncore = settings$n_cores,
  returnObject = TRUE,
  seed = settings$seed
)
saveRDS(result, output_file)
all_results <- result@misc$EpiLink$res_df
significant <- EpiLinkFilter(result, stat = settings$statistic,
                             fdrcutoff = settings$fdr,
                             correction = settings$correction,
                             conservative = settings$conservative,
                             verbose = FALSE)
write.table(all_results, file.path(output_dir, "all_pairs.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(significant, file.path(output_dir, "significant_pairs.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(result@misc$EpiLink$threshold_audit,
            file.path(output_dir, "fdr_threshold_audit.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
metadata_run <- data.frame(
  dataset = dataset_name, celltype = celltype, family = family,
  null_generation = "s", started = as.character(started),
  finished = as.character(Sys.time()), R = R.version.string,
  seed = settings$seed, statistic = settings$statistic, fdr = settings$fdr,
  correction = settings$correction,
  conservative = settings$conservative,
  threshold = result@misc$EpiLink$threshold, n_tested = nrow(all_results),
  n_significant = nrow(significant)
)
write.table(metadata_run, file.path(output_dir, "run_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Finished. Results: ", output_dir)
