#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(AnnotationDbi)
  library(EnsDb.Hsapiens.v86)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args)) args[[1L]] else "config/eqtl_validation.yml"
config <- yaml::read_yaml(config_path)
raw_dir <- config$reference$raw_dir
processed_dir <- config$reference$processed_dir
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

mapping_path <- file.path(processed_dir, "ensembl_to_symbol.rds")
if (file.exists(mapping_path)) {
  mapping <- as.data.table(readRDS(mapping_path))
} else {
  message("Building Ensembl-to-HGNC map from EnsDb.Hsapiens.v86...")
  mapping <- AnnotationDbi::select(
    EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    keys = AnnotationDbi::keys(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                              keytype = "GENEID"),
    keytype = "GENEID", columns = c("GENEID", "SYMBOL")
  )
  mapping <- as.data.table(mapping)[!is.na(SYMBOL) & SYMBOL != ""]
  mapping <- unique(mapping, by = "GENEID")
  saveRDS(mapping, mapping_path, compress = FALSE)
}

metadata_path <- file.path(raw_dir, "eqtl_catalogue_metadata.tsv")
metadata <- fread(metadata_path)
manifest <- list()

for (group in names(config$reference_groups)) {
  message("Preparing reference group: ", group)
  ids <- unlist(config$reference_groups[[group]], use.names = FALSE)
  output_path <- file.path(processed_dir, paste0("OneK1K_", group, ".rds"))
  if (file.exists(output_path)) {
    x <- readRDS(output_path)
    manifest[[length(manifest) + 1L]] <- data.table(
      resource = paste0("OneK1K_", group), path = output_path,
      n_rows = nrow(x), n_variants = uniqueN(x$variant_id),
      n_genes = uniqueN(x$gene), source_ids = paste(ids, collapse = ",")
    )
    next
  }
  tables <- lapply(ids, function(dataset_id) {
    path <- file.path(raw_dir, "eqtl_catalogue_onek1k",
                      paste0(dataset_id, ".credible_sets.tsv.gz"))
    if (!file.exists(path)) stop("Missing reference file: ", path)
    x <- fread(path, select = c("gene_id", "variant", "rsid", "pip"),
               showProgress = FALSE)
    x[, dataset_id := dataset_id]
    x
  })
  x <- rbindlist(tables)
  x[, gene_id := sub("\\..*$", "", gene_id)]
  x <- merge(x, mapping, by.x = "gene_id", by.y = "GENEID")
  parsed <- tstrsplit(x$variant, "_", fixed = TRUE, keep = 1:2)
  x[, chr := ifelse(grepl("^chr", parsed[[1L]]), parsed[[1L]],
                    paste0("chr", parsed[[1L]]))]
  x[, start := as.integer(parsed[[2L]])]
  x[, end := start]
  setnames(x, "SYMBOL", "gene")
  x <- x[!is.na(start) & !is.na(pip)]
  # A pair can appear in several cell states or credible sets.  The broad-cell
  # benchmark asks whether it is fine-mapped in any matched state.
  x <- x[order(-pip), .(
    chr = chr[[1L]], start = start[[1L]], end = end[[1L]],
    rsid = rsid[[1L]], pip = max(pip),
    source_datasets = paste(sort(unique(dataset_id)), collapse = ",")
  ), by = .(variant_id = variant, gene)]
  setcolorder(x, c("chr", "start", "end", "variant_id", "rsid", "gene",
                   "pip", "source_datasets"))
  saveRDS(as.data.frame(x), output_path, compress = FALSE)
  manifest[[length(manifest) + 1L]] <- data.table(
    resource = paste0("OneK1K_", group),
    path = output_path,
    n_rows = nrow(x),
    n_variants = uniqueN(x$variant_id),
    n_genes = uniqueN(x$gene),
    source_ids = paste(ids, collapse = ",")
  )
}

tss_path <- file.path(raw_dir, "pgboost", "TSS.hg38.tsv.gz")
tss <- fread(tss_path)
setnames(tss, c("gene", "chr", "tss"))
tss[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
tss <- unique(tss[!is.na(gene) & !is.na(tss), .(gene, chr, tss = as.integer(tss))])
saveRDS(as.data.frame(tss), file.path(processed_dir, "TSS.hg38.rds"),
        compress = "gzip")

background_output <- file.path(processed_dir, "1000G_MAC5_SNPs.hg38.rds")
if (!file.exists(background_output)) {
  message("Preparing binary cache for 1000 Genomes background SNPs...")
  background_path <- file.path(raw_dir, "pgboost",
                               "1000G_MAC5_SNPs.hg38.bed.gz")
  background <- fread(background_path, header = FALSE, showProgress = FALSE)
  background <- background[, .(
    chr = ifelse(grepl("^chr", V1), V1, paste0("chr", V1)),
    start = as.integer(V2) + 1L, end = as.integer(V3)
  )]
  background[, variant_id := paste(chr, start, end, sep = ":")]
  saveRDS(background, background_output, compress = FALSE)
}

raw_files <- c(
  Sys.glob(file.path(raw_dir, "eqtl_catalogue_onek1k", "*.gz")),
  Sys.glob(file.path(raw_dir, "pgboost", "*.gz")), metadata_path
)
checksums <- tools::md5sum(raw_files)
checksum_table <- data.table(path = names(checksums), md5 = unname(checksums))
fwrite(rbindlist(manifest), file.path(processed_dir, "reference_manifest.tsv"),
       sep = "\t")
fwrite(checksum_table, file.path(processed_dir, "raw_checksums.tsv"), sep = "\t")
message("Prepared ", length(config$reference_groups),
        " broad-cell fine-mapping references in ", processed_dir)
