#!/usr/bin/env Rscript

# Build Sakaue-style dataset-level annotations from the stored EpiLink objects.
# All available EpiLink cell types in a dataset are unioned into one BED.  The
# model is not re-fit: correction=TRUE/FALSE is applied by EpiLinkFilter() to
# the same stored pair-level object, using the exact finite-cell correction in
# R/EpiLink.R.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
fdr <- if (length(args) >= 2L) as.numeric(args[[2]]) else 0.1
gwas <- file.path(root, "results", "validation", "GWAS")
base <- file.path(gwas, "SCENT_style", "UKBB", "sakaue_style")
annotation_dir <- file.path(base, "annotations")
dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(root, "R", "EpiLink.R"))
datasets <- fread(file.path(gwas, "SCENT_style", "config", "datasets.tsv"))

import_chain_gz <- function(path) {
  tmp <- tempfile(fileext = ".over.chain")
  con <- gzfile(path, open = "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  writeLines(lines, tmp)
  import.chain(tmp)
}
chain_19_38 <- import_chain_gz(file.path(gwas, "reference", "liftover", "hg19ToHg38.over.chain.gz"))

safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

parse_peaks <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(data.table(chr = character(), start = integer(), end = integer()))
  z <- tstrsplit(x, "-", fixed = TRUE)
  if (length(z) != 3L) stop("Peak IDs must have chr-start-end form")
  out <- data.table(chr = z[[1]], start = as.integer(z[[2]]), end = as.integer(z[[3]]))
  out <- out[grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", chr) & !is.na(start) & !is.na(end) & start <= end]
  unique(out)
}

reduce_bed <- function(x) {
  if (!nrow(x)) return(x)
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand = TRUE)
  data.table(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
}

lift_bed <- function(x) {
  if (!nrow(x)) return(x)
  gr <- GRanges(x$chr, IRanges(x$start, x$end))
  lifted <- unlist(liftOver(gr, chain_19_38), use.names = FALSE)
  if (!length(lifted)) return(data.table(chr = character(), start = integer(), end = integer()))
  reduce_bed(data.table(chr = as.character(seqnames(lifted)), start = start(lifted), end = end(lifted)))
}

write_bed <- function(x, path) {
  y <- copy(x)
  if (nrow(y)) y[, start0 := pmax(start - 1L, 0L)] else y[, start0 := integer()]
  fwrite(y[, .(chr, start0, end)], path, sep = "\t", col.names = FALSE)
}

object_path <- function(dataset, celltype, family) {
  file.path(root, "results", "real", dataset, celltype,
            paste0("epilink_", family), "epilink_object.rds")
}

families <- c("nb", "poisson")
corrections <- c(TRUE, FALSE)
all_datasets <- unique(datasets$dataset)
manifest <- list()
celltype_manifest <- list()
k <- 0L
q <- 0L

for (ds in all_datasets) {
  ds_rows <- datasets[dataset == ds]
  # Use the primary cell types defined for the current validation analysis.
  # These are still unioned within each dataset, matching the Sakaue Fig. 3
  # dataset-level view while avoiding exploratory cell types.
  ds_rows <- ds_rows[interpretation == "primary"]
  available <- ds_rows[vapply(celltype, function(ct) {
    all(file.exists(object_path(ds, ct, families)))
  }, logical(1))]
  if (!nrow(available)) next
  message("Dataset ", ds, ": ", paste(available$celltype, collapse = ", "))

  for (family in families) {
    peak_by_correction <- setNames(vector("list", length(corrections)), paste0("correction_", c("T", "F")))
    pair_counts <- setNames(integer(length(corrections)), names(peak_by_correction))
    for (ct in available$celltype) {
      path <- object_path(ds, ct, family)
      obj <- readRDS(path)
      for (correction in corrections) {
        sig <- EpiLinkFilter(obj, stat = "zval", fdrcutoff = fdr,
                             correction = correction, conservative = FALSE,
                             verbose = FALSE)
        peaks <- if (nrow(sig) && "peak" %in% names(sig)) unique(sig$peak) else character()
        key <- paste0("correction_", if (correction) "T" else "F")
        peak_by_correction[[key]] <- c(peak_by_correction[[key]], peaks)
        pair_counts[[key]] <- pair_counts[[key]] + length(peaks)
      }
      q <- q + 1L
      celltype_manifest[[q]] <- data.table(
        dataset = ds, celltype = ct, family = family,
        object = path,
        n_pairs = nrow(obj@misc$EpiLink$res_df),
        n_sig_correction_T = length(unique(peak_by_correction[["correction_T"]])),
        n_sig_correction_F = length(unique(peak_by_correction[["correction_F"]]))
      )
      rm(obj)
      gc(verbose = FALSE)
    }

    for (correction in corrections) {
      corr_key <- paste0("correction_", if (correction) "T" else "F")
      native <- reduce_bed(parse_peaks(unique(peak_by_correction[[corr_key]])))
      hg38 <- if (all(available$build == "hg38")) native else {
        pieces <- list()
        for (ct in available$celltype) {
          obj <- readRDS(object_path(ds, ct, family))
          sig <- EpiLinkFilter(obj, stat = "zval", fdrcutoff = fdr,
                               correction = correction, conservative = FALSE,
                               verbose = FALSE)
          peaks <- if (nrow(sig) && "peak" %in% names(sig)) unique(sig$peak) else character()
          pieces[[ct]] <- if (available[celltype == ct, build] == "hg19") {
            lift_bed(parse_peaks(peaks))
          } else parse_peaks(peaks)
          rm(obj)
          gc(verbose = FALSE)
        }
        reduce_bed(rbindlist(pieces, fill = TRUE))
      }
      hg19 <- if (all(available$build == "hg19")) native else {
        # This is only an audit BED; the SCENT analysis itself uses hg38.
        reduce_bed(rbindlist(lapply(available$celltype, function(ct) {
          obj <- readRDS(object_path(ds, ct, family))
          sig <- EpiLinkFilter(obj, stat = "zval", fdrcutoff = fdr,
                               correction = correction, conservative = FALSE,
                               verbose = FALSE)
          peaks <- if (nrow(sig) && "peak" %in% names(sig)) unique(sig$peak) else character()
          rm(obj)
          if (available[celltype == ct, build] == "hg19") parse_peaks(peaks) else {
            # There is no need to create a reverse-lifted bed for the main
            # analysis; keep hg19 equal to native for mixed-build audit rows.
            data.table(chr = character(), start = integer(), end = integer())
          }
        }), fill = TRUE))
      }

      id <- paste(safe_id(ds), paste0("epilink_", family), corr_key, sep = "__")
      write_bed(hg38, file.path(annotation_dir, paste0(id, ".hg38.bed")))
      k <- k + 1L
      manifest[[k]] <- data.table(
        annotation_id = id, dataset = ds, celltype = "ALL_AVAILABLE",
        n_celltypes = nrow(available), celltypes = paste(available$celltype, collapse = ";"),
        method = paste0("epilink_", family), method_label = paste0("EpiLink-", toupper(family)),
        correction = correction, fdr = fdr, source_build = paste(unique(available$build), collapse = ";"),
        n_pairs = pair_counts[[corr_key]], n_native_intervals = nrow(native), n_hg38_intervals = nrow(hg38),
        n_hg19_intervals = nrow(hg19), lift_fraction_hg38 = if (nrow(native)) nrow(hg38) / nrow(native) else NA_real_
      )
    }
  }
}

fwrite(rbindlist(manifest, fill = TRUE), file.path(base, "sakaue_annotation_manifest.tsv"), sep = "\t")
fwrite(rbindlist(celltype_manifest, fill = TRUE), file.path(base, "sakaue_celltype_audit.tsv"), sep = "\t")
message("Wrote ", length(manifest), " dataset-level annotations under ", annotation_dir)
