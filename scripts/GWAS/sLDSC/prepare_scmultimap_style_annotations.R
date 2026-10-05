#!/usr/bin/env Rscript

# Build the linked-peak and all-tested-peak annotations used by the
# scMultiMap-style S-LDSC benchmark.  In Su et al., linked peaks were selected
# from peak-gene pairs with raw P < 0.05.  The all-tested-peak annotation is
# included as a conditioning covariate in addition to baseline-LD v2.2.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
raw_p_cutoff <- if (length(args) >= 2L) as.numeric(args[[2]]) else 0.05

source(file.path(root, "scripts", "GWAS", "sLDSC", "source_results_paths.R"))
paths <- sldsc_source_results_paths(root)
for (d in c(paths$celltype_annotations, paths$pooled_annotations,
            paths$celltype_backgrounds, paths$pooled_backgrounds)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

scent_root <- file.path(root, "results", "validation", "GWAS", "SCENT_style", "UKBB", "source_results")
scent_manifest <- fread(file.path(scent_root, "celltype", "manifest.tsv"))

method_order <- c("epilink_poisson", "epilink_nb", "scmultimap",
                  "signac_corrected", "signac_uncorrected")
method_labels <- c(
  epilink_poisson = "EpiLink-Poisson",
  epilink_nb = "EpiLink-NB",
  scmultimap = "scMultiMap",
  signac_corrected = "Signac (corrected)",
  signac_uncorrected = "Signac (uncorrected)"
)

keep_celltype <- function(dataset, celltype) {
  if (dataset %chin% c("pbmc", "bmmc")) return(celltype %chin% c("T", "Bcell", "Myeloid"))
  if (dataset %chin% c("Dogma_CTRL", "Dogma_STIM")) return(celltype == "T")
  if (dataset == "SHARE") return(celltype == "LCL")
  if (dataset == "Arthritis") return(celltype %chin% c("T", "B", "Myeloid", "Endothelial", "Fibroblast"))
  if (dataset == "Brain") return(celltype %chin% c("Excitatory neurons", "Astrocytes", "Oligodendroglia"))
  if (dataset %chin% c("Brain2_ALL", "Brain2_AD")) return(celltype %chin% c("Astrocytes", "Microglia"))
  if (dataset == "Pituitary") return(celltype %chin% c("Corticotropes", "Gonadotropes", "Lactotropes",
                                                        "Somatotropes", "Stem cells", "Thyrotropes"))
  FALSE
}

scope <- unique(scent_manifest[method %chin% method_order,
  .(dataset, celltype, method, source_build, result_dir, analysis_variant,
    annotation_source)])
scope <- scope[mapply(keep_celltype, dataset, celltype)]
scope[, method_order_index := match(method, method_order)]
setorder(scope, dataset, celltype, method_order_index)

safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

panel_for <- function(dataset, celltype = NULL, pooled = FALSE) {
  if (dataset %chin% c("Brain", "Brain2_ALL", "Brain2_AD")) return("brain")
  if (dataset == "Pituitary") return("pituitary")
  if (dataset == "Arthritis" && (pooled || celltype %chin% c("Endothelial", "Fibroblast"))) return("arthritis_context")
  "immune"
}

interpretation_for <- function(dataset, celltype = NULL) {
  if (dataset == "Pituitary") return("exploratory")
  if (dataset == "SHARE") return("exploratory")
  if (dataset == "Arthritis" && celltype %chin% c("Endothelial", "Fibroblast")) return("exploratory")
  "primary"
}

result_base <- function(dataset, celltype) {
  if (grepl("^Brain2_", dataset)) {
    condition <- sub("^Brain2_", "", dataset)
    file.path(root, "results", "real", "Brain2", celltype, condition)
  } else {
    file.path(root, "results", "real", dataset, celltype)
  }
}

method_dir <- function(row) {
  base <- result_base(row$dataset, row$celltype)
  if (row$method %chin% c("signac_corrected", "signac_uncorrected")) return(file.path(base, "signac"))
  if (row$method == "scmultimap") return(file.path(base, "scmultimap"))
  rd <- row$result_dir
  if (!is.na(rd) && nzchar(rd) && dir.exists(rd)) return(rd)
  file.path(base, row$method)
}

p_column <- function(method, columns) {
  wanted <- switch(method,
    epilink_poisson = c("pval"),
    epilink_nb = c("pval"),
    scmultimap = c("pval"),
    # Older Signac exports call the corrected two-sided value `pvalue`;
    # newer exports use the explicit `pvalue_corrected` name.
    signac_corrected = c("pvalue_corrected", "pvalue"),
    signac_uncorrected = c("pvalue_uncorrected")
  )
  hit <- wanted[wanted %chin% columns]
  if (!length(hit)) NA_character_ else hit[[1L]]
}

parse_peaks <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(data.table(chr=character(), start=integer(), end=integer()))
  z <- tstrsplit(x, "-", fixed = TRUE)
  if (length(z) != 3L) stop("Peak IDs must have chr-start-end form")
  ans <- data.table(chr=z[[1]], start=as.integer(z[[2]]), end=as.integer(z[[3]]))
  ans[grepl("^chr([1-9]|1[0-9]|2[0-2]|X)$", chr) & !is.na(start) & !is.na(end) & start <= end]
}

reduce_bed <- function(x) {
  if (!nrow(x)) return(data.table(chr=character(), start=integer(), end=integer()))
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand=TRUE)
  data.table(chr=as.character(seqnames(gr)), start=start(gr), end=end(gr))
}

chain_path <- file.path(paths$base, "reference", "raw", "hg38ToHg19.over.chain.gz")
tmp_chain <- tempfile(fileext = ".chain")
writeLines(readLines(gzfile(chain_path), warn=FALSE), tmp_chain)
chain38to19 <- import.chain(tmp_chain)

to_hg19 <- function(x, build) {
  x <- reduce_bed(x)
  if (!nrow(x) || identical(build, "hg19")) return(x)
  lifted <- unlist(liftOver(GRanges(x$chr, IRanges(x$start, x$end)), chain38to19), use.names=FALSE)
  if (!length(lifted)) return(data.table(chr=character(), start=integer(), end=integer()))
  reduce_bed(data.table(chr=as.character(seqnames(lifted)), start=start(lifted), end=end(lifted)))
}

write_bed <- function(x, path) {
  y <- copy(x)
  y[, start0 := pmax(start - 1L, 0L)]
  fwrite(y[, .(chr, start0, end)], path, sep="\t", col.names=FALSE)
}

celltype_rows <- list()
audit_rows <- list()
background_rows <- list()
k <- 0L; ka <- 0L; kb <- 0L

groups <- unique(scope[, .(dataset, celltype)])
for (gi in seq_len(nrow(groups))) {
  ds <- groups$dataset[gi]
  ct <- groups$celltype[gi]
  z <- scope[dataset == ds & celltype == ct]

  candidate_peaks <- character()
  candidate_paths <- character()
  for (i in seq_len(nrow(z))) {
    d <- method_dir(z[i])
    cf <- file.path(d, "candidate_pairs.rds")
    if (file.exists(cf)) {
      cp <- tryCatch(readRDS(cf), error=function(e) NULL)
      if (!is.null(cp) && "peak" %in% names(cp)) {
        candidate_peaks <- union(candidate_peaks, as.character(cp$peak))
        candidate_paths <- union(candidate_paths, cf)
      }
    }
  }
  builds <- unique(na.omit(z$source_build))
  build <- if (length(builds)) builds[[1L]] else "hg38"
  bg_id <- paste(safe_id(ds), safe_id(ct), "tested_peaks", sep="__")
  bg_hg19 <- to_hg19(parse_peaks(candidate_peaks), build)
  write_bed(bg_hg19, file.path(paths$celltype_backgrounds, paste0(bg_id, ".hg19.bed")))
  kb <- kb + 1L
  background_rows[[kb]] <- data.table(
    background_id=bg_id, dataset=ds, celltype=ct, panel=panel_for(ds, ct),
    source_build=build, n_candidate_peak_ids=uniqueN(candidate_peaks),
    n_hg19_intervals=nrow(bg_hg19), candidate_sources=paste(candidate_paths, collapse=";")
  )

  for (i in seq_len(nrow(z))) {
    a <- z[i]
    d <- method_dir(a)
    all_file <- file.path(d, "all_pairs.tsv.gz")
    target_id <- paste(safe_id(ds), safe_id(ct), a$method, "rawp05", sep="__")
    status <- "ok"; reason <- ""
    n_tested <- NA_integer_; n_sig <- NA_integer_; sig_peaks <- character()
    if (!file.exists(all_file)) {
      status <- "missing"; reason <- "all_pairs.tsv.gz not found"
    } else {
      cols <- names(fread(all_file, nrows=0L))
      pc <- p_column(a$method, cols)
      if (is.na(pc) || !"peak" %chin% cols) {
        status <- "missing"; reason <- "required peak/P-value column absent"
      } else {
        dt <- fread(all_file, select=c("peak", pc), showProgress=FALSE)
        setnames(dt, pc, "raw_p")
        n_tested <- nrow(dt)
        sig <- dt[is.finite(raw_p) & raw_p < raw_p_cutoff]
        n_sig <- nrow(sig)
        sig_peaks <- unique(as.character(sig$peak))
      }
    }
    target_hg19 <- to_hg19(parse_peaks(sig_peaks), build)
    nested_fraction <- if (length(sig_peaks)) mean(sig_peaks %chin% candidate_peaks) else NA_real_
    if (status == "ok" && (!nrow(target_hg19) || is.na(nested_fraction) || nested_fraction < 0.999)) {
      status <- "invalid"
      reason <- if (!nrow(target_hg19)) "no linked peaks at raw P cutoff" else "linked peaks not nested in tested-peak universe"
    }
    if (status == "ok") {
      write_bed(target_hg19, file.path(paths$celltype_annotations, paste0(target_id, ".hg19.bed")))
      k <- k + 1L
      celltype_rows[[k]] <- data.table(
        annotation_id=target_id, background_id=bg_id, dataset=ds, celltype=ct,
        panel=panel_for(ds, ct), interpretation=interpretation_for(ds, ct),
        method=a$method, method_label=unname(method_labels[a$method]),
        raw_p_cutoff=raw_p_cutoff, source_build=build, source_file=all_file,
        n_tested_pairs=n_tested, n_significant_pairs=n_sig,
        n_linked_peak_ids=uniqueN(sig_peaks), n_hg19_intervals=nrow(target_hg19),
        n_background_hg19_intervals=nrow(bg_hg19), nested_fraction=nested_fraction,
        annotation_source=a$annotation_source, analysis_variant=a$analysis_variant
      )
    }
    ka <- ka + 1L
    audit_rows[[ka]] <- data.table(
      view="celltype", annotation_id=target_id, background_id=bg_id,
      dataset=ds, celltype=ct, method=a$method, status=status, reason=reason,
      source_file=all_file, n_tested_pairs=n_tested, n_significant_pairs=n_sig,
      n_linked_peak_ids=uniqueN(sig_peaks), n_hg19_intervals=nrow(target_hg19),
      nested_fraction=nested_fraction
    )
  }
}

celltype_manifest <- rbindlist(celltype_rows, fill=TRUE)
celltype_backgrounds <- rbindlist(background_rows, fill=TRUE)
celltype_manifest[, method_order_index := match(method, method_order)]
setorder(celltype_manifest, dataset, celltype, method_order_index)
celltype_manifest[, method_order_index := NULL]
fwrite(celltype_manifest, file.path(paths$celltype, "manifest.tsv"), sep="\t")
fwrite(celltype_backgrounds, file.path(paths$celltype, "background_manifest.tsv"), sep="\t")

# Dataset-level unions, matching the pooled view used by the SCENT-style run.
pooled_rows <- list(); pooled_background_rows <- list(); kp <- 0L; kpb <- 0L
for (ds in unique(celltype_manifest$dataset)) {
  dct <- celltype_manifest[dataset == ds]
  bg_ids <- unique(dct$background_id)
  bg <- rbindlist(lapply(bg_ids, function(id) {
    f <- file.path(paths$celltype_backgrounds, paste0(id, ".hg19.bed"))
    if (!file.exists(f) || !file.info(f)$size) return(data.table(chr=character(), start=integer(), end=integer()))
    fread(f, header=FALSE, col.names=c("chr", "start0", "end"))[, .(chr, start=start0+1L, end)]
  }), fill=TRUE)
  bg <- reduce_bed(bg)
  pooled_bg_id <- paste(safe_id(ds), "ALL", "tested_peaks", sep="__")
  write_bed(bg, file.path(paths$pooled_backgrounds, paste0(pooled_bg_id, ".hg19.bed")))
  kpb <- kpb + 1L
  pooled_background_rows[[kpb]] <- data.table(
    background_id=pooled_bg_id, dataset=ds, celltype="ALL_REQUESTED",
    panel=panel_for(ds, pooled=TRUE), n_source_backgrounds=length(bg_ids),
    n_hg19_intervals=nrow(bg)
  )
  for (m in method_order) {
    z <- dct[method == m]
    if (!nrow(z)) next
    beds <- rbindlist(lapply(z$annotation_id, function(id) {
      f <- file.path(paths$celltype_annotations, paste0(id, ".hg19.bed"))
      if (!file.exists(f) || !file.info(f)$size) return(data.table(chr=character(), start=integer(), end=integer()))
      fread(f, header=FALSE, col.names=c("chr", "start0", "end"))[, .(chr, start=start0+1L, end)]
    }), fill=TRUE)
    beds <- reduce_bed(beds)
    target_id <- paste(safe_id(ds), "ALL", m, "rawp05", sep="__")
    write_bed(beds, file.path(paths$pooled_annotations, paste0(target_id, ".hg19.bed")))
    expected <- uniqueN(scope[dataset == ds, celltype])
    kp <- kp + 1L
    pooled_rows[[kp]] <- data.table(
      annotation_id=target_id, background_id=pooled_bg_id, dataset=ds,
      celltype="ALL_REQUESTED", panel=panel_for(ds, pooled=TRUE),
      interpretation=if (ds %chin% c("Pituitary", "SHARE")) "exploratory" else "primary",
      method=m, method_label=unname(method_labels[m]), raw_p_cutoff=raw_p_cutoff,
      source_build="mixed", n_source_celltypes=nrow(z), n_expected_celltypes=expected,
      complete_celltype_coverage=nrow(z) == expected,
      n_significant_pairs=sum(z$n_significant_pairs, na.rm=TRUE),
      n_linked_peak_ids=NA_integer_, n_hg19_intervals=nrow(beds),
      n_background_hg19_intervals=nrow(bg), nested_fraction=NA_real_
    )
    ka <- ka + 1L
    audit_rows[[ka]] <- data.table(
      view="pooled", annotation_id=target_id, background_id=pooled_bg_id,
      dataset=ds, celltype="ALL_REQUESTED", method=m,
      status=if (nrow(z) == expected) "ok" else "partial",
      reason=if (nrow(z) == expected) "" else paste0("available cell types ", nrow(z), "/", expected),
      n_significant_pairs=sum(z$n_significant_pairs, na.rm=TRUE),
      n_hg19_intervals=nrow(beds)
    )
  }
}

pooled_manifest <- rbindlist(pooled_rows, fill=TRUE)
pooled_backgrounds <- rbindlist(pooled_background_rows, fill=TRUE)
pooled_manifest[, method_order_index := match(method, method_order)]
setorder(pooled_manifest, dataset, method_order_index)
pooled_manifest[, method_order_index := NULL]
fwrite(pooled_manifest, file.path(paths$pooled, "manifest.tsv"), sep="\t")
fwrite(pooled_backgrounds, file.path(paths$pooled, "background_manifest.tsv"), sep="\t")
fwrite(rbindlist(audit_rows, fill=TRUE), file.path(paths$root, "annotation_audit.tsv"), sep="\t")

message("Celltype targets: ", nrow(celltype_manifest), "; backgrounds: ", nrow(celltype_backgrounds))
message("Pooled targets: ", nrow(pooled_manifest), "; backgrounds: ", nrow(pooled_backgrounds))
message("Audit: ", file.path(paths$root, "annotation_audit.tsv"))
