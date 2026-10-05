#!/usr/bin/env Rscript

# Sakaue et al.-style SCENT enrichment. Each annotation is one dataset-level
# union of all available EpiLink cell types. FinnGen and UKBB are evaluated
# separately. The final score is the arithmetic mean across traits and the
# 95% CI is a 1,000-replicate bootstrap over traits.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
n_boot <- if (length(args) >= 2L) as.integer(args[[2]]) else 1000L
seed <- if (length(args) >= 3L) as.integer(args[[3]]) else 20260930L
mode <- if (length(args) >= 4L) args[[4]] else "epilink"
source_arg <- if (length(args) >= 5L) args[[5]] else "both"
view <- if (length(args) >= 6L) args[[6]] else "legacy"
set.seed(seed)

gwas <- file.path(root, "results", "validation", "GWAS")
sakaue_root <- file.path(gwas, "SCENT_style", "UKBB", "sakaue_style")
source(file.path(root, "scripts", "GWAS", "SCENT_style", "source_results_paths.R"))
source_results <- source_results_paths(root)
base <- switch(view,
  pooled = source_results$pooled,
  celltype = source_results$celltype,
  pooled_main = source_results$pooled_main,
  celltype_main = source_results$celltype_main,
  pooled_brain2 = source_results$pooled_brain2,
  celltype_brain2 = source_results$celltype_brain2,
  if (identical(mode, "epilink")) sakaue_root else file.path(sakaue_root, mode)
)
out_dir <- file.path(base, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest_file <- if (view == "legacy") "sakaue_annotation_manifest.tsv" else "manifest.tsv"
manifest <- fread(file.path(base, manifest_file))

annotation_dirs <- switch(view,
  pooled = source_results$pooled_annotations,
  celltype = source_results$celltype_annotations,
  c(file.path(base, "annotations"))
)
annotation_path <- function(annotation_id) {
  candidates <- file.path(annotation_dirs, paste0(annotation_id, ".hg38.bed"))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Missing annotation BED for ", annotation_id)
  hit[[1L]]
}
common <- readRDS(file.path(gwas, "reference", "1000G", "1000G_EUR_MAFgt0.01.hg38.rds"))
common <- unique(common[, .(chr = as.character(chr), pos = as.integer(pos), SNP, MAF)], by = c("chr", "pos"))
common[, key := paste(chr, pos, sep = ":")]
common_gr <- GRanges(common$chr, IRanges(common$pos, common$pos))

fg <- readRDS(file.path(gwas, "reference", "FinnGen", "finngen_sakaue_reference.rds"))
uk_causal <- fread(file.path(gwas, "reference", "UKBB", "ukbb_94traits_causal.hg38.tsv.gz"))
uk_regions <- fread(file.path(gwas, "reference", "UKBB", "ukbb_94traits_regions.hg38.tsv.gz"))
uk_regions_list <- split(uk_regions[, .(chr, start, end)], uk_regions$trait)

fg_labels <- fread(file.path(gwas, "SCENT_style", "config", "traits.tsv"))[, .(trait, label)]
uk_labels <- fread(file.path(gwas, "reference", "UKBB", "UKBB_94traits_release1.traits"),
                  select = c("trait", "description"))
setnames(uk_labels, "description", "label")

read_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.table(chr = character(), start = integer(), end = integer()))
  }
  x <- fread(path, header = FALSE, col.names = c("chr", "start0", "end"))
  x[, .(chr, start = as.integer(start0) + 1L, end = as.integer(end))]
}

overlap_mask <- function(points, intervals) {
  hit <- rep(FALSE, nrow(points))
  if (!nrow(points) || !nrow(intervals)) return(hit)
  by_chr <- split(seq_len(nrow(points)), points$chr)
  for (cc in intersect(names(by_chr), unique(intervals$chr))) {
    rows <- by_chr[[cc]]
    iv <- intervals[chr == cc][order(start)]
    idx <- findInterval(points$pos[rows], iv$start)
    ok <- idx > 0L
    local <- rep(FALSE, length(rows))
    local[ok] <- points$pos[rows][ok] <= iv$end[idx[ok]]
    hit[rows] <- local
  }
  hit
}
count_overlap <- function(points, intervals) sum(overlap_mask(points, intervals))

common_in_regions <- function(regions) {
  if (!nrow(regions)) return(common[0])
  query <- GRanges(regions$chr, IRanges(regions$start, regions$end))
  hits <- findOverlaps(common_gr, query, ignore.strand = TRUE)
  if (!length(hits)) return(common[0])
  common[unique(queryHits(hits))]
}

reference_traits <- function(source) {
  if (source == "FinnGen") fg$trait_names else sort(unique(uk_causal$trait))
}

reference_trait_cache <- function(source) {
  traits <- reference_traits(source)
  out <- vector("list", length(traits))
  names(out) <- traits
  for (tr in traits) {
    if (source == "FinnGen") {
      regions <- fg$regions[[tr]]
      causal <- fg$causal[trait == tr, .(chr, pos)]
    } else {
      regions <- uk_regions_list[[tr]]
      causal <- uk_causal[trait == tr, .(chr, pos)]
    }
    if (is.null(regions)) regions <- data.table(chr = character(), start = integer(), end = integer())
    if (is.null(causal)) causal <- data.table(chr = character(), pos = integer())
    regions <- unique(regions)
    causal <- unique(causal)
    # Use the indexed GRanges overlap rather than scanning all 8.5M common
    # variants with a per-trait R loop. The statistic is identical.
    bg <- common_in_regions(regions)
    causal[, key := paste(chr, pos, sep = ":")]
    bg[, key := paste(chr, pos, sep = ":")]
    out[[tr]] <- list(regions = regions, causal = causal, background = bg,
                      background_noncausal = bg[!key %chin% causal$key])
    message(source, " ", tr, ": causal=", nrow(causal), "; common-locus=", nrow(bg))
  }
  out
}

sources_to_run <- if (source_arg %in% c("FinnGen", "UKBB")) source_arg else c("FinnGen", "UKBB")
cache <- setNames(lapply(sources_to_run, reference_trait_cache), sources_to_run)

# Combine all annotation intervals into one query set. For each trait we then
# do one overlap pass for causal variants and one for common variants, instead
# of rescanning the same points once per annotation.
annotation_tables <- lapply(seq_len(nrow(manifest)), function(i) {
  a <- manifest[i]
  z <- read_bed(annotation_path(a$annotation_id))
  z[, `:=`(annotation_id = a$annotation_id, annotation_bp = if (nrow(z)) sum(z$end - z$start + 1L) else 0)]
  z
})
annotation_intervals <- rbindlist(annotation_tables, fill = TRUE)
annotation_gr <- GRanges(annotation_intervals$chr,
                         IRanges(annotation_intervals$start, annotation_intervals$end),
                         annotation_id = annotation_intervals$annotation_id)

counts_by_annotation <- function(points) {
  if (!nrow(points)) return(setNames(integer(nrow(manifest)), manifest$annotation_id))
  point_gr <- GRanges(points$chr, IRanges(points$pos, points$pos))
  hits <- findOverlaps(point_gr, annotation_gr, ignore.strand = TRUE)
  ans <- setNames(integer(nrow(manifest)), manifest$annotation_id)
  if (length(hits)) {
    tab <- table(mcols(annotation_gr)$annotation_id[subjectHits(hits)])
    ans[names(tab)] <- as.integer(tab)
  }
  ans
}

rows <- list()
k <- 0L
for (source in sources_to_run) {
  for (tr in reference_traits(source)) {
    z <- cache[[source]][[tr]]
    n_ca <- nrow(z$causal)
    n_bg <- nrow(z$background)
    n_bg_nonca <- nrow(z$background_noncausal)
    causal_counts <- counts_by_annotation(z$causal)
    bg_counts <- counts_by_annotation(z$background)
    bg_nonca_counts <- counts_by_annotation(z$background_noncausal)
    for (i in seq_len(nrow(manifest))) {
      a <- manifest[i]
      id <- a$annotation_id
      n_ca_ann <- causal_counts[[id]]
      n_bg_ann <- bg_counts[[id]]
      n_bg_ann_nonca <- bg_nonca_counts[[id]]
      enrichment <- if (n_ca > 0L && n_bg > 0L && n_bg_ann > 0L)
        (n_ca_ann / n_bg_ann) / (n_ca / n_bg) else NA_real_
      tab <- matrix(c(n_ca_ann, n_ca - n_ca_ann,
                      n_bg_ann_nonca, n_bg_nonca - n_bg_ann_nonca), nrow = 2L, byrow = TRUE)
      ft <- if (sum(tab) > 0) fisher.test(tab, alternative = "greater") else NULL
      label <- if (source == "FinnGen") fg_labels[trait == tr, label][1] else uk_labels[trait == tr, label][1]
      k <- k + 1L
      rows[[k]] <- data.table(
        source = source, annotation_id = id, dataset = a$dataset,
        celltype_scope = a$celltype, n_celltypes = a$n_celltypes, celltypes = a$celltypes,
        method = a$method, method_label = a$method_label, correction = a$correction,
        fdr = a$fdr, trait = tr, trait_label = label,
        n_peak_intervals = sum(annotation_intervals$annotation_id == id),
        annotation_bp = annotation_intervals[annotation_id == id, annotation_bp[1]],
        n_causal_total = n_ca, n_causal_in_annotation = n_ca_ann,
        n_common_locus_total = n_bg, n_common_in_annotation = n_bg_ann,
        enrichment = enrichment, fisher_or = if (is.null(ft)) NA_real_ else unname(ft$estimate),
        fisher_p = if (is.null(ft)) NA_real_ else ft$p.value
      )
    }
    message("Finished trait ", source, " / ", tr)
  }
}

per_trait <- rbindlist(rows, fill = TRUE)
per_trait[, fisher_fdr := p.adjust(fisher_p, method = "BH")]
per_trait_file <- switch(view,
  pooled_brain2 = "brain2_pituitary_pooled_per_trait.tsv.gz",
  celltype_brain2 = "brain2_pituitary_celltype_per_trait.tsv.gz",
  "scent_sakaue_per_trait.tsv.gz"
)
summary_file <- switch(view,
  pooled_brain2 = "brain2_pituitary_pooled_summary.tsv",
  celltype_brain2 = "brain2_pituitary_celltype_summary.tsv",
  "scent_sakaue_summary.tsv"
)
fwrite(per_trait, file.path(out_dir, per_trait_file), sep = "\t")

summary_rows <- per_trait[, {
  vals <- enrichment[is.finite(enrichment)]
  boot <- if (length(vals) >= 2L) replicate(n_boot, mean(sample(vals, length(vals), replace = TRUE))) else NA_real_
  list(
    n_traits = length(vals), mean_enrichment = if (length(vals)) mean(vals) else NA_real_,
    median_enrichment = if (length(vals)) median(vals) else NA_real_,
    ci_low = if (length(vals) >= 2L) unname(quantile(boot, 0.025, na.rm = TRUE)) else NA_real_,
    ci_high = if (length(vals) >= 2L) unname(quantile(boot, 0.975, na.rm = TRUE)) else NA_real_,
    pooled_causal_in = sum(n_causal_in_annotation), pooled_causal_total = sum(n_causal_total),
    pooled_common_in = sum(n_common_in_annotation), pooled_common_total = sum(n_common_locus_total),
    pooled_enrichment = if (sum(n_common_in_annotation) > 0 && sum(n_common_locus_total) > 0)
      (sum(n_causal_in_annotation) / sum(n_common_in_annotation)) /
        (sum(n_causal_total) / sum(n_common_locus_total)) else NA_real_
  )
}, by = .(source, annotation_id, dataset, celltype_scope, n_celltypes, celltypes,
          method, method_label, correction, fdr)]
summary_rows[, `:=`(bootstrap_replicates = n_boot, bootstrap_unit = "trait", seed = seed,
                    reference_definition = fifelse(source == "FinnGen", "FinnGen R12 project subset", "UKBB 94 traits release 1.1"))]
fwrite(summary_rows, file.path(out_dir, summary_file), sep = "\t")

audit <- per_trait[, .(
  n_annotations = uniqueN(annotation_id), n_traits = uniqueN(trait),
  zero_causal_hits = sum(n_causal_in_annotation == 0L),
  zero_annotation_hits = sum(n_common_in_annotation == 0L),
  min_common_locus = min(n_common_locus_total), max_common_locus = max(n_common_locus_total)
), by = source]
fwrite(audit, file.path(out_dir, "scent_sakaue_audit.tsv"), sep = "\t")
message("Wrote ", nrow(per_trait), " per-trait rows and ", nrow(summary_rows), " summaries")
