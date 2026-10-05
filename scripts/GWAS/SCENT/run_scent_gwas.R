#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
n_boot <- if (length(args) >= 2L) as.integer(args[[2]]) else 1000L
seed <- if (length(args) >= 3L) as.integer(args[[3]]) else 20260925L
set.seed(seed)

gwas_dir <- file.path(root, "results", "validation", "GWAS")
manifest_path <- Sys.getenv("SCENT_MANIFEST", file.path(gwas_dir, "SCENT_style", "config", "annotation_manifest.tsv"))
annotation_dir <- Sys.getenv("SCENT_ANNOTATION_DIR", file.path(gwas_dir, "SCENT_style", "testing", "annotations"))
out_dir <- Sys.getenv("SCENT_OUT_DIR", file.path(gwas_dir, "SCENT_style", "testing", "results"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(manifest_path)
traits <- fread(file.path(gwas_dir, "SCENT_style", "config", "traits.tsv"))
fm <- readRDS(file.path(gwas_dir, "reference", "FinnGen", "finngen_R12_susie_all.rds"))
common <- readRDS(file.path(gwas_dir, "reference", "1000G", "1000G_EUR_MAFgt0.01.hg38.rds"))
fm <- fm[pip >= 0.2]
common[, key := paste(chr, pos, sep = ":")]

read_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(data.table(chr = character(), start = integer(), end = integer()))
  x <- fread(path, header = FALSE, col.names = c("chr", "start0", "end"))
  x[, .(chr, start = start0 + 1L, end)]
}

reduce_regions <- function(x) {
  if (!nrow(x)) return(x)
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand = TRUE)
  data.table(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
}

parse_regions <- function(x) {
  y <- tstrsplit(unique(x), "[:-]")
  reduce_regions(data.table(chr = y[[1]], start = as.integer(y[[2]]), end = as.integer(y[[3]])))
}

overlap_mask <- function(points, intervals) {
  hit <- rep(FALSE, nrow(points))
  if (!nrow(points) || !nrow(intervals)) return(hit)
  rows_by_chr <- split(seq_len(nrow(points)), points$chr)
  for (cc in intersect(names(rows_by_chr), unique(intervals$chr))) {
    rows <- rows_by_chr[[cc]]
    iv <- intervals[chr == cc][order(start)]
    idx <- findInterval(points$pos[rows], iv$start)
    ok <- idx > 0L
    local_hit <- rep(FALSE, length(rows))
    local_hit[ok] <- points$pos[rows][ok] <= iv$end[idx[ok]]
    hit[rows] <- local_hit
  }
  hit
}

subset_points <- function(points, intervals) points[overlap_mask(points, intervals)]
count_overlap <- function(points, intervals) sum(overlap_mask(points, intervals))

trait_cache <- list()
for (tr in unique(traits$trait)) {
  all_trait <- readRDS(file.path(gwas_dir, "reference", "FinnGen", "finngen_R12_susie_all.rds"))[trait == tr]
  regions <- parse_regions(all_trait$region)
  bg <- subset_points(common, regions)
  ca <- fm[trait == tr, .(chr, pos = position, pip, key = paste(chr, position, sep = ":"))]
  ca <- ca[!duplicated(key)]
  bg_nonca <- bg[!key %chin% ca$key]
  trait_cache[[tr]] <- list(regions = regions, causal = ca, background = bg, background_nonca = bg_nonca)
}

rows <- list()
k <- 0L
for (i in seq_len(nrow(manifest))) {
  a <- manifest[i]
  ann_path <- file.path(annotation_dir, paste0(a$annotation_id, ".hg38.bed"))
  ann <- read_bed(ann_path)
  panel_traits <- unique(traits[panel == a$panel, trait])
  if (!length(panel_traits)) next
  annotation_bp <- if (nrow(ann)) sum(ann$end - ann$start + 1) else 0
  for (tr in panel_traits) {
    z <- trait_cache[[tr]]
    n_ca <- nrow(z$causal)
    n_bg <- nrow(z$background)
    n_ca_ann <- count_overlap(z$causal, ann)
    n_bg_ann <- count_overlap(z$background, ann)
    n_bg_nonca <- nrow(z$background_nonca)
    n_bg_ann_nonca <- count_overlap(z$background_nonca, ann)
    enrichment <- if (n_ca > 0 && n_bg > 0 && n_bg_ann > 0) (n_ca_ann / n_bg_ann) / (n_ca / n_bg) else NA_real_
    mat <- matrix(c(n_ca_ann, n_ca - n_ca_ann,
                    n_bg_ann_nonca, n_bg_nonca - n_bg_ann_nonca), nrow = 2L, byrow = TRUE)
    ft <- if (all(dim(mat) == c(2L, 2L)) && sum(mat) > 0) fisher.test(mat, alternative = "greater") else NULL
    k <- k + 1L
    rows[[k]] <- data.table(
      annotation_id = a$annotation_id, dataset = a$dataset, celltype = a$celltype,
      panel = a$panel, interpretation = a$interpretation,
      method = a$method, method_label = a$method_label,
      trait = tr, trait_label = traits[trait == tr, label][1],
      n_peak_intervals = nrow(ann), annotation_bp = annotation_bp,
      n_causal_total = n_ca, n_causal_in_annotation = n_ca_ann,
      n_common_locus_total = n_bg, n_common_in_annotation = n_bg_ann,
      enrichment = enrichment,
      fisher_or = if (is.null(ft)) NA_real_ else unname(ft$estimate),
      fisher_p = if (is.null(ft)) NA_real_ else ft$p.value
    )
  }
}

per_trait <- rbindlist(rows, fill = TRUE)
per_trait[, fisher_fdr := p.adjust(fisher_p, method = "BH")]
fwrite(per_trait, file.path(out_dir, "scent_per_trait.tsv"), sep = "\t")

summaries <- per_trait[, {
  vals <- enrichment[is.finite(enrichment)]
  boots <- if (length(vals) >= 2L) replicate(n_boot, mean(sample(vals, length(vals), replace = TRUE))) else NA_real_
  list(
    n_traits = length(vals),
    mean_enrichment = if (length(vals)) mean(vals) else NA_real_,
    median_enrichment = if (length(vals)) median(vals) else NA_real_,
    ci_low = if (length(vals) >= 2L) unname(quantile(boots, 0.025, na.rm = TRUE)) else NA_real_,
    ci_high = if (length(vals) >= 2L) unname(quantile(boots, 0.975, na.rm = TRUE)) else NA_real_,
    pooled_causal_in = sum(n_causal_in_annotation),
    pooled_causal_total = sum(n_causal_total),
    pooled_common_in = sum(n_common_in_annotation),
    pooled_common_total = sum(n_common_locus_total),
    pooled_enrichment = (sum(n_causal_in_annotation) / sum(n_common_in_annotation)) /
      (sum(n_causal_total) / sum(n_common_locus_total))
  )
}, by = .(annotation_id, dataset, celltype, panel, interpretation, method, method_label)]

summaries[, `:=`(bootstrap_replicates = n_boot, bootstrap_unit = fifelse(n_traits >= 2, "trait", "not_available"), seed = seed)]
fwrite(summaries, file.path(out_dir, "scent_panel_summary.tsv"), sep = "\t")

audit <- per_trait[, .(
  n_annotations = uniqueN(annotation_id), n_traits = uniqueN(trait),
  min_causal_variants = min(n_causal_total), max_causal_variants = max(n_causal_total),
  zero_peak_annotations = sum(n_peak_intervals == 0L),
  zero_causal_hits = sum(n_causal_in_annotation == 0L)
), by = panel]
fwrite(audit, file.path(out_dir, "scent_audit.tsv"), sep = "\t")
message("Wrote SCENT-style results for ", uniqueN(per_trait$annotation_id), " annotations and ", uniqueN(per_trait$trait), " traits.")
