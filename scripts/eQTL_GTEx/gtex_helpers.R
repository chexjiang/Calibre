suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})
source("R/eqtl_validation.R")

gtex_panels <- function(config) {
  rbindlist(lapply(names(config$benchmarks), function(dataset) {
    data.table(dataset = dataset, celltype = unlist(config$benchmarks[[dataset]]))
  }))
}

gtex_method_specs <- function() {
  data.table(
    method = c("epilink_poisson", "epilink_nb", "scmultimap",
               "signac_corrected", "signac_uncorrected"),
    source_method = c("epilink_poisson", "epilink_nb", "scmultimap", "signac", "signac"),
    significance_file = c("significant_pairs.tsv.gz", "significant_pairs.tsv.gz",
                          "significant_pairs.tsv.gz", "significant_pairs_corrected.tsv.gz",
                          "significant_pairs_uncorrected.tsv.gz")
  )
}

read_scope_universe <- function(dataset, config, celltype = NULL) {
  if (!is.null(celltype)) return(as.data.table(read_candidate_universe(dataset, celltype)))
  cells <- config$benchmarks[[dataset]]
  unique(rbindlist(lapply(cells, function(ct) {
    as.data.table(read_candidate_universe(dataset, ct))
  }), fill = TRUE))
}

read_scope_selected <- function(dataset, method_spec, config, celltype = NULL) {
  cells <- if (is.null(celltype)) config$benchmarks[[dataset]] else celltype
  z <- rbindlist(lapply(cells, function(ct) {
    x <- read_method_pairs(dataset, ct, method_spec$source_method, config,
                           significant = TRUE,
                           significance_file = method_spec$significance_file)
    if (is.null(x)) return(NULL)
    x[is.finite(score), .(gene, peak)]
  }), fill = TRUE)
  unique(z)
}

read_scope_scores <- function(dataset, method, config, celltype = NULL) {
  cells <- if (is.null(celltype)) config$benchmarks[[dataset]] else celltype
  unique(rbindlist(lapply(cells, function(ct) {
    read_method_score_table(dataset, ct, method, config)
  }), fill = TRUE))
}

variant_peak_overlap_gtex <- function(universe, variants, include_tissue = FALSE) {
  peaks <- parse_peak_table(universe$peak)
  peak_gr <- table_to_granges(peaks, "peak")
  metadata <- c("variant_id", "gene", "pip")
  if (include_tissue) metadata <- c(metadata, "tissue")
  var_gr <- table_to_granges(variants, metadata)
  common <- intersect(GenomeInfoDb::seqlevelsInUse(peak_gr),
                      GenomeInfoDb::seqlevelsInUse(var_gr))
  peak_gr <- GenomeInfoDb::keepSeqlevels(peak_gr, common, pruning.mode = "coarse")
  var_gr <- GenomeInfoDb::keepSeqlevels(var_gr, common, pruning.mode = "coarse")
  hits <- findOverlaps(peak_gr, var_gr, ignore.strand = TRUE)
  if (!length(hits)) return(data.table())
  out <- data.table(
    peak = peak_gr$peak[queryHits(hits)],
    variant_id = var_gr$variant_id[subjectHits(hits)],
    eqtl_gene = var_gr$gene[subjectHits(hits)],
    pip = var_gr$pip[subjectHits(hits)]
  )
  if (include_tissue) out[, tissue := var_gr$tissue[subjectHits(hits)]]
  unique(out)
}

background_for_eqtl <- function(path) {
  x <- as.data.table(readRDS(path))
  unique(x[, .(chr = as.character(chr), start = as.integer(pos),
               end = as.integer(pos), variant_id = paste(chr, pos, sep = ":"))])
}

scope_cache <- function(universe, aggregate_eqtl, tissue_eqtl, background,
                        cache_dir, tss, window_bp) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(cache_dir, c("background_peak_hits.rds", "aggregate_eqtl_peak_hits.rds",
                                  "tissue_eqtl_peak_hits.rds", "background_cis_counts.rds"))
  if (all(file.exists(paths))) {
    return(list(background_hits = readRDS(paths[[1L]]),
                aggregate_hits = readRDS(paths[[2L]]),
                tissue_hits = readRDS(paths[[3L]]),
                background_cis = readRDS(paths[[4L]])))
  }
  message("    caching 1000G and GTEx overlaps")
  bg_hits <- background_peak_overlap(universe, background)
  agg_hits <- variant_peak_overlap_gtex(universe, aggregate_eqtl, FALSE)
  tissue_hits <- variant_peak_overlap_gtex(universe, tissue_eqtl, TRUE)
  genes <- unique(universe$gene)
  gene_tss <- tss[gene %chin% genes]
  gene_tss <- gene_tss[!duplicated(gene)]
  bg_cis <- count_positions_in_windows(background, gene_tss, window_bp)
  names(bg_cis) <- gene_tss$gene
  saveRDS(bg_hits, paths[[1L]], compress = FALSE)
  saveRDS(agg_hits, paths[[2L]], compress = FALSE)
  saveRDS(tissue_hits, paths[[3L]], compress = FALSE)
  saveRDS(bg_cis, paths[[4L]], compress = FALSE)
  list(background_hits = bg_hits, aggregate_hits = agg_hits,
       tissue_hits = tissue_hits, background_cis = bg_cis)
}

selected_hit_counts <- function(overlaps, selected, tissue = FALSE) {
  empty <- if (tissue) {
    data.table(tissue = character(), gene = character(), causal_annotated = integer())
  } else {
    data.table(gene = character(), causal_annotated = integer())
  }
  if (!nrow(overlaps) || !nrow(selected)) return(empty)
  links <- unique(as.data.table(selected)[, .(gene, peak)])
  z <- merge(overlaps, links, by = "peak", allow.cartesian = TRUE)
  z <- z[gene == eqtl_gene]
  if (!nrow(z)) return(empty)
  if (tissue) {
    unique(z[, .(tissue, gene, variant_id)])[, .(causal_annotated = .N),
                                             by = .(tissue, gene)]
  } else {
    unique(z[, .(gene, variant_id)])[, .(causal_annotated = .N), by = gene]
  }
}

background_hit_counts <- function(background_hits, selected) {
  empty <- data.table(gene = character(), common_annotated = integer())
  if (!nrow(background_hits) || !nrow(selected)) return(empty)
  z <- merge(background_hits, unique(as.data.table(selected)[, .(gene, peak)]),
             by = "peak", allow.cartesian = TRUE)
  if (!nrow(z)) return(empty)
  unique(z[, .(gene, variant_id)])[, .(common_annotated = .N), by = gene]
}

assemble_counts <- function(total, hits, bg_hits, bg_cis, genes, by_tissue = FALSE) {
  if (by_tissue) {
    out <- merge(total, hits, by = c("tissue", "gene"), all.x = TRUE)
  } else {
    out <- merge(total, hits, by = "gene", all.x = TRUE)
  }
  out <- out[gene %chin% genes]
  out <- merge(out, bg_hits, by = "gene", all.x = TRUE)
  out[, common_total := as.integer(bg_cis[gene])]
  for (nm in c("causal_annotated", "common_annotated")) out[is.na(get(nm)), (nm) := 0L]
  out <- out[!is.na(common_total) & common_total > 0]
  out[, `:=`(
    causal_not_annotated = pmax(0L, causal_total - causal_annotated),
    common_not_annotated = pmax(0L, common_total - common_annotated)
  )]
  out[]
}

summarize_scent_counts <- function(counts, n_boot, seed) {
  pg <- per_gene_scent_enrichment(as.data.frame(counts))
  pooled <- odds_ratio_from_counts(as.data.frame(counts))
  boot <- bootstrap_scent_metrics(as.data.frame(counts), n_boot, seed)
  q <- function(x, p) { x <- x[is.finite(x)]; if (length(x)) unname(quantile(x, p)) else NA_real_ }
  med <- function(x) { x <- x[is.finite(x)]; if (length(x)) median(x) else NA_real_ }
  mn <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
  list(
    per_gene = as.data.table(pg), bootstrap = as.data.table(boot),
    summary = data.table(
      n_genes_all = nrow(counts), n_genes_scent = nrow(pg),
      pooled_or_point = pooled$odds_ratio_ha,
      pooled_or_fisher = pooled$odds_ratio_fisher, pooled_or_pval = pooled$pval,
      pooled_or_bootstrap_median = med(boot$pooled_or),
      pooled_or_ci_low = q(boot$pooled_or, .025), pooled_or_ci_high = q(boot$pooled_or, .975),
      scent_mean_point = mn(pg$scent_enrichment),
      scent_mean_bootstrap_median = med(boot$scent_mean_enrichment),
      scent_mean_ci_low = q(boot$scent_mean_enrichment, .025),
      scent_mean_ci_high = q(boot$scent_mean_enrichment, .975)
    )
  )
}
