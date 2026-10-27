# Helpers for the pgBoost-style GRCh38 K562 CRISPR SNP-gene benchmark.

crispr_new_method_specs <- function() {
  data.frame(
    method = c("Calibre-Poisson", "Calibre-NB", "Poisson-GLM", "NB-GLM",
               "scMultiMap", "scMultiMap-noBatch", "Signac (corrected)",
               "Signac (uncorrected)"),
    method_id = c("calibre_poisson", "calibre_nb", "poisson_glm", "nb_glm",
                  "scmultimap", "scmultimap_nobatch", "signac",
                  "signac_uncorrected"),
    stringsAsFactors = FALSE
  )
}

read_crispr_new_benchmark <- function(path) {
  benchmark <- data.table::fread(path, data.table = FALSE,
                                 showProgress = FALSE)
  required <- c("snp", "gene", "chr", "snp_start", "snp_end",
                "positive", "distance_to_tss")
  if (!all(required %in% names(benchmark))) {
    stop("CRISPR_new benchmark lacks: ",
         paste(setdiff(required, names(benchmark)), collapse = ", "),
         call. = FALSE)
  }
  benchmark$positive <- as.logical(benchmark$positive)
  benchmark
}

map_crispr_new_to_pairs <- function(pairs, benchmark) {
  pairs <- unique(pairs[, c("gene", "peak"), drop = FALSE])
  if (!nrow(pairs)) return(data.frame())
  peak_ranges <- links_to_ranges(pairs)
  snp_ranges <- GenomicRanges::GRanges(
    seqnames = benchmark$chr,
    ranges = IRanges::IRanges(
      start = as.integer(benchmark$snp_start) + 1L,
      end = as.integer(benchmark$snp_end)
    )
  )
  common <- intersect(GenomeInfoDb::seqlevelsInUse(peak_ranges),
                      GenomeInfoDb::seqlevelsInUse(snp_ranges))
  if (!length(common)) return(data.frame())
  peak_ranges <- GenomeInfoDb::keepSeqlevels(
    peak_ranges, common, pruning.mode = "coarse"
  )
  snp_ranges$benchmark_index <- seq_along(snp_ranges)
  snp_ranges <- GenomeInfoDb::keepSeqlevels(
    snp_ranges, common, pruning.mode = "coarse"
  )
  hits <- GenomicRanges::findOverlaps(snp_ranges, peak_ranges,
                                      ignore.strand = TRUE)
  if (!length(hits)) return(data.frame())
  benchmark_index <- snp_ranges$benchmark_index[S4Vectors::queryHits(hits)]
  peak_index <- S4Vectors::subjectHits(hits)
  keep <- benchmark$gene[benchmark_index] == peak_ranges$gene[peak_index]
  if (!any(keep)) return(data.frame())
  benchmark_index <- benchmark_index[keep]
  peak_index <- peak_index[keep]
  unique(data.frame(
    snp = benchmark$snp[benchmark_index],
    gene = benchmark$gene[benchmark_index],
    peak = sub("^[^\r]*\r", "", peak_ranges$pair_key[peak_index]),
    positive = benchmark$positive[benchmark_index],
    distance_to_tss = benchmark$distance_to_tss[benchmark_index],
    stringsAsFactors = FALSE
  ))
}

crispr_new_score_column <- function(results, method_id = NULL) {
  if (identical(method_id, "signac_uncorrected") &&
      "pvalue_uncorrected" %in% names(results)) {
    return("pvalue_uncorrected")
  }
  if (identical(method_id, "signac") &&
      "pvalue_corrected" %in% names(results)) {
    return("pvalue_corrected")
  }
  column <- intersect(c("pval", "pvalue"), names(results))
  if (!length(column)) {
    stop("Method result has neither pval nor pvalue", call. = FALSE)
  }
  column[[1L]]
}

score_crispr_new_links <- function(mapping, results, pvalue_column = NULL) {
  if (!nrow(mapping)) return(data.frame())
  if (is.null(pvalue_column)) pvalue_column <- crispr_new_score_column(results)
  result_key <- paste(results$gene, results$peak, sep = "\r")
  map_key <- paste(mapping$gene, mapping$peak, sep = "\r")
  pvalue <- as.numeric(results[[pvalue_column]][match(map_key, result_key)])
  score <- -log10(pmax(pvalue, .Machine$double.xmin))
  score[!is.finite(score)] <- NA_real_
  scored <- data.frame(mapping, score = score, stringsAsFactors = FALSE)
  scored <- data.table::as.data.table(scored)[, .(
    positive = any(positive),
    distance_to_tss = min(distance_to_tss),
    score = if (all(is.na(score))) NA_real_ else max(score, na.rm = TRUE)
  ), by = .(snp, gene)]
  as.data.frame(scored)
}

crispr_new_rank_metrics <- function(scored, distance_threshold,
                                    top_fraction = 0.05) {
  usable <- !is.na(scored$positive) & is.finite(scored$distance_to_tss) &
    scored$distance_to_tss > distance_threshold
  x <- scored[usable, , drop = FALSE]
  if (!nrow(x) || !any(x$positive) || !any(!x$positive)) return(NULL)
  finite <- x$score[is.finite(x$score)]
  floor_score <- if (length(finite)) min(finite) - 1 else 0
  x$score[!is.finite(x$score)] <- floor_score

  grouped <- data.table::as.data.table(x)[, .(
    n = .N, positives = sum(positive)
  ), by = score][order(-score)]
  grouped[, selected := cumsum(n)]
  grouped[, tp := cumsum(positives)]
  n_positive <- sum(x$positive)
  n_negative <- nrow(x) - n_positive
  grouped[, recall := tp / n_positive]
  grouped[, precision := tp / selected]
  prevalence <- n_positive / nrow(x)
  grouped[, enrichment := precision / prevalence]
  recall_increment <- diff(c(0, grouped$recall))
  auprc <- sum(recall_increment * grouped$precision)
  average_enrichment <- auprc / prevalence

  cutoff <- as.numeric(stats::quantile(
    x$score, probs = 1 - top_fraction, names = FALSE, type = 1
  ))
  selected <- x$score >= cutoff
  tp <- sum(selected & x$positive)
  fp <- sum(selected & !x$positive)
  fn <- sum(!selected & x$positive)
  tn <- sum(!selected & !x$positive)
  contingency <- matrix(c(tp, fp, fn, tn), nrow = 2L, byrow = TRUE)
  defined <- all(rowSums(contingency) > 0L) && all(colSums(contingency) > 0L)
  fisher <- if (defined) {
    stats::fisher.test(contingency, alternative = "greater")
  } else NULL

  list(
    metrics = data.frame(
      distance_threshold = distance_threshold,
      n_links = nrow(x), n_positive = n_positive, n_negative = n_negative,
      prevalence = prevalence, auprc = auprc,
      average_enrichment = average_enrichment,
      top_fraction = top_fraction, score_cutoff = cutoff,
      TP = tp, FP = fp, FN = fn, TN = tn,
      odds_ratio = if (defined) unname(fisher$estimate) else NA_real_,
      CI_low = if (defined) fisher$conf.int[[1L]] else NA_real_,
      CI_high = if (defined) fisher$conf.int[[2L]] else NA_real_,
      pval = if (defined) fisher$p.value else NA_real_
    ),
    curve = as.data.frame(grouped[, .(score, recall, precision, enrichment)])
  )
}

bootstrap_crispr_new_metrics <- function(scored, distance_threshold,
                                         n_bootstrap = 200L, seed = 1L) {
  if (n_bootstrap <= 0L) {
    return(c(average_enrichment_low = NA_real_,
             average_enrichment_high = NA_real_,
             auprc_low = NA_real_, auprc_high = NA_real_))
  }
  usable <- !is.na(scored$positive) & is.finite(scored$distance_to_tss) &
    scored$distance_to_tss > distance_threshold
  x <- scored[usable, , drop = FALSE]
  if (!nrow(x) || !any(x$positive) || !any(!x$positive)) {
    return(c(average_enrichment_low = NA_real_,
             average_enrichment_high = NA_real_,
             auprc_low = NA_real_, auprc_high = NA_real_))
  }
  set.seed(seed)
  boot <- replicate(n_bootstrap, {
    sampled <- x[sample.int(nrow(x), nrow(x), replace = TRUE), , drop = FALSE]
    value <- crispr_new_rank_metrics(sampled, distance_threshold = -1L)
    if (is.null(value)) c(NA_real_, NA_real_) else
      c(value$metrics$average_enrichment, value$metrics$auprc)
  })
  c(
    average_enrichment_low = stats::quantile(boot[1L, ], 0.025, na.rm = TRUE),
    average_enrichment_high = stats::quantile(boot[1L, ], 0.975, na.rm = TRUE),
    auprc_low = stats::quantile(boot[2L, ], 0.025, na.rm = TRUE),
    auprc_high = stats::quantile(boot[2L, ], 0.975, na.rm = TRUE)
  )
}
