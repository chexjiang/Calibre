# Shared helpers for SCENT-style and ranked eQTL validation.

read_eqtl_validation_config <- function(path = "config/eqtl_validation.yml") {
  yaml::read_yaml(path)
}

parse_peak_table <- function(peaks) {
  peaks <- unique(as.character(peaks))
  match <- regexec("^(.*?)[-:]([0-9]+)[-:]([0-9]+)$", peaks)
  fields <- regmatches(peaks, match)
  valid <- lengths(fields) == 4L
  if (!all(valid)) {
    stop("Cannot parse peak coordinates: ",
         paste(head(peaks[!valid], 3L), collapse = ", "), call. = FALSE)
  }
  data.table::data.table(
    peak = peaks,
    chr = vapply(fields, `[[`, character(1L), 2L),
    start = as.integer(vapply(fields, `[[`, character(1L), 3L)),
    end = as.integer(vapply(fields, `[[`, character(1L), 4L))
  )
}

table_to_granges <- function(x, metadata = character()) {
  gr <- GenomicRanges::GRanges(
    seqnames = x$chr,
    ranges = IRanges::IRanges(as.integer(x$start), as.integer(x$end))
  )
  for (column in metadata) S4Vectors::mcols(gr)[[column]] <- x[[column]]
  gr
}

normalize_seqnames <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

read_candidate_universe <- function(dataset, celltype,
                                    results_root = "results/real") {
  cell_dir <- file.path(results_root, dataset, celltype)
  candidates <- Sys.glob(file.path(cell_dir, "*", "candidate_pairs.rds"))
  if (length(candidates)) {
    x <- readRDS(candidates[[1L]])
    x <- as.data.frame(x[, c("gene", "peak")])
  } else {
    all_pairs <- Sys.glob(file.path(cell_dir, "*", "all_pairs.tsv.gz"))
    if (!length(all_pairs)) stop("No candidate pairs for ", dataset, "/", celltype)
    x <- data.table::fread(all_pairs[[1L]], select = c("gene", "peak"),
                           showProgress = FALSE)
  }
  unique(x[stats::complete.cases(x), c("gene", "peak")])
}

read_method_pairs <- function(dataset, celltype, method, config,
                              significant = FALSE,
                              results_root = "results/real",
                              score_type = c("primary", "evidence"),
                              source_method = NULL,
                              significance_file = NULL) {
  score_type <- match.arg(score_type)
  if (is.null(source_method)) source_method <- method
  method_config <- config$methods[[source_method]]
  if (is.null(method_config)) stop("Unknown method: ", method)
  file_name <- if (significant) {
    if (is.null(significance_file)) "significant_pairs.tsv.gz" else significance_file
  } else "all_pairs.tsv.gz"
  path <- file.path(results_root, dataset, celltype, source_method, file_name)
  if (!file.exists(path)) return(NULL)
  score_column <- if (score_type == "primary") {
    method_config$primary_score
  } else {
    method_config$evidence_score
  }
  available <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  columns <- c("gene", "peak", score_column)
  if (!all(columns %in% available)) {
    stop(path, " lacks columns: ", paste(setdiff(columns, available), collapse = ", "))
  }
  x <- data.table::fread(path, select = columns, showProgress = FALSE)
  data.table::setnames(x, score_column, "score")
  unique(x)
}

eqtl1_method_specs <- function(config) {
  specs <- data.table::data.table(
    method = names(config$methods),
    source_method = names(config$methods),
    significance_file = "significant_pairs.tsv.gz"
  )
  if ("signac" %in% names(config$methods)) {
    specs <- data.table::rbindlist(list(
      specs[method != "signac"],
      data.table::data.table(
        method = c("signac_corrected", "signac_uncorrected"),
        source_method = "signac",
        significance_file = c(
          "significant_pairs_corrected.tsv.gz",
          "significant_pairs_uncorrected.tsv.gz"
        )
      )
    ), fill = TRUE)
  }
  specs
}

read_method_score_table <- function(dataset, celltype, method, config,
                                    results_root = "results/real") {
  method_config <- config$methods[[method]]
  path <- file.path(results_root, dataset, celltype, method, "all_pairs.tsv.gz")
  if (!file.exists(path)) return(NULL)
  primary <- method_config$primary_score
  evidence <- method_config$evidence_score
  columns <- unique(c("gene", "peak", primary, evidence))
  x <- data.table::fread(path, select = columns, showProgress = FALSE)
  data.table::setnames(x, primary, "primary_score")
  if (evidence == primary) {
    x[, evidence_score := primary_score]
  } else {
    data.table::setnames(x, evidence, "evidence_score")
  }
  unique(x)
}

load_eqtl_reference <- function(group,
                                processed_dir = "reference/eQTL/shared/processed") {
  path <- file.path(processed_dir, paste0("OneK1K_", group, ".rds"))
  if (!file.exists(path)) stop("Missing processed reference: ", path,
                               "; run prepare_reference.R")
  readRDS(path)
}

load_eqtl_tss <- function(processed_dir = "reference/eQTL/shared/processed") {
  path <- file.path(processed_dir, "TSS.hg38.rds")
  if (!file.exists(path)) stop("Missing processed TSS reference: ", path)
  readRDS(path)
}

load_background_variants <- function(
    path = "reference/eQTL/shared/raw/pgboost/1000G_MAC5_SNPs.hg38.bed.gz",
    processed_path = "reference/eQTL/shared/processed/1000G_MAC5_SNPs.hg38.rds") {
  if (file.exists(processed_path)) return(readRDS(processed_path))
  x <- data.table::fread(path, header = FALSE, showProgress = FALSE)
  if (ncol(x) < 3L) stop("Background SNP BED requires at least three columns")
  x <- x[, .(chr = normalize_seqnames(V1), start = as.integer(V2) + 1L,
             end = as.integer(V3))]
  x[, variant_id := paste(chr, start, end, sep = ":")]
  x
}

filter_distal_eqtls <- function(eqtl, tss, minimum_distance = 1000L) {
  x <- merge(eqtl, tss[, c("gene", "chr", "tss")], by = c("gene", "chr"))
  x$distance_to_tss <- abs(x$start - x$tss)
  x[x$distance_to_tss > minimum_distance, , drop = FALSE]
}

eqtl_peak_overlap <- function(universe, eqtl) {
  peaks <- parse_peak_table(universe$peak)
  peak_gr <- table_to_granges(peaks, "peak")
  eqtl_gr <- table_to_granges(eqtl, c("variant_id", "gene", "pip"))
  common <- intersect(GenomeInfoDb::seqlevelsInUse(peak_gr),
                      GenomeInfoDb::seqlevelsInUse(eqtl_gr))
  peak_gr <- GenomeInfoDb::keepSeqlevels(peak_gr, common, pruning.mode = "coarse")
  eqtl_gr <- GenomeInfoDb::keepSeqlevels(eqtl_gr, common, pruning.mode = "coarse")
  hits <- GenomicRanges::findOverlaps(peak_gr, eqtl_gr, ignore.strand = TRUE)
  if (!length(hits)) return(data.table::data.table())
  data.table::data.table(
    peak = peak_gr$peak[S4Vectors::queryHits(hits)],
    variant_id = eqtl_gr$variant_id[S4Vectors::subjectHits(hits)],
    eqtl_gene = eqtl_gr$gene[S4Vectors::subjectHits(hits)],
    pip = eqtl_gr$pip[S4Vectors::subjectHits(hits)]
  )
}

make_ranked_validation_units <- function(universe, eqtl, pip_threshold = 0.5) {
  eqtl <- eqtl[eqtl$pip >= pip_threshold, , drop = FALSE]
  overlaps <- eqtl_peak_overlap(universe, eqtl)
  if (!nrow(overlaps)) return(data.table::data.table())
  pairs <- data.table::as.data.table(unique(universe[, c("gene", "peak")]))
  units <- merge(overlaps, pairs, by = "peak", allow.cartesian = TRUE)
  units[, positive := gene == eqtl_gene]
  unique(units[, .(variant_id, peak, gene, eqtl_gene, pip, positive)])
}

# EpiLink's significance calibration is two-sided: EpiLinkFilter() applies
# abs(beta) or abs(zval).  The eQTL validation labels also do not encode effect
# direction, so both directions must be eligible for ranking.  Keep the
# transform method-specific because the other methods expose different score
# semantics.
eqtl_rank_score <- function(score, method) {
  if (method %in% c("epilink_poisson", "epilink_nb",
                   "poisson_glm", "nb_glm")) {
    abs(score)
  } else {
    score
  }
}

ranked_curve <- function(scored_pairs, units, top_fractions) {
  scored <- data.table::as.data.table(scored_pairs)[is.finite(score)]
  scored <- scored[, .(score = max(score)), by = .(gene, peak)]
  x <- merge(units, scored, by = c("gene", "peak"))
  if (!nrow(x)) return(data.frame())
  # One prediction score is inherited by every variant inside its peak, as in
  # pgBoost SNP-gene candidate construction.
  data.table::setorder(x, -score, gene, peak, variant_id)
  n_positive <- sum(x$positive)
  base_rate <- mean(x$positive)
  out <- lapply(top_fractions, function(fraction) {
    n <- max(1L, min(nrow(x), ceiling(nrow(x) * fraction)))
    selected <- x[seq_len(n)]
    tp <- sum(selected$positive)
    data.frame(
      top_fraction = fraction,
      n_selected = n,
      TP = tp,
      FP = n - tp,
      recall = if (n_positive) tp / n_positive else NA_real_,
      precision = tp / n,
      baseline_precision = base_rate,
      enrichment = if (base_rate > 0) (tp / n) / base_rate else NA_real_,
      score_threshold = selected$score[[n]],
      n_units = nrow(x),
      n_positive = n_positive
    )
  })
  data.table::rbindlist(out)
}

count_positions_in_windows <- function(variants, genes, window_bp) {
  split_positions <- lapply(split(variants$start, variants$chr), sort)
  counts <- integer(nrow(genes))
  for (i in seq_len(nrow(genes))) {
    positions <- split_positions[[genes$chr[[i]]]]
    if (is.null(positions)) next
    lower <- genes$tss[[i]] - window_bp
    upper <- genes$tss[[i]] + window_bp
    counts[[i]] <- findInterval(upper, positions) - findInterval(lower - 1L, positions)
  }
  counts
}

annotated_variant_gene_keys <- function(variants, selected_pairs,
                                        variant_gene_column = NULL) {
  if (!nrow(variants) || !nrow(selected_pairs)) return(character())
  links <- merge(parse_peak_table(selected_pairs$peak),
                 unique(selected_pairs[, c("gene", "peak")]), by = "peak")
  link_gr <- table_to_granges(links, c("peak", "gene"))
  variant_gr <- table_to_granges(variants, "variant_id")
  if (!is.null(variant_gene_column)) {
    variant_gr$variant_gene <- variants[[variant_gene_column]]
  }
  common <- intersect(GenomeInfoDb::seqlevelsInUse(link_gr),
                      GenomeInfoDb::seqlevelsInUse(variant_gr))
  link_gr <- GenomeInfoDb::keepSeqlevels(link_gr, common, pruning.mode = "coarse")
  variant_gr <- GenomeInfoDb::keepSeqlevels(variant_gr, common, pruning.mode = "coarse")
  hits <- GenomicRanges::findOverlaps(link_gr, variant_gr, ignore.strand = TRUE)
  if (!length(hits)) return(character())
  link_gene <- link_gr$gene[S4Vectors::queryHits(hits)]
  keep <- rep(TRUE, length(hits))
  if (!is.null(variant_gene_column)) {
    keep <- link_gene == variant_gr$variant_gene[S4Vectors::subjectHits(hits)]
  }
  unique(paste(
    variant_gr$variant_id[S4Vectors::subjectHits(hits)[keep]],
    link_gene[keep], sep = "\r"
  ))
}

background_peak_overlap <- function(universe, background) {
  peaks <- parse_peak_table(universe$peak)
  peak_gr <- table_to_granges(peaks, "peak")
  background_gr <- table_to_granges(background, "variant_id")
  common <- intersect(GenomeInfoDb::seqlevelsInUse(peak_gr),
                      GenomeInfoDb::seqlevelsInUse(background_gr))
  peak_gr <- GenomeInfoDb::keepSeqlevels(peak_gr, common, pruning.mode = "coarse")
  background_gr <- GenomeInfoDb::keepSeqlevels(background_gr, common,
                                                pruning.mode = "coarse")
  hits <- GenomicRanges::findOverlaps(peak_gr, background_gr,
                                      ignore.strand = TRUE)
  if (!length(hits)) return(data.table::data.table())
  unique(data.table::data.table(
    peak = peak_gr$peak[S4Vectors::queryHits(hits)],
    variant_id = background_gr$variant_id[S4Vectors::subjectHits(hits)]
  ))
}

scent_gene_counts <- function(universe, selected, eqtl, background, tss,
                              pip_threshold = 0.2, window_bp = 500000L,
                              background_peak_hits = NULL,
                              background_cis_counts = NULL) {
  genes <- intersect(unique(universe$gene), unique(eqtl$gene))
  gene_tss <- tss[tss$gene %in% genes, , drop = FALSE]
  gene_tss <- gene_tss[!duplicated(gene_tss$gene), , drop = FALSE]
  eqtl <- eqtl[eqtl$pip >= pip_threshold & eqtl$gene %in% gene_tss$gene, , drop = FALSE]
  eqtl$key <- paste(eqtl$variant_id, eqtl$gene, sep = "\r")
  selected <- selected[selected$gene %in% gene_tss$gene, , drop = FALSE]

  causal_annotated <- annotated_variant_gene_keys(
    eqtl, selected, variant_gene_column = "gene"
  )
  if (is.null(background_peak_hits)) {
    background_annotated <- annotated_variant_gene_keys(background, selected)
  } else if (!nrow(selected)) {
    background_annotated <- character()
  } else {
    selected_dt <- data.table::as.data.table(unique(selected[, c("gene", "peak")]))
    background_annotated_table <- merge(
      data.table::as.data.table(background_peak_hits), selected_dt,
      by = "peak", allow.cartesian = TRUE
    )
    background_annotated <- unique(paste(
      background_annotated_table$variant_id,
      background_annotated_table$gene, sep = "\r"
    ))
  }

  causal_total <- table(factor(eqtl$gene, levels = gene_tss$gene))
  causal_hit_gene <- sub("^.*\r", "", causal_annotated)
  causal_hit <- table(factor(causal_hit_gene, levels = gene_tss$gene))
  background_hit_gene <- sub("^.*\r", "", background_annotated)
  background_hit <- table(factor(background_hit_gene, levels = gene_tss$gene))
  background_total <- if (is.null(background_cis_counts)) {
    count_positions_in_windows(background, gene_tss, window_bp)
  } else {
    cached_counts <- as.integer(background_cis_counts[gene_tss$gene])
    missing_counts <- is.na(cached_counts)
    if (any(missing_counts)) {
      cached_counts[missing_counts] <- count_positions_in_windows(
        background, gene_tss[missing_counts, , drop = FALSE], window_bp
      )
    }
    cached_counts
  }

  data.frame(
    gene = gene_tss$gene,
    causal_annotated = as.integer(causal_hit),
    causal_not_annotated = pmax(0L, as.integer(causal_total) - as.integer(causal_hit)),
    common_annotated = as.integer(background_hit),
    common_not_annotated = pmax(0L, background_total - as.integer(background_hit)),
    stringsAsFactors = FALSE
  )
}

odds_ratio_from_counts <- function(counts, correction = 0.5) {
  totals <- colSums(counts[, c("causal_annotated", "causal_not_annotated",
                               "common_annotated", "common_not_annotated"),
                           drop = FALSE])
  a <- totals[[1L]]; b <- totals[[2L]]; c <- totals[[3L]]; d <- totals[[4L]]
  corrected <- c(a, b, c, d) + correction
  margins_defined <- (a + b) > 0 && (c + d) > 0 &&
    (a + c) > 0 && (b + d) > 0
  fisher <- if (margins_defined) {
    stats::fisher.test(matrix(c(a, b, c, d), nrow = 2L), alternative = "greater")
  } else NULL
  data.frame(
    causal_annotated = a, causal_not_annotated = b,
    common_annotated = c, common_not_annotated = d,
    odds_ratio_ha = if (margins_defined) corrected[[1L]] * corrected[[4L]] /
      (corrected[[2L]] * corrected[[3L]]) else NA_real_,
    odds_ratio_fisher = if (is.null(fisher)) NA_real_ else unname(fisher$estimate),
    pval = if (is.null(fisher)) NA_real_ else fisher$p.value
  )
}

bootstrap_scent_odds <- function(counts, replicates = 1000L, seed = 1L) {
  eligible <- counts[(counts$causal_annotated + counts$causal_not_annotated) > 0, ]
  if (!nrow(eligible)) return(numeric())
  set.seed(seed)
  replicate(replicates, {
    sampled <- eligible[sample.int(nrow(eligible), nrow(eligible), replace = TRUE), ]
    totals <- colSums(sampled[, c("causal_annotated", "causal_not_annotated",
                                  "common_annotated", "common_not_annotated"),
                              drop = FALSE])
    a <- totals[[1L]]; b <- totals[[2L]]; c <- totals[[3L]]; d <- totals[[4L]]
    if ((a + b) == 0 || (c + d) == 0 || (a + c) == 0 || (b + d) == 0) {
      NA_real_
    } else {
      (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))
    }
  })
}

# Per-gene causal-variant enrichment used for the SCENT-style summary.
#
# The numerator is the fraction of causal eQTL variants that fall in the
# method's linked peaks.  The denominator is the corresponding fraction of
# common background variants.  This follows the SCENT paper's formula:
# (causal-in-annotation / common-in-annotation) /
# (causal-in-cis / common-in-cis).  Genes with no common background variant in
# the annotation have an undefined denominator and are omitted from the mean.
# The final SCENT-style statistic is the arithmetic mean of valid per-gene
# enrichments.
per_gene_scent_enrichment <- function(counts) {
  causal_total <- counts$causal_annotated + counts$causal_not_annotated
  common_total <- counts$common_annotated + counts$common_not_annotated
  eligible <- is.finite(causal_total) & is.finite(common_total) &
    causal_total > 0 & common_total > 0 & counts$common_annotated > 0
  out <- counts[eligible, c("gene", "causal_annotated", "causal_not_annotated",
                            "common_annotated", "common_not_annotated"),
                drop = FALSE]
  if (!nrow(out)) {
    out$causal_fraction <- numeric()
    out$common_fraction <- numeric()
    out$scent_enrichment <- numeric()
    return(out)
  }
  causal_total <- out$causal_annotated + out$causal_not_annotated
  common_total <- out$common_annotated + out$common_not_annotated
  out$causal_fraction <- out$causal_annotated / causal_total
  out$common_fraction <- out$common_annotated / common_total
  out$scent_enrichment <- out$causal_fraction / out$common_fraction
  out
}

# Run both summaries on the same gene-bootstrap draws.  This is useful for
# comparing the current pooled 2x2-table OR with the SCENT-style mean of
# per-gene enrichments without introducing Monte Carlo differences between the
# two metrics.
bootstrap_scent_metrics <- function(counts, replicates = 1000L, seed = 1L) {
  causal_total <- counts$causal_annotated + counts$causal_not_annotated
  common_total <- counts$common_annotated + counts$common_not_annotated
  pooled_eligible <- counts[
    is.finite(causal_total) & is.finite(common_total) &
      causal_total > 0 & common_total > 0,
    , drop = FALSE
  ]
  if (!nrow(pooled_eligible)) {
    return(data.frame(
      bootstrap_replicate = seq_len(replicates),
      pooled_or = NA_real_,
      scent_mean_enrichment = NA_real_
    ))
  }
  set.seed(seed)
  out <- lapply(seq_len(replicates), function(i) {
    sampled <- pooled_eligible[
      sample.int(nrow(pooled_eligible), nrow(pooled_eligible), replace = TRUE),
      , drop = FALSE
    ]
    pooled <- odds_ratio_from_counts(sampled)
    sampled_scent <- per_gene_scent_enrichment(sampled)
    scent_values <- sampled_scent$scent_enrichment[
      is.finite(sampled_scent$scent_enrichment)
    ]
    scent_mean <- if (length(scent_values)) mean(scent_values) else NA_real_
    data.frame(
      bootstrap_replicate = i,
      pooled_or = pooled$odds_ratio_ha[[1L]],
      scent_mean_enrichment = scent_mean
    )
  })
  data.table::rbindlist(out)
}

summarize_ranked_curve <- function(curve, recall_grid) {
  if (!nrow(curve)) return(data.frame())
  curve <- curve[order(curve$recall, curve$top_fraction), ]
  values <- vapply(recall_grid, function(target) {
    eligible <- which(curve$recall >= target)
    if (!length(eligible)) return(NA_real_)
    curve$enrichment[[eligible[[1L]]]]
  }, numeric(1L))
  data.frame(
    mean_enrichment_at_recall = mean(values, na.rm = TRUE),
    max_recall = max(curve$recall, na.rm = TRUE),
    enrichment_at_10pct_recall = values[[which.min(abs(recall_grid - 0.1))]],
    n_recall_grid_observed = sum(is.finite(values))
  )
}
