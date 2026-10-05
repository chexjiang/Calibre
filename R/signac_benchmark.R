# Signac LinkPeaks adapter for the simulation benchmark.

run_signac_benchmark <- function(object, pairs_df, q = 0.1,
                                 distance = 20000, n_sample = 200L,
                                 verbose = FALSE) {
  if (!requireNamespace("Signac", quietly = TRUE) ||
      !requireNamespace("Seurat", quietly = TRUE)) {
    stop("Signac and Seurat are required for the Signac benchmark")
  }
  object <- Seurat::NormalizeData(
    object, assay = "RNA", normalization.method = "LogNormalize",
    verbose = FALSE
  )

  peaks <- rownames(object[["ATAC"]])
  peak_ranges <- Signac::StringToGRanges(peaks, sep = InferPeakSeparator(peaks))
  peak_metadata <- data.frame(
    GC.percent = 40 + (seq_along(peaks) %% 21) / 2,
    sequence.length = as.numeric(GenomicRanges::width(peak_ranges)),
    row.names = peaks
  )
  atac_assay <- Seurat::AddMetaData(object[["ATAC"]], metadata = peak_metadata)
  object[["ATAC"]] <- atac_assay

  gene_ranges <- peak_ranges[match(pairs_df$peak, peaks)]
  gene_ranges$gene_name <- pairs_df$gene
  gene_ranges$gene_id <- pairs_df$gene

  object <- suppressMessages(Signac::LinkPeaks(
    object = object,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    peak.slot = "counts",
    expression.slot = "data",
    method = "pearson",
    gene.coords = gene_ranges,
    distance = distance,
    min.cells = 0,
    genes.use = unique(pairs_df$gene),
    n_sample = as.integer(n_sample),
    pvalue_cutoff = 1,
    score_cutoff = 0,
    verbose = verbose
  ))
  links <- as.data.frame(Signac::Links(object[["ATAC"]]))
  if (!nrow(links)) {
    return(data.frame(gene = pairs_df$gene, peak = pairs_df$peak,
                      pval = NA_real_, padj = NA_real_))
  }
  link_key <- paste(links$gene, links$peak, sep = "\r")
  pair_key <- paste(pairs_df$gene, pairs_df$peak, sep = "\r")
  matched <- links[match(pair_key, link_key), , drop = FALSE]
  pval <- pmin(2 * matched$pvalue, 1)
  data.frame(
    gene = pairs_df$gene, peak = pairs_df$peak,
    pval = pval, padj = stats::p.adjust(pval, method = "BH"),
    stringsAsFactors = FALSE
  )
}

# Candidate-pair implementation of the same normalized-expression,
# matched-background statistic used by Signac::LinkPeaks. It avoids testing and
# garbage-collecting every incidental cis pair, which is prohibitive in a
# high-replicate simulation where the candidate pairs are fixed in advance.
run_signac_candidate_benchmark <- function(object, pairs_df, q = 0.1,
                                           n_sample = 200L,
                                           correct_pvalues = TRUE) {
  object <- Seurat::NormalizeData(
    object, assay = "RNA", normalization.method = "LogNormalize",
    verbose = FALSE
  )
  peak_counts <- Seurat::GetAssayData(object, assay = "ATAC", layer = "counts")
  expression <- Seurat::GetAssayData(object, assay = "RNA", layer = "data")
  peaks <- rownames(peak_counts)
  peak_ranges <- Signac::StringToGRanges(peaks, sep = InferPeakSeparator(peaks))
  meta_features <- data.frame(
    GC.percent = 40 + (seq_along(peaks) %% 21) / 2,
    count = as.numeric(Matrix::rowSums(peak_counts)),
    sequence.length = as.numeric(GenomicRanges::width(peak_ranges)),
    row.names = peaks
  )
  peak_matrix <- as.matrix(Matrix::t(peak_counts))
  gene_matrix <- as.matrix(Matrix::t(expression[pairs_df$gene, , drop = FALSE]))
  peak_index <- match(pairs_df$peak, peaks)
  chromosomes <- as.character(GenomicRanges::seqnames(peak_ranges))

  pval <- vapply(seq_len(nrow(pairs_df)), function(index) {
    observed <- suppressWarnings(stats::cor(
      peak_matrix[, peak_index[index]], gene_matrix[, index]
    ))
    trans_peaks <- peaks[chromosomes != chromosomes[peak_index[index]]]
    background <- suppressWarnings(Signac::MatchRegionStats(
      meta.feature = meta_features[trans_peaks, , drop = FALSE],
      query.feature = meta_features[pairs_df$peak[index], , drop = FALSE],
      features.match = c("GC.percent", "count", "sequence.length"),
      n = min(as.integer(n_sample), length(trans_peaks)), verbose = FALSE
    ))
    background_cor <- suppressWarnings(stats::cor(
      peak_matrix[, background, drop = FALSE], gene_matrix[, index]
    ))
    background_sd <- stats::sd(background_cor)
    if (!is.finite(observed) || !is.finite(background_sd) || background_sd <= 0) {
      return(NA_real_)
    }
    zscore <- (observed - mean(background_cor)) / background_sd
    pvalue <- stats::pnorm(-abs(zscore))
    if (isTRUE(correct_pvalues)) pvalue <- pmin(2 * pvalue, 1)
    pvalue
  }, numeric(1))
  data.frame(
    gene = pairs_df$gene, peak = pairs_df$peak,
    pval = pval, padj = stats::p.adjust(pval, method = "BH"),
    stringsAsFactors = FALSE
  )
}
