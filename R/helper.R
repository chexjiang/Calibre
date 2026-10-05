# Helper functions retained from incoming/helper.R.

#' Infer the chromosome/start/end separator used by a peak matrix.
#' Signac commonly sees either chr-start-end or chr:start-end feature names.
#' @param peaks Character vector of peak names.
#' @return A two-element separator vector for Signac::StringToGRanges.
InferPeakSeparator <- function(peaks) {
  peaks <- as.character(peaks)
  probe <- peaks[seq_len(min(length(peaks), 100L))]
  colon_format <- grepl("^[^:]+:[0-9]+-[0-9]+$", probe)
  hyphen_format <- grepl("^[^-]+-[0-9]+-[0-9]+$", probe)
  if (length(probe) && all(colon_format)) return(c(":", "-"))
  if (length(probe) && all(hyphen_format)) return(c("-", "-"))
  stop("Unsupported or mixed peak coordinate format; examples: ",
       paste(utils::head(peaks, 3L), collapse = ", "), call. = FALSE)
}

#' Convert peak IDs to the canonical format used by Signac assays.
#' @param peaks Character vector of peak names.
#' @return Peak names in chr-start-end form.
CanonicalPeakNames <- function(peaks) {
  peaks <- as.character(peaks)
  colon_format <- grepl("^[^:]+:[0-9]+-[0-9]+$", peaks)
  peaks[colon_format] <- sub(
    "^([^:]+):([0-9]+)-([0-9]+)$", "\\1-\\2-\\3", peaks[colon_format]
  )
  peaks
}

#' Collapse transcript ranges to one range per gene
#' @param ranges Gene annotation ranges.
#' @return A GRanges object.
#' @export
CollapseToLongestTranscript <- function(ranges) {
  names(ranges) <- NULL
  range.df <- as.data.frame(ranges)
  range.df$strand <- as.character(range.df$strand)
  range.df$strand <- ifelse(range.df$strand == "*", "+", range.df$strand)
  by_gene <- split(range.df, range.df$gene_id)
  collapsed <- do.call(rbind, lapply(by_gene, function(x) data.frame(
    gene_id = x$gene_id[[1]],
    seqnames = as.character(x$seqnames[[1]]),
    start = min(x$start),
    end = max(x$end),
    strand = x$strand[[1]],
    gene_biotype = x$gene_biotype[[1]],
    gene_name = x$gene_name[[1]],
    stringsAsFactors = FALSE
  )))
  rownames(collapsed) <- NULL
  collapsed$gene_name <- make.unique(collapsed$gene_name)
  GenomicRanges::makeGRangesFromDataFrame(collapsed, keep.extra.columns = TRUE)
}

#' Generate candidate gene-peak pairs
#' @param object Seurat multiome object.
#' @param annotation Optional gene annotation ranges.
#' @param expression.assay RNA assay name.
#' @param peak.assay ATAC assay name.
#' @param distance Window around each gene.
#' @param genes.use Optional genes to retain.
#' @return Candidate-pair data frame.
#' @export
GetGenePeakPairs <- function(object, annotation = NULL,
                             expression.assay = "RNA", peak.assay = "ATAC",
                             distance = 5e5, genes.use = NULL) {
  if (is.null(Signac::Annotation(object[[peak.assay]]))) {
    if (is.null(annotation)) {
      stop("annotation must be supplied when the ATAC assay has no annotation")
    }
    Signac::Annotation(object[[peak.assay]]) <- annotation
  } else {
    message("Using existing annotation from assay '", peak.assay, "'.")
  }
  annot <- Signac::Annotation(object[[peak.assay]])
  gene.range <- CollapseToLongestTranscript(annot)
  gene.window <- Signac::Extend(gene.range, upstream = distance, downstream = distance)
  peaks <- Signac::granges(object[[peak.assay]])
  if (is.null(names(peaks))) names(peaks) <- Signac::GRangesToString(peaks)
  common.chr <- intersect(GenomeInfoDb::seqlevels(peaks), GenomeInfoDb::seqlevels(gene.window))
  peaks <- GenomeInfoDb::keepSeqlevels(peaks, common.chr, pruning.mode = "coarse")
  gene.window <- GenomeInfoDb::keepSeqlevels(gene.window, common.chr, pruning.mode = "coarse")
  hits <- GenomicRanges::findOverlaps(peaks, gene.window)
  df <- data.frame(
    gene = gene.window$gene_name[S4Vectors::subjectHits(hits)],
    peak = names(peaks[S4Vectors::queryHits(hits)])
  )
  if (is.null(genes.use)) genes.use <- rownames(object[[expression.assay]])
  df <- df[df$gene %in% genes.use, , drop = FALSE]
  rownames(df) <- NULL
  message("Generated ", nrow(df), " gene-peak pairs (distance = +/-", distance / 1e3, " kb).")
  df
}

#' Filter candidate pairs using detection rate or abundance rank
#' @param object Seurat multiome object.
#' @param gene_peak_df Candidate-pair data frame.
#' @param expression.assay RNA assay name.
#' @param peak.assay ATAC assay name.
#' @param selection Either `"rate"` or `"rank"`.
#' @param expression.rank Number of genes for rank filtering.
#' @param peak.rank Number of peaks for rank filtering.
#' @param expression.rate Minimum RNA detection fraction.
#' @param peak.rate Minimum ATAC detection fraction.
#' @param verbose Report the number of retained pairs.
#' @return Filtered candidate-pair data frame.
#' @export
FilterGenePeakPairs <- function(object, gene_peak_df,
                                expression.assay = "RNA", peak.assay = "ATAC",
                                selection = "rate", expression.rank = 2000,
                                peak.rank = 20000, expression.rate = 0.05,
                                peak.rate = 0.05, verbose = TRUE) {
  mrna <- Seurat::GetAssayData(object, assay = expression.assay, layer = "counts")
  atac <- Seurat::GetAssayData(object, assay = peak.assay, layer = "counts")
  if (selection == "rate") {
    expr_detect <- Matrix::rowMeans(mrna > 0)
    atac_detect <- Matrix::rowMeans(atac > 0)
    gene_peak_df$expr_detect <- expr_detect[gene_peak_df$gene]
    gene_peak_df$atac_detect <- atac_detect[gene_peak_df$peak]
    keep_pairs <- with(gene_peak_df, expr_detect >= expression.rate & atac_detect >= peak.rate)
    gp_filt <- gene_peak_df[keep_pairs & !is.na(keep_pairs), , drop = FALSE]
  } else if (selection == "rank") {
    total_expr <- Matrix::rowSums(mrna)
    total_atac <- Matrix::rowSums(atac)
    top_genes <- names(sort(total_expr, decreasing = TRUE))[seq_len(min(expression.rank, length(total_expr)))]
    top_peaks <- names(sort(total_atac, decreasing = TRUE))[seq_len(min(peak.rank, length(total_atac)))]
    gp_filt <- gene_peak_df[gene_peak_df$gene %in% top_genes &
                              gene_peak_df$peak %in% top_peaks, , drop = FALSE]
  } else stop("selection must be either 'rate' or 'rank'")
  if (verbose) {
    message("Kept ", nrow(gp_filt), " out of ", nrow(gene_peak_df), " gene-peak pairs.")
  }
  gp_filt
}

#' Validate links against a genomic assay
#' @param res_df All tested links.
#' @param sig_assay_gr Reference assay GRanges.
#' @param pvar Significance column.
#' @param p_cutoff Significance cutoff.
#' @return Named enrichment statistics.
#' @export
validate_with_assay <- function(res_df, sig_assay_gr, pvar = "padj", p_cutoff = 0.2) {
  grange_peaks <- Signac::StringToGRanges(res_df$peak, sep = c("-", "-"))
  sig_res_inds <- res_df[[pvar]] < p_cutoff
  validate_with_assay1(res_df, res_df[sig_res_inds, , drop = FALSE], sig_assay_gr)
}

#' Validate an explicitly selected link set against a genomic assay
#' @param res_df All tested links.
#' @param sig_df Selected links.
#' @param sig_assay_gr Reference assay GRanges.
#' @return Named enrichment statistics.
#' @export
validate_with_assay1 <- function(res_df, sig_df, sig_assay_gr) {
  grange_peaks <- Signac::StringToGRanges(res_df$peak, sep = c("-", "-"))
  sig_res_inds <- seq_len(nrow(res_df)) %in% as.integer(rownames(sig_df))
  bovp <- as.data.frame(GenomicRanges::findOverlaps(grange_peaks, sig_assay_gr))
  sovp <- as.data.frame(GenomicRanges::findOverlaps(grange_peaks[sig_res_inds], sig_assay_gr))
  n_assay_sig <- sum(res_df$gene[sig_res_inds][sovp$queryHits] ==
                       sig_assay_gr$gene[sovp$subjectHits], na.rm = TRUE)
  n_assay <- sum(res_df$gene[bovp$queryHits] == sig_assay_gr$gene[bovp$subjectHits], na.rm = TRUE)
  res_gene_inds <- res_df$gene %in% sig_assay_gr$gene
  n_not_assay <- sum(res_gene_inds) - n_assay
  n_sig <- sum(sig_res_inds & res_gene_inds)
  pval <- stats::phyper(n_assay_sig, n_assay, n_not_assay, n_sig, lower.tail = FALSE)
  ratio <- n_assay_sig / n_sig
  enrichment <- ratio / (n_assay / (n_assay + n_not_assay))
  A <- n_assay_sig; B <- n_sig - A; C <- n_assay - A
  D <- (n_assay + n_not_assay) - (A + B + C)
  oddsratio <- (A * D) / (B * C)
  c(pval = pval, n_overlap = n_assay_sig, oddsratio = oddsratio,
    ratio = ratio, enr = enrichment)
}
