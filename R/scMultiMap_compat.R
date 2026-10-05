# scMultiMap 1.0.1 compatibility adapter.
#
# The upstream irls=FALSE path leaves theta_median undefined. 
scMultiMap_IRLS_compat <- function(X, seq_depth, bsample = NULL,
                                   irls = FALSE, verbose = FALSE) {
  if (irls) {
    return(scMultiMap::scMultiMap_IRLS(
      X, seq_depth, bsample = bsample, irls = TRUE, verbose = verbose
    ))
  }
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  if (nrow(X) != length(seq_depth)) {
    stop("The sequencing-depth length must match the number of cells.")
  }
  n_cell <- nrow(X)
  if (is.null(bsample)) {
    bsample_mtx <- matrix(1, nrow = n_cell, ncol = 1)
  } else {
    bsample <- factor(bsample, levels = unique(bsample))
    bsample_mtx <- stats::model.matrix(~ bsample - 1)
  }
  bind_mtx <- matrix(as.logical(bsample_mtx), nrow = n_cell,
                     ncol = ncol(bsample_mtx))
  seq_depth_sq <- seq_depth^2
  mu_mat <- t(crossprod(bind_mtx, X * seq_depth) /
                colSums(seq_depth_sq * bind_mtx))
  M <- tcrossprod(seq_depth * bind_mtx, mu_mat)
  centered <- X - M
  sigma_sq <- colSums((centered^2 - M) * seq_depth_sq) / sum(seq_depth^4)
  rownames(mu_mat) <- names(sigma_sq) <- colnames(X)
  theta_median <- rep(Inf, n_cell)
  if (ncol(bsample_mtx) == 1L) mu_mat <- mu_mat[, 1, drop = TRUE]
  if (ncol(bsample_mtx) == 1L) {
    list(mu = mu_mat, sigma_sq = sigma_sq, theta = unique(theta_median))
  } else {
    rownames(M) <- rownames(X)
    colnames(M) <- colnames(X)
    list(mu = mu_mat, sigma_sq = sigma_sq, theta = theta_median, M = M)
  }
}

run_scMultiMap_compat <- function(object, pairs_df, gene_assay = "RNA",
                                  peak_assay = "ATAC", bsample = NULL,
                                  irls = FALSE, verbose = FALSE) {
  if (!requireNamespace("scMultiMap", quietly = TRUE)) {
    stop("scMultiMap is not installed")
  }
  assay_names <- c(gene_assay, peak_assay)
  count_list <- list()
  depth_list <- list()
  for (index in seq_along(assay_names)) {
    assay <- assay_names[index]
    counts <- Seurat::GetAssayData(object, assay = assay, layer = "counts")
    depth_list[[assay]] <- Matrix::colSums(counts)
    features <- unique(pairs_df[[index]])
    present <- rownames(counts)[rownames(counts) %in% features]
    pairs_df <- pairs_df[pairs_df[[index]] %in% present, , drop = FALSE]
    count_list[[assay]] <- as.matrix(Matrix::t(counts[present, , drop = FALSE]))
  }
  if (!nrow(pairs_df)) stop("No configured pairs are present in the object")
  sample_vector <- NULL
  if (!is.null(bsample)) {
    sample_vector <- if (length(bsample) == ncol(object)) bsample else object@meta.data[[bsample]]
  }
  irls_results <- lapply(assay_names, function(assay) {
    scMultiMap_IRLS_compat(
      count_list[[assay]], depth_list[[assay]], bsample = sample_vector,
      irls = irls, verbose = verbose
    )
  })
  names(irls_results) <- assay_names
  results <- scMultiMap::scMultiMap_WLS(
    count_list, depth_list, pairs_df, irls_results, verbose = verbose
  )
  pair_key <- paste(pairs_df$gene, pairs_df$peak, sep = ".")
  result_key <- paste(results$gene, results$peak, sep = ".")
  results <- results[match(pair_key, result_key), , drop = FALSE]
  results$padj <- stats::p.adjust(results$pval, method = "BH")
  gene_var <- irls_results[[1]]$sigma_sq[
    match(pairs_df$gene, names(irls_results[[1]]$sigma_sq))]
  peak_var <- irls_results[[2]]$sigma_sq[
    match(pairs_df$peak, names(irls_results[[2]]$sigma_sq))]
  gene_var[gene_var <= 0] <- NA
  peak_var[peak_var <= 0] <- NA
  results$cor <- results$covar / sqrt(gene_var * peak_var)
  results$cor[is.na(results$cor)] <- 0
  results$cor <- pmax(pmin(results$cor, 1), -1)
  rownames(results) <- NULL
  results
}

