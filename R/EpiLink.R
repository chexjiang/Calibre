# EpiLink method implementation
#
# This file deliberately preserves the structure of incoming/EpiLink.R:
# binary_search() -> EpiLink_single() -> EpiLink() -> EpiLinkFilter().

# -----------------------------------------------------------
# Helper function: exhaustive FDR threshold search
# -----------------------------------------------------------
#' Calibrate the shared EpiLink FDR threshold
#' @param coef_real Non-negative observed statistics.
#' @param coef_syn Non-negative synthetic-null statistics.
#' @param q Target FDR.
#' @param conservative Use an offset of one instead of `q / 2`.
#' @param return_audit Return the full threshold table with the threshold.
#' @return A threshold, or a list containing threshold and audit table.
#' @export
binary_search <- function(coef_real, coef_syn, q, conservative = FALSE,
                          return_audit = FALSE) {
  if (!is.numeric(coef_real) || !is.numeric(coef_syn)) {
    stop("coef_real and coef_syn must be numeric")
  }
  if (length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1) {
    stop("q must be strictly between 0 and 1")
  }
  coef_real <- coef_real[is.finite(coef_real)]
  coef_syn <- coef_syn[is.finite(coef_syn)]
  if (!length(coef_real)) {
    out <- list(threshold = Inf, audit = data.frame())
    return(if (return_audit) out else out$threshold)
  }
  thresholds <- sort(unique(coef_real))
  # Count values >= every observed threshold without rescanning the full
  # vectors for each threshold. This is exactly equivalent to the former
  # vapply(..., sum(x >= threshold)) implementation, but O(n log n), not O(n^2).
  real_sorted <- sort(coef_real)
  syn_sorted <- sort(coef_syn)
  num_pos <- length(real_sorted) - findInterval(
    thresholds, real_sorted, left.open = TRUE
  )
  num_neg <- length(syn_sorted) - findInterval(
    thresholds, syn_sorted, left.open = TRUE
  )
  offset <- if (conservative) 1 else q / 2
  fdp <- (num_neg + offset) / pmax(num_pos, 1L)
  valid <- thresholds[fdp <= q]
  threshold <- if (length(valid)) min(valid) else Inf
  audit <- data.frame(threshold = thresholds, n_real = num_pos,
                      n_null = num_neg, fdp_hat = fdp)
  out <- list(threshold = threshold, audit = audit)
  if (return_audit) out else threshold
}

# -----------------------------------------------------------
# Fit one gene-peak pair and generate one synthetic null
# -----------------------------------------------------------
#' Fit one EpiLink gene-peak model
#' @param df Cell-level data containing `exprs`, `atac`, and covariates.
#' @param covariates Optional model terms.
#' @param null_generation Null generator (`"s"` is the analysis default).
#' @param family Either `"poisson"` or `"nb"`.
#' @param theory Return only the observed GLM fit.
#' @param verbose Print the model formula.
#' @param seed Pair-specific seed.
#' @return A one-row result data frame, or `NULL` for a failed fit.
#' @export
EpiLink_single <- function(
    df,
    covariates = NULL,
    null_generation = "s",
    family = "poisson",
    theory = FALSE,
    verbose = FALSE,
    seed = 123
) {
  if (!family %in% c("poisson", "nb")) stop("family must be 'poisson' or 'nb'")
  if (!null_generation %in% c("p", "n", "s0", "s")) {
    stop("Unknown null_generation type. Use 'p', 'n', 's0', or 's'.")
  }
  all_terms <- c("atac", covariates)
  formula_str <- paste("exprs ~", paste(all_terms, collapse = " + "))
  formula <- stats::as.formula(formula_str)
  if (verbose) message("Formula: ", formula_str)

  fit <- try(
    if (family == "poisson") {
      stats::glm(formula, family = stats::poisson(), data = df)
    } else {
      MASS::glm.nb(formula, data = df)
    },
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NULL)
  theta <- if (family == "nb") fit$theta else NA_real_

  coef_sum <- summary(fit)$coefficients
  if (!"atac" %in% rownames(coef_sum)) return(NULL)
  beta_real <- coef_sum["atac", "Estimate"]
  zval_real <- coef_sum["atac", "z value"]
  pval_real <- coef_sum["atac", "Pr(>|z|)"]
  fit_p <- stats::pchisq(fit$deviance, df = fit$df.residual, lower.tail = FALSE)
  disp <- sum(stats::residuals(fit, type = "pearson")^2) / fit$df.residual

  if (theory) {
    return(data.frame(beta = beta_real, zval = zval_real, pval = pval_real,
                      fit_p = fit_p, disp = disp, theta = theta,
                      stringsAsFactors = FALSE))
  }

  n_cell <- nrow(df)
  coefs <- stats::coef(fit)
  generate_synthetic <- function() {
    if (family == "poisson") {
      if (null_generation == "p") {
        coefs_null <- coefs
        coefs_null["atac"] <- 0
        X <- stats::model.matrix(formula, data = df)
        mu_null <- pmax(exp(drop(X %*% coefs_null)), 1e-8)
        Y_syn <- stats::rpois(n_cell, lambda = mu_null)
      } else if (null_generation == "n") {
        df$atac <- sample(df$atac)
        Y_syn <- df$exprs
      } else if (null_generation == "s0") {
        mu_full <- fit$fitted.values
        resid <- (fit$y - mu_full) / sqrt(mu_full)
        resid <- resid[is.finite(resid)]
        if (length(resid) < 2) return(NULL)
        mu_null <- pmax(mu_full / exp(df$atac * coefs["atac"]), 1e-8)
        resampled_idx <- sample.int(length(resid), size = n_cell, replace = TRUE)
        Y_syn <- pmax(round(mu_null + resid[resampled_idx] * sqrt(mu_null)), 0)
      } else if (null_generation == "s") {
        # (s) Standardized Pearson residual bootstrap (Moulton & Zeger, 1991)
        # Manual Pearson residuals for Poisson
        mu_full <- fit$fitted.values
        resid <- (fit$y - mu_full) / sqrt(mu_full)

        if (length(resid) < 2 || !all(is.finite(resid))) return(NULL)
        # Compute leverage values (h_i)
        h <- pmin(stats::lm.influence(fit, do.coef = FALSE)$hat, 0.999)
        sqrt_1_minus_h <- sqrt(1 - h)
        mu_null <- pmax(mu_full / exp(df$atac * coefs["atac"]), 1e-8)
        # Standardized Pearson residuals
        resid_std <- resid / sqrt_1_minus_h
        # Fast integer sampling
        resampled_idx <- sample.int(n_cell, replace = TRUE)
        # Generate synthetic counts
        Y_syn <- pmax(round(mu_null + resid_std[resampled_idx] * sqrt(mu_null) *
                             sqrt_1_minus_h), 0)
      }
    } else {
      if (is.na(theta) || theta <= 0) return(NULL)
      if (null_generation == "p") {
        coefs_null <- coefs
        coefs_null["atac"] <- 0
        X <- stats::model.matrix(formula, data = df)
        mu_null <- pmax(exp(drop(X %*% coefs_null)), 1e-8)
        Y_syn <- stats::rnbinom(n_cell, mu = mu_null, size = theta)
      } else if (null_generation == "n") {
        df$atac <- sample(df$atac)
        Y_syn <- df$exprs
      } else if (null_generation == "s0") {
        return(NULL)
      } else if (null_generation == "s") {
        # (s) Standardized Pearson residual bootstrap for NB
        mu_full <- fit$fitted.values
        var_full <- mu_full + mu_full^2 / theta
        resid <- (fit$y - mu_full) / sqrt(var_full)
        if (length(resid) < 2 || !all(is.finite(resid))) return(NULL)
        h <- pmin(stats::lm.influence(fit, do.coef = FALSE)$hat, 0.999)
        sqrt_1_minus_h <- sqrt(1 - h)
        mu_null <- pmax(mu_full / exp(df$atac * coefs["atac"]), 1e-8)
        resid_std <- resid / sqrt_1_minus_h
        resampled_idx <- sample.int(n_cell, replace = TRUE)
        var_null <- mu_null + mu_null^2 / theta
        Y_syn <- pmax(round(mu_null + resid_std[resampled_idx] * sqrt(var_null) *
                             sqrt_1_minus_h), 0)
      }
    }
    list(data = df, response = Y_syn)
  }

  synthetic <- with_epilink_seed(as.integer(seed), generate_synthetic())
  if (is.null(synthetic)) return(NULL)
  df <- synthetic$data
  df$exprs_syn <- synthetic$response
  formula_syn <- stats::update(formula, exprs_syn ~ .)
  fit_syn <- try(
    if (family == "poisson") {
      stats::glm(formula_syn, family = stats::poisson(), data = df)
    } else {
      MASS::glm.nb(formula_syn, data = df)
    },
    silent = TRUE
  )
  if (inherits(fit_syn, "try-error")) return(NULL)
  coef_syn <- summary(fit_syn)$coefficients
  if (!"atac" %in% rownames(coef_syn)) return(NULL)

  data.frame(
    beta = beta_real, zval = zval_real, pval = pval_real,
    fit_p = fit_p, disp = disp, theta = theta,
    beta_syn = coef_syn["atac", "Estimate"],
    zval_syn = coef_syn["atac", "z value"],
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------
# Main wrapper: run EpiLink on all pairs in a Seurat object
# -----------------------------------------------------------
#' Run EpiLink on a Seurat multiome object
#' @param object Seurat object with RNA and ATAC assays.
#' @param gene_peak_df Candidate pairs with `gene` and `peak` columns.
#' @param celltype_col Metadata cell-type column.
#' @param family Either `"poisson"` or `"nb"`.
#' @param celltype Optional cell type; `NULL` uses all cells.
#' @param covariates Optional cell-level model terms.
#' @param null_generation Null generator; formal analyses use `"s"`.
#' @param stat Selection statistic.
#' @param bin Binarize ATAC counts.
#' @param selection Candidate-pair filtering strategy.
#' @param returnObject Return the Seurat object instead of significant rows.
#' @param fdrcutoff Target FDR.
#' @param ncore Number of forked workers.
#' @param theory Return BH-ready observed GLM results without null generation.
#' @param conservative Use the conservative FDR numerator offset.
#' @param correction Apply the finite-cell correction in [EpiLinkFilter()].
#' @param verbose Report progress.
#' @param seed Master seed used to generate pair-specific seeds.
#' @return A Seurat object, significant-pair data frame, or theory result table.
#' @export
EpiLink <- function(
    object,
    gene_peak_df,
    celltype_col = "newCT",
    family = "poisson",
    celltype = NULL,
    covariates = NULL,
    null_generation = "s",
    stat = "zval",
    bin = TRUE,
    selection = "rate",
    returnObject = FALSE,
    fdrcutoff = 0.05,
    ncore = 10,
    theory = FALSE,
    conservative = FALSE,
    correction = FALSE,
    verbose = FALSE,
    seed = 123
) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (!("RNA" %in% Seurat::Assays(object))) stop("object must contain RNA assay")
  if (!("ATAC" %in% Seurat::Assays(object))) stop("object must contain ATAC assay")
  meta <- object@meta.data

  if (!is.null(celltype)) {
    if (!celltype_col %in% colnames(meta)) stop("celltype_col not found: ", celltype_col)
    keep_cells <- rownames(meta)[meta[[celltype_col]] == celltype]
    if (!length(keep_cells)) stop("No cells found for celltype = ", celltype)
    if (verbose) message("Using ", length(keep_cells), " cells from celltype = ", celltype)
  } else {
    keep_cells <- colnames(object)
    if (verbose) message("No celltype specified; using all ", length(keep_cells), " cells.")
  }
  n_cell <- length(keep_cells)
  meta_sub <- meta[keep_cells, , drop = FALSE]
  mrna <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")[, keep_cells, drop = FALSE]
  atac <- Seurat::GetAssayData(object, assay = "ATAC", layer = "counts")[, keep_cells, drop = FALSE]

  gene_peak_df$.pair_index_original <- seq_len(nrow(gene_peak_df))
  filter_object <- if (identical(keep_cells, colnames(object))) object else object[, keep_cells]
  gp_filt <- FilterGenePeakPairs(filter_object, gene_peak_df, selection = selection,
                                 verbose = verbose)
  n_pair <- nrow(gp_filt)
  if (!n_pair) stop("No gene-peak pairs remain after filtering")
  if (verbose) message("Kept ", n_pair, " gene-peak pairs after filtering.")
  pair_seeds <- epilink_pair_seeds(n_pair, seed)

  fit_one_pair <- function(j) {
    gene_name <- as.character(gp_filt[j, "gene"])
    peak_name <- as.character(gp_filt[j, "peak"])
    if (!gene_name %in% rownames(mrna) || !peak_name %in% rownames(atac)) return(NULL)
    expr_vec <- as.numeric(mrna[gene_name, ])
    atac_vec <- as.numeric(atac[peak_name, ])
    if (bin) atac_vec <- as.numeric(atac_vec > 0)
    df <- data.frame(exprs = expr_vec, atac = atac_vec, meta_sub, check.names = FALSE)
    res <- try(EpiLink_single(
      df, covariates, null_generation, family = family, theory = theory,
      seed = pair_seeds[j]
    ), silent = TRUE)
    if (inherits(res, "try-error") || is.null(res)) return(NULL)
    res$gene <- gene_name
    res$peak <- peak_name
    res$pair_index <- gp_filt$.pair_index_original[j]
    res$pair_seed <- pair_seeds[j]
    res
  }

  if (verbose) message("Fitting models in parallel (", ncore, " cores)...")
  if (ncore > 1 && .Platform$OS.type != "windows") {
    res_list <- pbmcapply::pbmclapply(seq_len(n_pair), fit_one_pair,
                                     mc.cores = ncore, mc.set.seed = FALSE)
  } else {
    res_list <- lapply(seq_len(n_pair), fit_one_pair)
  }
  res_df <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(res_df) || !nrow(res_df)) {
    warning("No valid results.")
    return(NULL)
  }
  res_df$padj <- stats::p.adjust(res_df$pval, method = "BH")
  object@misc$EpiLink <- list(
    res_df = res_df, n_cell = n_cell,
    celltype = if (is.null(celltype)) "AllCells" else celltype,
    null_generation = if (theory) "none" else null_generation,
    family = family, theory = theory, seed = seed, correction = correction
  )

  if (theory) return(res_df)
  EpiLinkFilter(object, stat = stat, fdrcutoff = fdrcutoff,
                correction = correction, returnObject = returnObject,
                conservative = conservative,
                verbose = verbose)
}

# -----------------------------------------------------------
# Differential peak-gene association
# -----------------------------------------------------------

EpiLinkDiff_cluster_vcov <- function(fit, cluster) {
  cluster <- droplevels(factor(cluster))
  n_cluster <- nlevels(cluster)
  n_obs <- length(cluster)
  design <- stats::model.matrix(fit)
  coefficients <- stats::coef(fit)
  if (n_cluster < 2L || nrow(design) != n_obs ||
      anyNA(coefficients) || !identical(colnames(design), names(coefficients))) {
    return(NULL)
  }

  y <- stats::model.response(stats::model.frame(fit))
  mu <- stats::fitted(fit)
  score_scale <- if (inherits(fit, "negbin")) {
    (y - mu) / (1 + mu / fit$theta)
  } else {
    y - mu
  }
  cluster_scores <- rowsum(design * score_scale, cluster, reorder = FALSE)
  bread <- try(stats::vcov(fit), silent = TRUE)
  if (inherits(bread, "try-error") || any(!is.finite(bread))) return(NULL)
  correction <- (n_cluster / (n_cluster - 1)) *
    ((n_obs - 1) / max(n_obs - ncol(design), 1))
  covariance <- correction * bread %*% crossprod(cluster_scores) %*% bread
  dimnames(covariance) <- list(names(coefficients), names(coefficients))
  covariance
}

#' Fit one differential EpiLink gene-peak model
#' @param df Cell-level data containing `exprs`, `atac`, `condition`, and
#'   `atac_condition`.
#' @param covariates Optional model terms.
#' @param null_generation Null generator (`"s"` by default; also `"p"` or
#'   `"s0"`).
#' @param family Either `"poisson"` or `"nb"`.
#' @param theory Return only the observed GLM fit.
#' @param individual Whether the individual fixed effect is present.
#' @param cluster Optional cluster labels used for cluster-robust variance;
#'   `NULL` or `FALSE` disables it.
#' @param seed Pair-specific seed.
#' @return A one-row result data frame, or `NULL` for a failed fit.
#' @noRd
EpiLinkDiff_single <- function(
    df,
    covariates = NULL,
    null_generation = "s",
    family = "poisson",
    theory = FALSE,
    individual = FALSE,
    cluster = NULL,
    seed = 123
) {
  if (!family %in% c("poisson", "nb")) stop("family must be 'poisson' or 'nb'")
  if (!null_generation %in% c("p", "s0", "s")) {
    stop(
      "EpiLinkDiff null_generation must be 'p', 's0', or 's'; ",
      "permutation ('n') does not preserve the shared peak-gene association"
    )
  }
  cluster_disabled <- isFALSE(cluster)
  if (cluster_disabled) cluster <- NULL
  if (!is.null(cluster) && length(cluster) != nrow(df)) {
    stop("cluster must have one label per row of df")
  }
  cluster <- if (is.null(cluster) && individual && !cluster_disabled) {
    df$individual_id
  } else {
    cluster
  }

  base_terms <- c(if (individual) "individual_id" else "condition_id",
                  covariates, "atac")
  full_formula <- stats::reformulate(c(base_terms, "atac_condition"), response = "exprs")
  null_formula <- stats::reformulate(base_terms, response = "exprs")

  full_fit <- try(if (family == "poisson") {
    stats::glm(full_formula, family = stats::poisson(), data = df)
  } else {
    MASS::glm.nb(full_formula, data = df)
  }, silent = TRUE)
  if (inherits(full_fit, "try-error")) return(NULL)

  get_diff <- function(fit) {
    coefficients <- stats::coef(fit)
    if (anyNA(coefficients[c("atac", "atac_condition")])) return(NULL)
    covariance <- if (!is.null(cluster)) {
      EpiLinkDiff_cluster_vcov(fit, cluster)
    } else {
      stats::vcov(fit)
    }
    if (is.null(covariance) || inherits(covariance, "try-error")) return(NULL)
    se_diff <- sqrt(covariance["atac_condition", "atac_condition"])
    if (!is.finite(se_diff) || se_diff <= 0) return(NULL)
    beta_control <- unname(coefficients["atac"])
    beta_diff <- unname(coefficients["atac_condition"])
    zval <- beta_diff / se_diff
    list(beta_control = beta_control, beta_diff = beta_diff,
         beta_disease = beta_control + beta_diff, se_diff = se_diff,
         zval = zval, pval = 2 * stats::pnorm(-abs(zval)))
  }

  observed <- get_diff(full_fit)
  if (is.null(observed)) return(NULL)

  fit_p <- stats::pchisq(
    full_fit$deviance, df = full_fit$df.residual, lower.tail = FALSE
  )
  disp <- sum(stats::residuals(full_fit, type = "pearson")^2) /
    full_fit$df.residual
  theta <- if (family == "nb") full_fit$theta else NA_real_
  result <- data.frame(
    beta = observed$beta_diff,
    zval = observed$zval,
    pval = observed$pval,
    fit_p = fit_p,
    disp = disp,
    theta = theta,
    beta_control = observed$beta_control,
    beta_diff = observed$beta_diff,
    beta_disease = observed$beta_disease,
    se_diff = observed$se_diff,
    wald_z = observed$zval,
    wald_p = observed$pval,
    stringsAsFactors = FALSE
  )
  if (theory) return(result)

  # The differential null removes only the interaction. The common ATAC slope,
  # condition/donor mean structure, covariates, and cell labels remain fixed.
  # With cluster labels, the residual bootstrap below uses one wild-bootstrap
  # weight per cluster so donor-level residual dependence is retained.
  null_fit <- try(if (family == "poisson") {
    stats::glm(null_formula, family = stats::poisson(), data = df)
  } else {
    MASS::glm.nb(null_formula, data = df)
  }, silent = TRUE)
  if (inherits(null_fit, "try-error")) return(NULL)
  synthetic_response <- with_epilink_seed(as.integer(seed), {
    mu_null <- pmax(stats::fitted(null_fit), 1e-8)
    theta_null <- if (family == "nb") null_fit$theta else NA_real_
    if (null_generation == "p") {
      if (family == "poisson") {
        stats::rpois(nrow(df), lambda = mu_null)
      } else {
        stats::rnbinom(nrow(df), mu = mu_null, size = theta_null)
      }
    } else {
      variance <- if (family == "poisson") mu_null else
        mu_null + mu_null^2 / theta_null
      residual <- (null_fit$y - mu_null) / sqrt(variance)
      scale <- if (null_generation == "s") {
        sqrt(1 - pmin(stats::lm.influence(null_fit, do.coef = FALSE)$hat, 0.999))
      } else rep(1, nrow(df))
      residual <- residual / scale
      if (!is.null(cluster)) {
        # Wild cluster bootstrap: all cells from one donor receive the same
        # Rademacher weight, preserving the donor-level residual structure.
        cluster_factor <- droplevels(factor(cluster))
        cluster_weight <- sample(c(-1, 1), nlevels(cluster_factor), replace = TRUE)
        sampled <- residual * cluster_weight[cluster_factor]
      } else {
        sampled <- residual[sample.int(nrow(df), nrow(df), replace = TRUE)]
      }
      pmax(round(mu_null + sampled * sqrt(variance) * scale), 0)
    }
  })
  synthetic_df <- df
  synthetic_df$exprs <- synthetic_response
  synthetic_fit <- try(if (family == "poisson") {
    stats::glm(full_formula, family = stats::poisson(), data = synthetic_df)
  } else {
    MASS::glm.nb(full_formula, data = synthetic_df)
  }, silent = TRUE)
  if (inherits(synthetic_fit, "try-error")) return(NULL)
  synthetic <- get_diff(synthetic_fit)
  if (is.null(synthetic)) return(NULL)
  result$beta_syn <- synthetic$beta_diff
  result$zval_syn <- synthetic$zval
  result
}

#' Test differential peak-gene associations between two conditions
#'
#' Fits a common peak-gene slope and a peak-by-condition interaction. The
#' interaction is the differential-association estimand. Its synthetic null is
#' generated from a reduced model that removes only that interaction.
#'
#' @inheritParams EpiLink
#' @param condition_col Metadata column containing exactly two conditions. The
#'   first factor level is the reference; for a non-factor column, the first
#'   observed value is the reference.
#' @param null_generation Null generator; formal analyses use `"s"`, matching
#'   [EpiLink()].
#' @param donor_col Optional metadata column containing donor/subject labels.
#'   It is required whenever `donor_effect` or `donor_robust` is enabled.
#' @param donor_effect Whether to include donor fixed effects. When enabled,
#'   ATAC is centered within donor and the condition main effect is omitted
#'   because donor is nested in condition.
#' @param donor_robust Whether to use donor-cluster-robust variance. The default
#'   `NULL` is automatic: it resolves to `TRUE` when `donor_col` is supplied
#'   and to `FALSE` otherwise. Set it explicitly to `FALSE` for model-based
#'   variance.
#' @return A Seurat object, significant-pair data frame, or theory result table.
#' @export
EpiLinkDiff <- function(
    object,
    gene_peak_df,
    celltype_col = "newCT",
    family = "poisson",
    celltype = NULL,
    condition_col,
    donor_col = NULL,
    donor_effect = FALSE,
    donor_robust = NULL,
    covariates = NULL,
    null_generation = "s",
    stat = "zval",
    bin = TRUE,
    selection = "rate",
    returnObject = FALSE,
    fdrcutoff = 0.05,
    ncore = 10,
    theory = FALSE,
    conservative = FALSE,
    correction = FALSE,
    verbose = FALSE,
    seed = 123
) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (!("RNA" %in% Seurat::Assays(object))) stop("object must contain RNA assay")
  if (!("ATAC" %in% Seurat::Assays(object))) stop("object must contain ATAC assay")
  meta <- object@meta.data
  if (!is.character(condition_col) || length(condition_col) != 1L ||
      !condition_col %in% colnames(meta)) {
    stop("condition_col must name one metadata column")
  }
  if (!is.logical(donor_effect) || length(donor_effect) != 1L ||
      is.na(donor_effect)) {
    stop("donor_effect must be TRUE or FALSE")
  }
  if (!is.null(donor_robust) &&
      (!is.logical(donor_robust) || length(donor_robust) != 1L ||
       is.na(donor_robust))) {
    stop("donor_robust must be TRUE, FALSE, or NULL")
  }
  if (!is.null(donor_col) &&
      (!is.character(donor_col) || length(donor_col) != 1L ||
       is.na(donor_col) || !donor_col %in% colnames(meta))) {
    stop("donor_col must be NULL or the name of a metadata column")
  }
  if (is.null(donor_robust)) donor_robust <- !is.null(donor_col)
  if (is.null(donor_col) && (isTRUE(donor_effect) || isTRUE(donor_robust))) {
    stop("donor_col is required when donor_effect or donor_robust is TRUE")
  }
  donor_effect <- isTRUE(donor_effect)
  donor_robust <- isTRUE(donor_robust)
  donor_effect_mode <- donor_effect
  if (condition_col %in% covariates ||
      (!is.null(donor_col) && donor_col %in% covariates)) {
    stop("Do not repeat condition or donor columns in covariates")
  }

  if (!is.null(celltype)) {
    if (!celltype_col %in% colnames(meta)) stop("celltype_col not found: ", celltype_col)
    keep_cells <- rownames(meta)[meta[[celltype_col]] == celltype]
    if (!length(keep_cells)) stop("No cells found for celltype = ", celltype)
  } else {
    keep_cells <- colnames(object)
  }
  n_cell <- length(keep_cells)
  meta_sub <- meta[keep_cells, , drop = FALSE]
  if (anyNA(meta_sub[[condition_col]])) {
    stop("condition_col contains missing values in the selected cells")
  }
  condition_raw <- meta_sub[[condition_col]]
  condition_factor <- if (is.factor(condition_raw)) {
    droplevels(condition_raw)
  } else {
    factor(condition_raw, levels = unique(condition_raw))
  }
  if (nlevels(condition_factor) != 2L) {
    stop("condition_col must contain exactly two conditions in the selected cells")
  }
  condition_levels <- levels(condition_factor)
  condition01 <- as.numeric(condition_factor) - 1

  if (!is.null(donor_col)) {
    if (anyNA(meta_sub[[donor_col]])) {
      stop("donor_col contains missing values in the selected cells")
    }
    donor_condition_count <- tapply(
      condition01, meta_sub[[donor_col]], function(x) length(unique(x))
    )
    if (any(donor_condition_count != 1L)) {
      stop("Each donor must belong to exactly one condition")
    }
    donors_by_condition <- tapply(
      meta_sub[[donor_col]], condition_factor,
      function(x) length(unique(x))
    )
    if (any(donors_by_condition < 1L)) {
      stop("Each condition must contain at least one donor")
    }
  }
  cluster_mode <- donor_robust
  cluster_sub <- if (cluster_mode) meta_sub[[donor_col]] else NULL
  if (cluster_mode) {
    if (anyNA(cluster_sub)) {
      stop("donor_col contains missing values in the selected cells")
    }
    if (length(unique(cluster_sub)) < 2L) {
      stop("donor_col must contain at least two donors")
    }
  }
  if (verbose) {
    message(
      "Testing ", condition_levels[[2L]], " versus ", condition_levels[[1L]],
      if (donor_effect_mode) paste0(" with ", donor_col, " fixed effects")
      else " with a condition main effect",
      " in ", n_cell, " cells."
    )
  }

  mrna <- Seurat::GetAssayData(
    object, assay = "RNA", layer = "counts"
  )[, keep_cells, drop = FALSE]
  atac <- Seurat::GetAssayData(
    object, assay = "ATAC", layer = "counts"
  )[, keep_cells, drop = FALSE]

  gene_peak_df$.pair_index_original <- seq_len(nrow(gene_peak_df))
  filter_object <- if (identical(keep_cells, colnames(object))) object else object[, keep_cells]
  gp_filt <- FilterGenePeakPairs(
    filter_object, gene_peak_df, selection = selection, verbose = verbose
  )
  n_pair <- nrow(gp_filt)
  if (!n_pair) stop("No gene-peak pairs remain after filtering")
  pair_seeds <- epilink_pair_seeds(n_pair, seed)

  fit_one_pair <- function(j) {
    gene_name <- as.character(gp_filt[j, "gene"])
    peak_name <- as.character(gp_filt[j, "peak"])
    if (!gene_name %in% rownames(mrna) || !peak_name %in% rownames(atac)) return(NULL)
    atac_vec <- as.numeric(atac[peak_name, ])
    if (bin) atac_vec <- as.numeric(atac_vec > 0)
    if (donor_effect_mode) {
      atac_vec <- atac_vec - ave(
        atac_vec, meta_sub[[donor_col]], FUN = function(x) mean(x)
      )
    }
    df <- data.frame(
      exprs = as.numeric(mrna[gene_name, ]),
      atac = atac_vec,
      condition_id = condition01,
      atac_condition = atac_vec * condition01,
      meta_sub,
      check.names = FALSE
    )
    if (donor_effect_mode) {
      df$individual_id <- factor(meta_sub[[donor_col]])
    }
    res <- try(
      EpiLinkDiff_single(
        df, covariates = covariates, null_generation = null_generation,
        family = family, theory = theory, individual = donor_effect_mode,
        cluster = if (cluster_mode) cluster_sub else FALSE,
        seed = pair_seeds[j]
      ),
      silent = TRUE
    )
    if (inherits(res, "try-error") || is.null(res)) return(NULL)
    res$gene <- gene_name
    res$peak <- peak_name
    res$reference_condition <- condition_levels[[1L]]
    res$comparison_condition <- condition_levels[[2L]]
    res$pair_index <- gp_filt$.pair_index_original[j]
    res$pair_seed <- pair_seeds[j]
    res
  }

  if (verbose) message("Fitting differential models (", ncore, " cores)...")
  if (ncore > 1L && .Platform$OS.type != "windows") {
    res_list <- pbmcapply::pbmclapply(
      seq_len(n_pair), fit_one_pair, mc.cores = ncore, mc.set.seed = FALSE
    )
  } else {
    res_list <- lapply(seq_len(n_pair), fit_one_pair)
  }
  res_df <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(res_df) || !nrow(res_df)) {
    warning("No valid differential results.")
    return(NULL)
  }
  res_df$padj <- stats::p.adjust(res_df$pval, method = "BH")
  if (!theory) {
    observed_order <- order(abs(res_df$zval), decreasing = TRUE)
    ordered_stat <- abs(res_df$zval[observed_order])
    null_stat <- abs(res_df$zval_syn[is.finite(res_df$zval_syn)])
    false_discoveries <- vapply(
      ordered_stat, function(threshold) sum(null_stat >= threshold), numeric(1)
    )
    raw_fdr <- pmin(1, false_discoveries / seq_along(ordered_stat))
    ordered_q <- rev(cummin(rev(raw_fdr)))
    res_df$calibrated_q <- NA_real_
    res_df$calibrated_q[observed_order] <- ordered_q
  }
  object@misc$EpiLink <- list(
    res_df = res_df,
    n_cell = n_cell,
    celltype = if (is.null(celltype)) "AllCells" else celltype,
    analysis = "differential",
    condition_col = condition_col,
    condition_levels = condition_levels,
    donor_col = donor_col,
    donor_effect = donor_effect,
    donor_robust = donor_robust,
    atac_centering = if (donor_effect_mode) "within_donor" else "none",
    variance = if (cluster_mode) "donor_cluster_robust_CR1" else "model_based",
    null_generation = if (theory) "none" else null_generation,
    null_bootstrap = if (theory || null_generation == "p") {
      if (theory) "none" else "parametric"
    } else if (cluster_mode) {
      "wild_cluster_pearson"
    } else {
      "cellwise_pearson"
    },
    family = family,
    theory = theory,
    seed = seed,
    correction = correction
  )

  if (theory) return(res_df)
  filtered_object <- EpiLinkFilter(
    object, stat = stat, fdrcutoff = fdrcutoff,
    correction = correction, returnObject = TRUE,
    conservative = conservative, verbose = verbose
  )
  threshold <- filtered_object@misc$EpiLink$threshold
  statistic <- switch(
    stat,
    zval = abs(filtered_object@misc$EpiLink$res_df$zval),
    beta = abs(filtered_object@misc$EpiLink$res_df$beta),
    pval = -log2(pmax(
      2 * stats::pnorm(-abs(filtered_object@misc$EpiLink$res_df$zval)),
      .Machine$double.xmin
    ))
  )
  filtered_object@misc$EpiLink$res_df$significant <-
    is.finite(statistic) & statistic >= threshold
  if (returnObject) return(filtered_object)
  filtered_object@misc$EpiLink$res_df[
    filtered_object@misc$EpiLink$res_df$significant, , drop = FALSE
  ]
}

# -----------------------------------------------------------
# Thresholding and significance filter
# -----------------------------------------------------------
#' Filter an EpiLink Seurat result
#' @param object Seurat object returned by [EpiLink()].
#' @param stat Selection statistic.
#' @param fdrcutoff Target FDR.
#' @param correction Apply the optional finite-cell correction.
#' @param returnObject Return the updated Seurat object.
#' @param conservative Use the conservative FDR numerator offset.
#' @param verbose Report the number of discoveries.
#' @return Updated Seurat object or significant-pair data frame.
#' @export
EpiLinkFilter <- function(
    object,
    stat = "pval",
    fdrcutoff = 0.05,
    correction = FALSE,
    returnObject = FALSE,
    conservative = FALSE,
    verbose = TRUE
) {
  res_df <- object@misc$EpiLink$res_df
  n_cell <- object@misc$EpiLink$n_cell
  if (is.null(res_df) || !nrow(res_df)) {
    warning("No results found in object@misc$EpiLink$res_df")
    return(if (returnObject) object else res_df)
  }
  if (stat == "zval") {
    coef_real <- abs(res_df$zval)
    coef_syn <- abs(res_df$zval_syn)
  } else if (stat == "beta") {
    coef_real <- abs(res_df$beta)
    coef_syn <- abs(res_df$beta_syn)
  } else if (stat == "pval") {
    coef_real <- -log2(pmax(2 * stats::pnorm(-abs(res_df$zval)), .Machine$double.xmin))
    coef_syn <- -log2(pmax(2 * stats::pnorm(-abs(res_df$zval_syn)), .Machine$double.xmin))
  } else {
    stop("stat must be either 'zval', 'pval', or 'beta'")
  }
  correct_factor <- if (correction) {
    1 / (1 + sqrt(log2(nrow(res_df)) / n_cell))
  } else 1
  calibrated <- binary_search(
    coef_real, coef_syn, q = fdrcutoff * correct_factor,
    conservative = conservative, return_audit = TRUE
  )
  sig_idx <- which(is.finite(coef_real) & coef_real >= calibrated$threshold)
  sig_df <- res_df[sig_idx, , drop = FALSE]
  if (verbose) message("Found ", length(sig_idx), " significant pairs (FDR <= ",
                       fdrcutoff, ") using ", stat, ".")
  object@misc$EpiLink$threshold <- calibrated$threshold
  object@misc$EpiLink$threshold_audit <- calibrated$audit
  object@misc$EpiLink$fdrcutoff <- fdrcutoff
  object@misc$EpiLink$stat <- stat
  object@misc$EpiLink$n_sig <- length(sig_idx)
  if (returnObject) object else sig_df
}
