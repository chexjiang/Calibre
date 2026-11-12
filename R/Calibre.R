# Calibre method implementation
#
# This file preserves the implementation structure of the legacy method
# while exposing the Calibre API.

# -----------------------------------------------------------
# Helper function: exhaustive FDR threshold search
# -----------------------------------------------------------
#' Calibrate the shared Calibre FDR threshold
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
#' Fit one Calibre gene-peak model
#' @param df Cell-level data containing `exprs`, `atac`, and covariates.
#' @param covariates Optional model terms.
#' @param null_generation Null generator (`"s"` is the analysis default).
#' @param family Either `"poisson"` or `"nb"`.
#' @param theory Return only the observed GLM fit.
#' @param verbose Print the model formula.
#' @param seed Pair-specific seed.
#' @return A one-row result data frame, or `NULL` for a failed fit.
#' @export
Calibre_single <- function(
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

  synthetic <- with_calibre_seed(as.integer(seed), generate_synthetic())
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
# Main wrapper: run Calibre on all pairs in a Seurat object
# -----------------------------------------------------------
#' Run Calibre on a Seurat multiome object
#' @param object Seurat object with `RNA` and `ATAC` assays; required, with no
#'   default.
#' @param gene_peak_df Data frame containing `gene` and `peak` columns;
#'   required, with no default.
#' @param celltype_col Metadata column used to select cells by cell type. The
#'   default is `"newCT"`.
#' @param celltype Cell-type label to analyze. The default is `NULL`, which
#'   uses all cells.
#' @param covariates Optional cell-level covariates included in the GLM. The
#'   default is `NULL`.
#' @param fdrcutoff Target FDR level. The default is `0.1`.
#' @param ncore Number of parallel workers on supported platforms. The default
#'   is `10`.
#' @param stat Test statistic used to rank pairs for FDR thresholding. The
#'   default is `"zval"`; `"pval"` and `"beta"` are also supported.
#' @param family GLM family. The default is `"poisson"`, with `"nb"` as the
#'   negative binomial alternative.
#' @param null_generation Synthetic-null data generation scheme. The default,
#'   `"s"`, uses standardized Pearson-residual resampling.
#' @param bin Whether to binarize ATAC counts before model fitting. The default
#'   is `TRUE`.
#' @param selection Candidate-pair filtering strategy. The default is
#'   `"rate"`, which keeps pairs with both gene expression and CRE
#'   accessibility detected in at least 5% of cells.
#' @param returnObject Whether to return the updated Seurat object. The default
#'   is `FALSE`, which returns the significant CRE--gene data frame.
#' @param theory Whether to return observed GLM results without generating a
#'   synthetic null. The default is `FALSE`.
#' @param conservative Whether to use the more conservative FDR numerator
#'   offset. The default is `FALSE`.
#' @param correction Whether to apply the finite-cell correction during
#'   threshold calibration. The default is `FALSE`.
#' @param verbose Whether to report fitting progress and the number of
#'   discoveries. The default is `FALSE`.
#' @param seed Master seed for reproducible pair-specific seeds. The default is
#'   `123`.
#' @param progress_file Optional file to which parent-process progress messages
#'   are appended. The default is `NULL`.
#' @param progress_every Number of pairs processed between progress updates.
#'   The default is `10000L`.
#' @return A Seurat object, significant-pair data frame, or theory result table.
#' @export
Calibre <- function(
    object,
    gene_peak_df,
    celltype_col = "newCT",
    celltype = NULL,
    covariates = NULL,
    fdrcutoff = 0.1,
    ncore = 10,
    stat = "zval",
    family = "poisson",
    null_generation = "s",
    bin = TRUE,
    selection = "rate",
    returnObject = FALSE,
    theory = FALSE,
    conservative = FALSE,
    correction = FALSE,
    verbose = FALSE,
    seed = 123,
    progress_file = NULL,
    progress_every = 10000L
) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (!("RNA" %in% Seurat::Assays(object))) stop("object must contain RNA assay")
  if (!("ATAC" %in% Seurat::Assays(object))) stop("object must contain ATAC assay")
  if (!is.null(progress_file) &&
      (!is.character(progress_file) || length(progress_file) != 1L ||
       is.na(progress_file) || !nzchar(progress_file))) {
    stop("progress_file must be NULL or a non-empty file path")
  }
  if (!is.numeric(progress_every) || length(progress_every) != 1L ||
      !is.finite(progress_every) || progress_every < 1) {
    stop("progress_every must be a positive number")
  }
  progress_every <- as.integer(progress_every)
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
  pair_seeds <- calibre_pair_seeds(n_pair, seed)

  fit_one_pair <- function(j) {
    gene_name <- as.character(gp_filt[j, "gene"])
    peak_name <- as.character(gp_filt[j, "peak"])
    if (!gene_name %in% rownames(mrna) || !peak_name %in% rownames(atac)) return(NULL)
    expr_vec <- as.numeric(mrna[gene_name, ])
    atac_vec <- as.numeric(atac[peak_name, ])
    if (bin) atac_vec <- as.numeric(atac_vec > 0)
    df <- data.frame(exprs = expr_vec, atac = atac_vec, meta_sub, check.names = FALSE)
    res <- try(Calibre_single(
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

  progress_message <- function(done, total, prefix = "Completed") {
    line <- paste0(
      "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
      prefix, " celltype models: ", done, "/", total
    )
    if (verbose) message(line)
    if (!is.null(progress_file)) {
      cat(line, "\n", file = progress_file, append = TRUE)
    }
  }
  if (verbose) message("Fitting celltype models in parallel (", ncore, " cores)...")
  progress_message(0L, n_pair, prefix = "Starting")
  batch_size <- min(progress_every, n_pair)
  batch_starts <- seq.int(1L, n_pair, by = batch_size)
  batch_ends <- pmin(batch_starts + batch_size - 1L, n_pair)
  batch_results <- vector("list", length(batch_starts))
  for (i_batch in seq_along(batch_starts)) {
    batch_indices <- seq.int(batch_starts[[i_batch]], batch_ends[[i_batch]])
    if (ncore > 1L && .Platform$OS.type != "windows") {
      batch_results[[i_batch]] <- pbmcapply::pbmclapply(
        batch_indices, fit_one_pair, mc.cores = ncore,
        mc.set.seed = FALSE, ignore.interactive = TRUE
      )
    } else {
      batch_results[[i_batch]] <- lapply(batch_indices, fit_one_pair)
    }
    progress_message(batch_ends[[i_batch]], n_pair)
  }
  res_list <- unlist(batch_results, recursive = FALSE)
  res_df <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(res_df) || !nrow(res_df)) {
    warning("No valid results.")
    return(NULL)
  }
  res_df$padj <- stats::p.adjust(res_df$pval, method = "BH")
  object@misc$Calibre <- list(
    res_df = res_df, n_cell = n_cell,
    celltype = if (is.null(celltype)) "AllCells" else celltype,
    null_generation = if (theory) "none" else null_generation,
    family = family, theory = theory, seed = seed, correction = correction
  )

  if (theory) return(res_df)
  CalibreFilter(object, stat = stat, fdrcutoff = fdrcutoff,
                correction = correction, returnObject = returnObject,
                conservative = conservative,
                verbose = verbose)
}

# -----------------------------------------------------------
# Dynamic peak-gene association
# -----------------------------------------------------------

#' Fit one dynamic Calibre gene-peak model
#' @param df Cell-level data containing `exprs`, `atac`, and `pseudotime`.
#' @param pseudotime_formula Pseudotime model mode or smooth formula. Supported
#'   modes are `"linear"`, `"quadratic"`, and an `s(pseudotime, ...)` smooth.
#' @param covariates Optional model terms.
#' @param null_generation Null generator (`"s"` by default; also `"p"` or
#'   `"s0"`).
#' @param family Either `"poisson"` or `"nb"`.
#' @param theory Return only the observed GAM fit.
#' @param seed Pair-specific seed.
#' @param rescue_dispersion If non-NULL and the quadratic fit has dispersion
#'   above this value, refit the pair using `rescue_pseudotime_formula`.
#' @param rescue_pseudotime_formula Formula used for dispersion rescue.
#' @return A one-row result data frame, or `NULL` for a failed fit.
#' @noRd
CalibreDynamic_single <- function(
    df,
    pseudotime_formula = "quadratic",
    covariates = NULL,
    null_generation = "s",
    family = "poisson",
    theory = FALSE,
    seed = 123,
    rescue_dispersion = 1000,
    rescue_pseudotime_formula = "s(pseudotime, k = 5, bs = 'tp')",
    return_dynamic = FALSE,
    dynamic_grid = NULL
) {
  if (!family %in% c("poisson", "nb")) {
    stop("family must be 'poisson' or 'nb'")
  }
  if (!null_generation %in% c("p", "s0", "s")) {
    stop(
      "CalibreDynamic null_generation must be 's', 's0', or 'p'; ",
      "formal analyses use standardized residual resampling ('s')"
    )
  }
  if (!is.character(pseudotime_formula) || length(pseudotime_formula) != 1L ||
      is.na(pseudotime_formula) || !nzchar(trimws(pseudotime_formula))) {
    stop("pseudotime_formula must be 'linear', 'quadratic', or an s(...) term")
  }
  pseudotime_formula <- trimws(pseudotime_formula)
  if (!is.null(rescue_dispersion) &&
      (!is.numeric(rescue_dispersion) || length(rescue_dispersion) != 1L ||
       !is.finite(rescue_dispersion) || rescue_dispersion <= 0)) {
    stop("rescue_dispersion must be NULL or a positive number")
  }
  if (!is.null(rescue_pseudotime_formula) &&
      (!is.character(rescue_pseudotime_formula) ||
       length(rescue_pseudotime_formula) != 1L ||
       is.na(rescue_pseudotime_formula) ||
       !nzchar(trimws(rescue_pseudotime_formula)))) {
    stop("rescue_pseudotime_formula must be NULL or one non-empty formula")
  }
  pseudotime_mode <- if (pseudotime_formula %in% c("linear", "quadratic")) {
    pseudotime_formula
  } else {
    "spline"
  }
  if (identical(pseudotime_mode, "spline")) {
    if (!grepl("^s\\s*\\(", pseudotime_formula) ||
        !grepl("\\)\\s*$", pseudotime_formula) ||
        !grepl("pseudotime", pseudotime_formula, fixed = TRUE)) {
      stop("pseudotime_formula must be 'linear', 'quadratic', or an s(...) term")
    }
    if (grepl("\\bby\\s*=", pseudotime_formula)) {
      stop("pseudotime_formula must not already contain by=; CalibreDynamic adds by=atac")
    }
  }
  if (!is.numeric(df$pseudotime) || any(!is.finite(df$pseudotime))) {
    stop("pseudotime must be a finite numeric column")
  }

  # The common ATAC slope is identified at the mean observed pseudotime. For
  # linear and quadratic modes, center pseudotime explicitly so that `atac`
  # remains the association at that anchor. The dynamic null removes only the
  # accessibility-by-pseudotime terms.
  pseudotime_anchor_numeric <- mean(df$pseudotime)
  pseudotime_anchor <- format(
    pseudotime_anchor_numeric, digits = 17, scientific = TRUE
  )
  if (identical(pseudotime_mode, "spline")) {
    pseudotime_by_formula <- sub(
      "\\)\\s*$",
      paste0(", by = atac, pc = ", pseudotime_anchor, ")"),
      pseudotime_formula
    )
    full_terms <- c("atac", pseudotime_formula, pseudotime_by_formula, covariates)
    null_terms <- c("atac", pseudotime_formula, covariates)
  } else {
    df$pseudotime_centered <- df$pseudotime - pseudotime_anchor_numeric
    main_terms <- "pseudotime_centered"
    dynamic_terms <- "atac:pseudotime_centered"
    if (identical(pseudotime_mode, "quadratic")) {
      df$pseudotime_centered_sq <- df$pseudotime_centered^2
      main_terms <- c(main_terms, "pseudotime_centered_sq")
      dynamic_terms <- c(dynamic_terms, "atac:pseudotime_centered_sq")
    }
    full_terms <- c("atac", main_terms, dynamic_terms, covariates)
    null_terms <- c("atac", main_terms, covariates)
  }
  full_formula <- stats::as.formula(
    paste("exprs ~", paste(full_terms, collapse = " + "))
  )
  null_formula <- stats::as.formula(
    paste("exprs ~", paste(null_terms, collapse = " + "))
  )

  fit_model <- function(formula, data) {
    try(
      mgcv::gam(
        formula,
        family = if (family == "poisson") stats::poisson() else mgcv::nb(),
        data = data,
        method = "REML"
      ),
      silent = TRUE
    )
  }

  get_dynamic <- function(fit) {
    coefficients <- stats::coef(fit)
    if (!"atac" %in% names(coefficients) || is.na(coefficients[["atac"]])) {
      return(NULL)
    }
    if (identical(pseudotime_mode, "spline")) {
      fit_summary <- try(summary(fit), silent = TRUE)
      if (inherits(fit_summary, "try-error") || is.null(fit_summary$s.table)) {
        return(NULL)
      }
      smooth_table <- as.matrix(fit_summary$s.table)
      smooth_names <- rownames(smooth_table)
      dynamic_idx <- which(
        grepl("pseudotime", smooth_names, fixed = TRUE) &
          grepl("atac", smooth_names, fixed = TRUE)
      )
      if (length(dynamic_idx) != 1L ||
          !all(c("edf", "Chi.sq", "p-value") %in% colnames(smooth_table))) {
        return(NULL)
      }
      chi_sq <- unname(smooth_table[dynamic_idx, "Chi.sq"])
      dynamic_edf <- unname(smooth_table[dynamic_idx, "edf"])
      pval <- unname(smooth_table[dynamic_idx, "p-value"])
      dynamic_coefficient_names <- NULL
    } else {
      dynamic_idx <- which(
        grepl("pseudotime_centered", names(coefficients), fixed = TRUE) &
          grepl("atac", names(coefficients), fixed = TRUE)
      )
      if (!length(dynamic_idx)) return(NULL)
      covariance <- try(stats::vcov(fit), silent = TRUE)
      if (inherits(covariance, "try-error") ||
          !is.matrix(covariance) ||
          any(!is.finite(covariance[dynamic_idx, dynamic_idx, drop = FALSE]))) {
        return(NULL)
      }
      dynamic_coefficients <- coefficients[dynamic_idx]
      dynamic_covariance <- covariance[dynamic_idx, dynamic_idx, drop = FALSE]
      chi_sq <- try(
        as.numeric(crossprod(dynamic_coefficients,
                             solve(dynamic_covariance, dynamic_coefficients))),
        silent = TRUE
      )
      if (inherits(chi_sq, "try-error") || length(chi_sq) != 1L ||
          !is.finite(chi_sq) || chi_sq < 0) {
        return(NULL)
      }
      dynamic_edf <- length(dynamic_idx)
      pval <- stats::pchisq(chi_sq, df = dynamic_edf, lower.tail = FALSE)
      dynamic_coefficient_names <- names(dynamic_coefficients)
    }
    if (!all(is.finite(c(chi_sq, dynamic_edf, pval)))) return(NULL)
    pval <- min(max(pval, .Machine$double.xmin), 1)
    # CalibreFilter expects a z-like statistic. This is the two-sided normal
    # equivalent of the GAM smooth p-value; ranking is identical to p-value.
    # Use the upper-tail form to avoid 1 - p/2 rounding to 1 for small p-values.
    zval <- stats::qnorm(pval / 2, lower.tail = FALSE)
    list(
      beta = unname(coefficients[["atac"]]),
      zval = zval,
      pval = pval,
      chi_sq = chi_sq,
      dynamic_edf = dynamic_edf,
      dynamic_coefficient_names = dynamic_coefficient_names
    )
  }

  full_fit <- fit_model(full_formula, df)
  if (inherits(full_fit, "try-error")) return(NULL)
  observed <- get_dynamic(full_fit)
  if (is.null(observed)) return(NULL)

  variance <- try(full_fit$family$variance(full_fit$fitted.values), silent = TRUE)
  if (inherits(variance, "try-error") || any(!is.finite(variance)) ||
      any(variance <= 0)) {
    return(NULL)
  }
  fit_p <- stats::pchisq(
    full_fit$deviance, df = full_fit$df.residual, lower.tail = FALSE
  )
  disp <- sum((full_fit$y - full_fit$fitted.values)^2 / variance) /
    full_fit$df.residual
  theta <- if (family == "nb") {
    tryCatch(full_fit$family$getTheta(TRUE), error = function(e) NA_real_)
  } else {
    NA_real_
  }
  if (identical(pseudotime_mode, "quadratic") &&
      !is.null(rescue_dispersion) &&
      !is.null(rescue_pseudotime_formula) &&
      is.finite(disp) && disp > rescue_dispersion) {
    rescued <- CalibreDynamic_single(
      df = df,
      pseudotime_formula = rescue_pseudotime_formula,
      covariates = covariates,
      null_generation = null_generation,
      family = family,
      theory = theory,
      seed = seed,
      rescue_dispersion = NULL,
      rescue_pseudotime_formula = NULL,
      return_dynamic = return_dynamic,
      dynamic_grid = dynamic_grid
    )
    if (!is.null(rescued)) {
      attr(rescued, "rescue_used") <- TRUE
      return(rescued)
    }
  }
  result <- data.frame(
    beta = observed$beta,
    zval = observed$zval,
    pval = observed$pval,
    chi_sq = observed$chi_sq,
    # Legacy aliases retained for compatibility with earlier dynamic outputs.
    dynamic_stat = observed$chi_sq,
    dynamic_edf = observed$dynamic_edf,
    fit_p = fit_p,
    disp = disp,
    theta = theta,
    stringsAsFactors = FALSE
  )

  attach_dynamic_output <- function(result) {
    if (!isTRUE(return_dynamic)) return(result)
    dynamic_output_failure <- function(message) {
      attr(result, "dynamic_output_error") <- message
      result
    }

    coefficient <- stats::coef(full_fit)
    if (identical(pseudotime_mode, "spline")) {
      # Keep only the fitted dynamic smooth and the small covariance block
      # needed to reconstruct the effect at any future pseudotime.
      dynamic_smooth_index <- which(vapply(
        full_fit$smooth,
        function(smooth) {
          identical(as.character(smooth$by), "atac") &&
            grepl("pseudotime", smooth$label, fixed = TRUE)
        },
        logical(1)
      ))
      if (length(dynamic_smooth_index) != 1L) {
        return(dynamic_output_failure(
          "could not identify the pseudotime-by-ATAC smooth"
        ))
      }
      smooth <- full_fit$smooth[[dynamic_smooth_index]]
      smooth_indices <- seq.int(smooth$first.para, smooth$last.para)
      dynamic_coefficients <- coefficient[smooth_indices]
      parameter_names <- c("atac", names(dynamic_coefficients))
      smooth_label <- smooth$label
    } else {
      smooth <- NULL
      smooth_label <- NULL
      if (is.null(observed$dynamic_coefficient_names) ||
          !length(observed$dynamic_coefficient_names)) {
        return(dynamic_output_failure(
          "could not identify parametric dynamic coefficients"
        ))
      }
      parameter_names <- c("atac", observed$dynamic_coefficient_names)
      coefficient_indices <- match(parameter_names, names(coefficient))
      if (anyNA(coefficient_indices)) {
        return(dynamic_output_failure(
          "parametric dynamic coefficients are missing from the fitted model"
        ))
      }
      dynamic_coefficients <- coefficient[coefficient_indices[-1L]]
    }
    covariance <- full_fit$Vp
    covariance <- if (is.matrix(covariance) &&
                      nrow(covariance) == length(coefficient) &&
                      ncol(covariance) == length(coefficient)) {
      covariance_indices <- match(parameter_names, names(coefficient))
      covariance[covariance_indices, covariance_indices, drop = FALSE]
    } else {
      NULL
    }
    if (!is.finite(observed$beta) ||
        any(!is.finite(dynamic_coefficients))) {
      return(dynamic_output_failure(
        "dynamic component contains non-finite coefficients"
      ))
    }
    attr(result, "dynamic_parameters") <- list(
      beta = unname(observed$beta),
      coefficients = unname(dynamic_coefficients),
      coefficient_names = names(dynamic_coefficients),
      covariance = covariance,
      smooth = smooth,
      smooth_label = smooth_label,
      dynamic_basis = pseudotime_mode,
      pseudotime_formula = pseudotime_formula,
      pseudotime_anchor = pseudotime_anchor,
      family = family,
      link = full_fit$family$link,
      version = 1L
    )
    result
  }

  if (theory) return(attach_dynamic_output(result))

  # The dynamic null removes only the accessibility-by-pseudotime smooth.
  # The common ATAC slope, gene pseudotime trend, covariates, and cell labels
  # remain fixed. This is the key distinction from a null that removes ATAC.
  null_fit <- fit_model(null_formula, df)
  if (inherits(null_fit, "try-error")) return(NULL)
  synthetic_response <- with_calibre_seed(as.integer(seed), {
    mu_null <- pmax(stats::fitted(null_fit), 1e-8)
    theta_null <- if (family == "nb") {
      tryCatch(null_fit$family$getTheta(TRUE), error = function(e) NA_real_)
    } else {
      NA_real_
    }
    if (family == "nb" && (!is.finite(theta_null) || theta_null <= 0)) {
      return(NULL)
    }
    if (null_generation == "p") {
      if (family == "poisson") {
        stats::rpois(nrow(df), lambda = mu_null)
      } else {
        stats::rnbinom(nrow(df), mu = mu_null, size = theta_null)
      }
    } else {
      variance_null <- try(
        null_fit$family$variance(mu_null), silent = TRUE
      )
      if (inherits(variance_null, "try-error") ||
          any(!is.finite(variance_null)) || any(variance_null <= 0)) {
        return(NULL)
      }
      # Match the static Calibre residual Nullstrap: estimate residuals from
      # the full model, but recenter the synthetic response at the reduced-null
      # mean. This removes the fitted dynamic signal before resampling while
      # retaining the null model's baseline and pseudotime trend.
      mu_full <- pmax(stats::fitted(full_fit), 1e-8)
      variance_full <- try(
        full_fit$family$variance(mu_full), silent = TRUE
      )
      if (inherits(variance_full, "try-error") ||
          any(!is.finite(variance_full)) || any(variance_full <= 0)) {
        return(NULL)
      }
      residual <- (full_fit$y - mu_full) / sqrt(variance_full)
      leverage <- pmin(as.numeric(full_fit$hat), 0.999)
      scale <- if (null_generation == "s") {
        sqrt(pmax(1 - leverage, .Machine$double.eps))
      } else {
        rep(1, nrow(df))
      }
      residual_std <- residual / scale
      if (length(residual_std) < 2L ||
          any(!is.finite(residual_std))) return(NULL)
      sampled <- residual_std[sample.int(
        nrow(df), nrow(df), replace = TRUE
      )]
      pmax(round(mu_null + sampled * sqrt(variance_null) * scale), 0)
    }
  })
  if (is.null(synthetic_response)) return(NULL)

  synthetic_df <- df
  synthetic_df$exprs <- synthetic_response
  synthetic_fit <- fit_model(full_formula, synthetic_df)
  if (inherits(synthetic_fit, "try-error")) return(NULL)
  synthetic <- get_dynamic(synthetic_fit)
  if (is.null(synthetic)) return(NULL)
  result$beta_syn <- synthetic$beta
  result$zval_syn <- synthetic$zval
  result$pval_syn <- synthetic$pval
  result$chi_sq_syn <- synthetic$chi_sq
  # Legacy alias retained for compatibility with earlier dynamic outputs.
  result$dynamic_stat_syn <- synthetic$chi_sq
  result$dynamic_edf_syn <- synthetic$dynamic_edf
  attach_dynamic_output(result)
}

#' Run Calibre-dynamic on a Seurat multiome object
#' The input object is assumed to contain cells from one lineage; all cells in
#' the input object are used.
#' @param object A Seurat object containing `RNA` and `ATAC` assays, already
#'   subsetted to the pre-specified lineage.
#' @param gene_peak_df A data frame of candidate gene--CRE pairs.
#' @param pseudotime_col Metadata column containing finite numeric pseudotime
#'   values. The default is `"pseudotime"`.
#' @param pseudotime_formula Model for the dynamic effect. The default is
#'   `"quadratic"`.
#' @param rescue_pseudotime_formula Formula used for automatic dispersion
#'   rescue. The default is `"s(pseudotime, k = 5, bs = 'tp')"`.
#' @param covariates Optional cell-level covariates included in the GAM. The
#'   default is `NULL`.
#' @param fdrcutoff Target FDR level. The default is `0.1`.
#' @param ncore Number of parallel workers on supported platforms. The default
#'   is `10`.
#' @param stat Test statistic used to rank dynamic associations for FDR
#'   thresholding. The default is `"chisq"`.
#' @param family GLM family. The default is `"poisson"`, with `"nb"` as the
#'   negative binomial alternative.
#' @param null_generation Synthetic-null data generation scheme. The default,
#'   `"s"`, uses standardized Pearson-residual resampling.
#' @param bin Whether to binarize ATAC counts before model fitting. The default
#'   is `TRUE`.
#' @param selection Candidate-pair filtering strategy. The default is
#'   `"rate"`.
#' @param returnObject Whether to return the updated Seurat object. The default
#'   is `FALSE`, which returns the significant CRE--gene data frame.
#' @param theory Return only observed GAM results without null generation.
#'   Defaults to `FALSE`.
#' @param conservative Use the conservative FDR numerator offset. Defaults to
#'   `FALSE`.
#' @param correction Apply the finite-cell correction in [CalibreFilter()].
#'   Defaults to `FALSE`.
#' @param verbose Report progress. Defaults to `FALSE`.
#' @param seed Master seed used to generate pair-specific seeds. Defaults to
#'   `123`.
#' @param pair_seeds Optional pair-specific seeds aligned to the input rows of
#'   `gene_peak_df`. Defaults to `NULL`. This is useful when refitting a subset
#'   while preserving the seeds from a completed analysis.
#' @param rescue_dispersion For quadratic fits, pair-specific
#'   fits above this dispersion are refit using `rescue_pseudotime_formula`.
#'   Defaults to `1000`; set to `NULL` to disable automatic rescue.
#' @param progress_file Optional path to append parent-process progress logs.
#'   This is useful for long parallel runs because worker-level progress is not
#'   reliably visible in the main log. Defaults to `NULL`.
#' @param progress_every Number of pairs per progress batch. Defaults to
#'   10,000.
#' @param checkpoint_dir Optional directory for per-batch checkpoint files.
#'   Completed batches are reused automatically when the same batch file is
#'   present. Each checkpoint contains the fitted results and dynamic
#'   parameters for one progress batch. When supplied, completed batches are
#'   released from the parent process and reconstructed from disk only after
#'   all fitting batches have finished. Defaults to `NULL`.
#' @param return_dynamic Return compact dynamic parameters in addition
#'   to the usual Calibre result. These parameters can be evaluated at any
#'   pseudotime after fitting; no pair-by-cell table is created. Defaults to
#'   `FALSE`.
#' @param dynamic_grid Deprecated compatibility argument; it is ignored.
#'   Defaults to `100L`.
#' @return A Seurat object, significant-pair data frame, or theory result table.
#' @export
CalibreDynamic <- function(
    object,
    gene_peak_df,
    pseudotime_col = "pseudotime",
    pseudotime_formula = "quadratic",
    rescue_pseudotime_formula = "s(pseudotime, k = 5, bs = 'tp')",
    covariates = NULL,
    fdrcutoff = 0.1,
    ncore = 10,
    stat = "chisq",
    family = "poisson",
    null_generation = "s",
    bin = TRUE,
    selection = "rate",
    returnObject = FALSE,
    theory = FALSE,
    conservative = FALSE,
    correction = FALSE,
    verbose = FALSE,
    seed = 123,
    pair_seeds = NULL,
    rescue_dispersion = 1000,
    progress_file = NULL,
    progress_every = 10000L,
    checkpoint_dir = NULL,
    return_dynamic = FALSE,
    dynamic_grid = 100L
) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (!("RNA" %in% Seurat::Assays(object))) stop("object must contain RNA assay")
  if (!("ATAC" %in% Seurat::Assays(object))) stop("object must contain ATAC assay")
  if (!stat %in% c("chisq", "zval", "pval")) {
    stop("CalibreDynamic stat must be 'chisq', 'zval', or 'pval'")
  }
  meta <- object@meta.data
  if (!is.character(pseudotime_col) || length(pseudotime_col) != 1L ||
      is.na(pseudotime_col) || !pseudotime_col %in% colnames(meta)) {
    stop("pseudotime_col must name one metadata column")
  }
  if (!is.character(pseudotime_formula) || length(pseudotime_formula) != 1L ||
      is.na(pseudotime_formula) || !nzchar(trimws(pseudotime_formula))) {
    stop("pseudotime_formula must be 'linear', 'quadratic', or an s(...) term")
  }
  pseudotime_formula <- trimws(pseudotime_formula)
  if (!pseudotime_formula %in% c("linear", "quadratic")) {
    if (!grepl("^s\\s*\\(", pseudotime_formula) ||
        !grepl("\\)\\s*$", pseudotime_formula) ||
        !grepl("pseudotime", pseudotime_formula, fixed = TRUE)) {
      stop("pseudotime_formula must be 'linear', 'quadratic', or an s(...) term")
    }
    if (grepl("\\bby\\s*=", pseudotime_formula)) {
      stop("pseudotime_formula must not already contain by=; CalibreDynamic adds by=atac")
    }
  }
  if (pseudotime_col %in% covariates) {
    stop("Do not repeat pseudotime_col in covariates")
  }

  if (!is.logical(return_dynamic) || length(return_dynamic) != 1L ||
      is.na(return_dynamic)) {
    stop("return_dynamic must be TRUE or FALSE")
  }
  if (!is.null(progress_file) &&
      (!is.character(progress_file) || length(progress_file) != 1L ||
       is.na(progress_file) || !nzchar(progress_file))) {
    stop("progress_file must be NULL or a non-empty file path")
  }
  if (!is.numeric(progress_every) || length(progress_every) != 1L ||
      !is.finite(progress_every) || progress_every < 1) {
    stop("progress_every must be a positive number")
  }
  progress_every <- as.integer(progress_every)
  if (!is.null(checkpoint_dir) &&
      (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1L ||
       is.na(checkpoint_dir) || !nzchar(checkpoint_dir))) {
    stop("checkpoint_dir must be NULL or a non-empty directory path")
  }
  if (!is.null(pair_seeds) &&
      (!is.numeric(pair_seeds) || length(pair_seeds) != nrow(gene_peak_df) ||
       any(!is.finite(pair_seeds)) || any(pair_seeds < 1) ||
       any(pair_seeds != as.integer(pair_seeds)))) {
    stop("pair_seeds must be NULL or one positive integer seed per input pair")
  }
  if (!is.null(rescue_dispersion) &&
      (!is.numeric(rescue_dispersion) || length(rescue_dispersion) != 1L ||
       !is.finite(rescue_dispersion) || rescue_dispersion <= 0)) {
    stop("rescue_dispersion must be NULL or a positive number")
  }
  if (!is.null(rescue_pseudotime_formula) &&
      (!is.character(rescue_pseudotime_formula) ||
       length(rescue_pseudotime_formula) != 1L ||
       is.na(rescue_pseudotime_formula) ||
       !nzchar(trimws(rescue_pseudotime_formula)))) {
    stop("rescue_pseudotime_formula must be NULL or one non-empty formula")
  }

  keep_cells <- colnames(object)
  if (verbose) message("Using all ", length(keep_cells), " cells in the input object.")
  n_cell <- length(keep_cells)
  meta_sub <- meta[keep_cells, , drop = FALSE]
  if (!is.numeric(meta_sub[[pseudotime_col]]) ||
      any(!is.finite(meta_sub[[pseudotime_col]]))) {
    stop("pseudotime_col must contain finite numeric values in selected cells")
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
  if (verbose) message("Kept ", n_pair, " gene-peak pairs after filtering.")
  pair_seeds <- if (is.null(pair_seeds)) {
    calibre_pair_seeds(n_pair, seed)
  } else {
    as.integer(pair_seeds[gp_filt$.pair_index_original])
  }

  fit_one_pair <- function(j) {
    gene_name <- as.character(gp_filt[j, "gene"])
    peak_name <- as.character(gp_filt[j, "peak"])
    if (!gene_name %in% rownames(mrna) || !peak_name %in% rownames(atac)) return(NULL)
    atac_vec <- as.numeric(atac[peak_name, ])
    if (bin) atac_vec <- as.numeric(atac_vec > 0)
    df <- data.frame(
      exprs = as.numeric(mrna[gene_name, ]),
      atac = atac_vec,
      meta_sub,
      check.names = FALSE
    )
    df$pseudotime <- meta_sub[[pseudotime_col]]
    res <- try(
      CalibreDynamic_single(
        df,
        pseudotime_formula = pseudotime_formula,
        covariates = covariates,
        null_generation = null_generation,
        family = family,
        theory = theory,
        seed = pair_seeds[j],
        rescue_dispersion = rescue_dispersion,
        rescue_pseudotime_formula = rescue_pseudotime_formula,
        return_dynamic = return_dynamic
      ),
      silent = TRUE
    )
    if (inherits(res, "try-error") || is.null(res)) return(NULL)
    rescue_used <- isTRUE(attr(res, "rescue_used"))
    res$gene <- gene_name
    res$peak <- peak_name
    res$pair_index <- gp_filt$.pair_index_original[j]
    res$pair_seed <- pair_seeds[j]
    res$fit_mode <- if (rescue_used) {
      "spline"
    } else if (identical(pseudotime_formula, "quadratic")) {
      "quadratic"
    } else {
      "spline"
    }
    res$rescue_used <- rescue_used
    list(
      result = res,
      dynamic_parameters = attr(res, "dynamic_parameters")
    )
  }

  progress_message <- function(done, total, prefix = "Completed") {
    line <- paste0(
      "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
      prefix, " dynamic models: ", done, "/", total
    )
    if (verbose) message(line)
    if (!is.null(progress_file)) {
      cat(line, "\n", file = progress_file, append = TRUE)
    }
  }
  if (verbose) message("Fitting dynamic models in parallel (", ncore, " cores)...")
  progress_message(0L, n_pair, prefix = "Starting")
  batch_size <- min(progress_every, n_pair)
  batch_starts <- seq.int(1L, n_pair, by = batch_size)
  batch_ends <- pmin(batch_starts + batch_size - 1L, n_pair)
  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }
  batch_results <- if (is.null(checkpoint_dir)) {
    vector("list", length(batch_starts))
  } else {
    NULL
  }
  batch_files <- if (is.null(checkpoint_dir)) {
    NULL
  } else {
    rep(NA_character_, length(batch_starts))
  }
  for (i_batch in seq_along(batch_starts)) {
    batch_start <- batch_starts[[i_batch]]
    batch_end <- batch_ends[[i_batch]]
    checkpoint_file <- if (is.null(checkpoint_dir)) {
      NULL
    } else {
      file.path(
        checkpoint_dir,
        sprintf("dynamic_batch_%06d_%06d.rds", batch_start, batch_end)
      )
    }
    if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
      checkpoint <- try(readRDS(checkpoint_file), silent = TRUE)
      valid_checkpoint <- !inherits(checkpoint, "try-error") &&
        is.list(checkpoint) &&
        identical(checkpoint$version, checkpoint_version) &&
        identical(as.integer(checkpoint$n_pair), as.integer(n_pair)) &&
        identical(as.integer(checkpoint$first), as.integer(batch_start)) &&
        identical(as.integer(checkpoint$last), as.integer(batch_end)) &&
        is.list(checkpoint$results)
      if (valid_checkpoint) {
        batch_files[[i_batch]] <- checkpoint_file
        rm(checkpoint)
        invisible(gc())
        progress_message(batch_end, n_pair, prefix = "Loaded checkpoint for")
        next
      }
      warning("Ignoring invalid dynamic checkpoint: ", checkpoint_file)
      rm(checkpoint)
    }
    batch_indices <- seq.int(batch_start, batch_end)
    batch_result <- NULL
    if (ncore > 1 && .Platform$OS.type != "windows") {
      batch_result <- pbmcapply::pbmclapply(
        batch_indices, fit_one_pair,
        mc.cores = ncore, mc.set.seed = FALSE,
        ignore.interactive = TRUE
      )
    } else {
      batch_result <- lapply(batch_indices, fit_one_pair)
    }
    if (!is.null(checkpoint_file)) {
      saveRDS(
        list(
          version = checkpoint_version,
          n_pair = n_pair,
          first = batch_start,
          last = batch_end,
          results = batch_result
        ),
        checkpoint_file,
        compress = "gzip"
      )
      batch_files[[i_batch]] <- checkpoint_file
      rm(batch_result)
      invisible(gc())
    } else {
      batch_results[[i_batch]] <- batch_result
    }
    progress_message(batch_end, n_pair)
  }

  if (!is.null(checkpoint_dir)) {
    if (anyNA(batch_files) || !all(file.exists(batch_files))) {
      stop("one or more dynamic checkpoint files are missing")
    }
    result_tables <- vector("list", length(batch_files))
    parameter_list <- if (return_dynamic) vector("list", n_pair) else NULL
    n_parameter <- 0L
    for (i_batch in seq_along(batch_files)) {
      checkpoint <- readRDS(batch_files[[i_batch]])
      batch_result <- Filter(Negate(is.null), checkpoint$results)
      if (length(batch_result)) {
        result_tables[[i_batch]] <- do.call(
          rbind, lapply(batch_result, `[[`, "result")
        )
        if (return_dynamic) {
          for (one_result in batch_result) {
            parameters <- one_result$dynamic_parameters
            if (is.null(parameters)) next
            result <- one_result$result
            parameters$gene <- result$gene[[1L]]
            parameters$peak <- result$peak[[1L]]
            parameters$pair_index <- result$pair_index[[1L]]
            parameters$pair_seed <- result$pair_seed[[1L]]
            n_parameter <- n_parameter + 1L
            parameter_list[[n_parameter]] <- parameters
          }
        }
      }
      rm(checkpoint, batch_result)
      invisible(gc())
    }
    result_tables <- Filter(Negate(is.null), result_tables)
    res_df <- if (length(result_tables)) {
      do.call(rbind, result_tables)
    } else {
      NULL
    }
    if (return_dynamic) {
      parameter_list <- if (n_parameter) {
        parameter_list[seq_len(n_parameter)]
      } else {
        list()
      }
    }
  } else {
    res_list <- unlist(batch_results, recursive = FALSE)
    res_list <- Filter(Negate(is.null), res_list)
    res_df <- do.call(rbind, lapply(res_list, `[[`, "result"))
    if (return_dynamic) {
      parameter_list <- lapply(seq_along(res_list), function(i) {
        parameters <- res_list[[i]]$dynamic_parameters
        if (is.null(parameters)) return(NULL)
        result <- res_list[[i]]$result
        parameters$gene <- result$gene[[1L]]
        parameters$peak <- result$peak[[1L]]
        parameters$pair_index <- result$pair_index[[1L]]
        parameters$pair_seed <- result$pair_seed[[1L]]
        parameters
      })
      parameter_list <- Filter(Negate(is.null), parameter_list)
    }
  }
  if (is.null(res_df) || !nrow(res_df)) {
    warning("No valid dynamic results.")
    return(NULL)
  }
  dynamic_payload <- NULL
  if (return_dynamic) {
    dynamic_payload <- list(
      parameters = parameter_list,
      scale = "log_mean / linear_predictor",
      description = list(
        dynamic_component = "h(t) reconstructed from saved dynamic coefficients",
        total_effect = "beta + h(t)",
        dynamic_contribution = "not stored; evaluate h(t) for selected cells if needed",
        total_contribution = "not stored; evaluate beta + h(t) for selected cells if needed"
      ),
      version = 1L
    )
  }
  res_df$padj <- stats::p.adjust(res_df$pval, method = "BH")
  object@misc$Calibre <- list(
    res_df = res_df,
    n_cell = n_cell,
    celltype = "AllCells",
    analysis = "dynamic",
    pseudotime_col = pseudotime_col,
    pseudotime_formula = pseudotime_formula,
    rescue_dispersion = rescue_dispersion,
    rescue_pseudotime_formula = rescue_pseudotime_formula,
    null_generation = if (theory) "none" else null_generation,
    null_bootstrap = if (theory || null_generation == "p") {
      if (theory) "none" else "parametric"
    } else if (null_generation == "s") {
      "cellwise_standardized_pearson"
    } else {
      "cellwise_pearson"
    },
    family = family,
    theory = theory,
    seed = seed,
    correction = correction
  )

  if (return_dynamic) {
    object@misc$Calibre$dynamic_output <- dynamic_payload
  }
  if (theory) {
    if (return_dynamic) {
      return(list(
        results = res_df,
        dynamic_parameters = dynamic_payload$parameters
      ))
    }
    return(res_df)
  }
  filtered <- CalibreFilter(
    object, stat = stat, fdrcutoff = fdrcutoff,
    correction = correction, returnObject = returnObject,
    conservative = conservative, verbose = verbose
  )
  if (!return_dynamic) return(filtered)
  if (returnObject) {
    filtered@misc$Calibre$dynamic_output <- dynamic_payload
    filtered@misc$EpiLink <- filtered@misc$Calibre
    return(filtered)
  }
  list(
    selected = filtered,
    dynamic_parameters = dynamic_payload$parameters
  )
}

#' Evaluate a saved Calibre-dynamic spline parameter object
#'
#' @param dynamic_parameters One element of `dynamic_parameters` returned by
#'   [CalibreDynamic()].
#' @param pseudotime Numeric pseudotime values at which to evaluate the model.
#' @return A data frame containing the dynamic and total log-scale effects,
#'   their standard errors when available, and the fold change.
#' @export
CalibreDynamicPredict <- function(dynamic_parameters, pseudotime) {
  if (!is.list(dynamic_parameters) ||
      is.null(dynamic_parameters$coefficients) ||
      is.null(dynamic_parameters$beta)) {
    stop("dynamic_parameters is not a valid Calibre-dynamic parameter object")
  }
  pseudotime <- as.numeric(pseudotime)
  if (!length(pseudotime) || any(!is.finite(pseudotime))) {
    stop("pseudotime must contain at least one finite numeric value")
  }
  newdata <- data.frame(
    pseudotime = pseudotime,
    atac = 1
  )
  dynamic_basis <- if (is.null(dynamic_parameters$dynamic_basis)) {
    "spline"
  } else {
    as.character(dynamic_parameters$dynamic_basis)
  }
  if (identical(dynamic_basis, "spline")) {
    if (is.null(dynamic_parameters$smooth)) {
      stop("saved spline dynamic parameters do not contain a smooth object")
    }
    # PredictMat applies the same identifiability reparameterization that
    # mgcv::predict.gam uses. Predict.matrix() alone can return the unconstrained
    # basis (bs.dim columns) rather than the fitted coefficient basis (df cols).
    basis <- mgcv::PredictMat(dynamic_parameters$smooth, newdata)
  } else if (dynamic_basis %in% c("linear", "quadratic")) {
    anchor <- suppressWarnings(as.numeric(dynamic_parameters$pseudotime_anchor))
    if (length(anchor) != 1L || !is.finite(anchor)) {
      stop("parametric dynamic parameters do not contain a valid pseudotime anchor")
    }
    centered <- pseudotime - anchor
    basis <- matrix(centered, ncol = 1L)
    if (identical(dynamic_basis, "quadratic")) {
      basis <- cbind(basis, centered^2)
    }
  } else {
    stop("unknown saved dynamic basis: ", dynamic_basis)
  }
  coefficients <- as.numeric(dynamic_parameters$coefficients)
  if (ncol(basis) != length(coefficients)) {
    stop("saved spline basis and coefficient dimensions do not match")
  }
  dynamic_component <- drop(basis %*% coefficients)
  total_effect <- unname(dynamic_parameters$beta) + dynamic_component
  covariance <- dynamic_parameters$covariance
  effect_se <- rep(NA_real_, length(pseudotime))
  if (is.matrix(covariance) && nrow(covariance) == ncol(basis) + 1L) {
    contrast <- cbind(1, basis)
    effect_se <- sqrt(pmax(rowSums((contrast %*% covariance) * contrast), 0))
  }
  data.frame(
    pseudotime = pseudotime,
    dynamic_component = dynamic_component,
    total_effect = total_effect,
    effect_se = effect_se,
    effect_lower = total_effect - 1.96 * effect_se,
    effect_upper = total_effect + 1.96 * effect_se,
    fold_change = exp(total_effect),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------
# Differential peak-gene association
# -----------------------------------------------------------

CalibreDiff_cluster_vcov <- function(fit, cluster) {
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

#' Fit one differential Calibre gene-peak model
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
CalibreDiff_single <- function(
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
      "CalibreDiff null_generation must be 'p', 's0', or 's'; ",
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
      CalibreDiff_cluster_vcov(fit, cluster)
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
  synthetic_response <- with_calibre_seed(as.integer(seed), {
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

#' Test differential CRE-gene associations between two conditions or cell types
#'
#' Fits a common CRE-gene effect and a group interaction. The interaction is
#' the differential-association estimand. Its synthetic null is generated from
#' a reduced model that removes only that interaction. Use either
#' `condition_diff` or `celltype_diff`, but not both. The placeholder defaults
#' must be replaced with the labels in the metadata, and the unused comparison
#' vector must be set to `NULL`.
#'
#' @inheritParams Calibre
#' @param condition_col Metadata column containing condition labels. The
#'   default is `"condition"`.
#' @param condition_diff The two condition labels to compare. The second label
#'   is compared with the first. The default is
#'   `c("condition1", "condition2")`; set to `NULL` for a cell-type
#'   comparison.
#' @param celltype_col Metadata column containing cell-type labels. The default
#'   is `"newCT"`.
#' @param celltype Optional single cell type used to subset cells. The default
#'   is `NULL`, which uses all cell types.
#' @param condition Optional condition label or labels used to subset cells.
#'   The default is `NULL`, which uses all conditions.
#' @param celltype_diff The two cell-type labels to compare. The second label
#'   is compared with the first. The default is
#'   `c("celltype1", "celltype2")`; set to `NULL` for a condition
#'   comparison.
#' @param null_generation Null generator; formal analyses use `"s"`, matching
#'   [Calibre()].
#' @param donor_col Optional metadata column containing donor/subject labels.
#'   It is required whenever `donor_effect` or `donor_robust` is enabled.
#' @param donor_effect Whether to include donor fixed effects. When enabled,
#'   the condition main effect is omitted because donor is nested in condition.
#' @param donor_robust Whether to use donor-cluster-robust variance. The default
#'   `NULL` is automatic: it resolves to `TRUE` when `donor_col` is supplied
#'   and to `FALSE` otherwise. Set it explicitly to `FALSE` for model-based
#'   variance.
#' @param progress_file Optional path to append parent-process progress logs.
#'   This is useful for long parallel runs because worker-level progress is not
#'   reliably visible in the main log.
#' @param progress_every Number of pairs per progress batch. Defaults to
#'   10,000.
#' @return A Seurat object, significant-pair data frame, or theory result table.
#' @export
CalibreDiff <- function(
    object,
    gene_peak_df,
    condition_col = "condition",
    condition_diff = c("condition1", "condition2"),
    celltype_col = "newCT",
    celltype = NULL,
    condition = NULL,
    celltype_diff = c("celltype1", "celltype2"),
    covariates = NULL,
    fdrcutoff = 0.1,
    ncore = 10,
    stat = "zval",
    family = "poisson",
    null_generation = "s",
    bin = TRUE,
    selection = "rate",
    returnObject = FALSE,
    theory = FALSE,
    conservative = FALSE,
    correction = FALSE,
    verbose = FALSE,
    seed = 123,
    donor_col = NULL,
    donor_effect = FALSE,
    donor_robust = NULL,
    progress_file = NULL,
    progress_every = 10000L
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
  if (!is.null(progress_file) &&
      (!is.character(progress_file) || length(progress_file) != 1L ||
       is.na(progress_file) || !nzchar(progress_file))) {
    stop("progress_file must be NULL or a non-empty file path")
  }
  if (!is.numeric(progress_every) || length(progress_every) != 1L ||
      !is.finite(progress_every) || progress_every < 1) {
    stop("progress_every must be a positive number")
  }
  progress_every <- as.integer(progress_every)
  if (!is.null(condition) &&
      (!is.character(condition) || !length(condition) || anyNA(condition))) {
    stop("condition must be NULL or a non-empty character vector")
  }
  if (!is.null(celltype_diff) &&
      (!is.character(celltype_diff) || length(celltype_diff) != 2L ||
       anyNA(celltype_diff) || identical(celltype_diff[[1L]], celltype_diff[[2L]]))) {
    stop("celltype_diff must be NULL or a character vector of two different cell types")
  }
  if (!is.null(condition_diff) &&
      (!is.character(condition_diff) || length(condition_diff) != 2L ||
       anyNA(condition_diff) || identical(condition_diff[[1L]], condition_diff[[2L]]))) {
    stop("condition_diff must be NULL or a character vector of two different conditions")
  }
  if (!is.null(celltype_diff) && !is.null(celltype)) {
    stop("celltype and celltype_diff cannot be used together")
  }
  if (!is.null(celltype_diff) && !is.null(condition_diff)) {
    stop("celltype_diff and condition_diff cannot be used together")
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

  celltype_diff_mode <- !is.null(celltype_diff)
  if (!is.null(celltype)) {
    if (!celltype_col %in% colnames(meta)) stop("celltype_col not found: ", celltype_col)
    keep_cells <- rownames(meta)[meta[[celltype_col]] == celltype]
    if (!length(keep_cells)) stop("No cells found for celltype = ", celltype)
  } else if (celltype_diff_mode) {
    if (!celltype_col %in% colnames(meta)) stop("celltype_col not found: ", celltype_col)
    keep_mask <- as.character(meta[[celltype_col]]) %in% celltype_diff
    if (!is.null(condition)) {
      keep_mask <- keep_mask & as.character(meta[[condition_col]]) %in% condition
    }
    keep_cells <- rownames(meta)[keep_mask]
  } else {
    keep_mask <- rep(TRUE, nrow(meta))
    if (!is.null(condition)) {
      keep_mask <- as.character(meta[[condition_col]]) %in% condition
    }
    if (!is.null(condition_diff)) {
      keep_mask <- keep_mask &
        as.character(meta[[condition_col]]) %in% condition_diff
    }
    keep_cells <- rownames(meta)[keep_mask]
  }
  if (!length(keep_cells)) stop("No cells remain after the requested filtering")
  n_cell <- length(keep_cells)
  meta_sub <- meta[keep_cells, , drop = FALSE]
  if (anyNA(meta_sub[[condition_col]])) {
    stop("condition_col contains missing values in the selected cells")
  }
  comparison_col <- if (celltype_diff_mode) celltype_col else condition_col
  comparison_raw <- meta_sub[[comparison_col]]
  requested_levels <- if (celltype_diff_mode) {
    celltype_diff
  } else if (!is.null(condition_diff)) {
    condition_diff
  } else {
    NULL
  }
  condition_factor <- if (!is.null(requested_levels)) {
    factor(as.character(comparison_raw), levels = requested_levels)
  } else if (is.factor(comparison_raw)) {
    droplevels(comparison_raw)
  } else {
    factor(comparison_raw, levels = unique(comparison_raw))
  }
  if (nlevels(condition_factor) != 2L) {
    stop(
      if (celltype_diff_mode) {
        "celltype_diff must contain cells from both requested cell types"
      } else {
        "condition_col must contain exactly two conditions in the selected cells"
      }
    )
  }
  condition_levels <- levels(condition_factor)
  condition01 <- as.numeric(condition_factor) - 1

  if (!is.null(donor_col)) {
    if (anyNA(meta_sub[[donor_col]])) {
      stop("donor_col contains missing values in the selected cells")
    }
    if (!celltype_diff_mode) {
      donor_condition_count <- tapply(
        condition01, meta_sub[[donor_col]], function(x) length(unique(x))
      )
      if (any(donor_condition_count != 1L)) {
        stop("Each donor must belong to exactly one condition")
      }
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
      else if (celltype_diff_mode) " with a celltype main effect"
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
  pair_seeds <- calibre_pair_seeds(n_pair, seed)

  fit_one_pair <- function(j) {
    gene_name <- as.character(gp_filt[j, "gene"])
    peak_name <- as.character(gp_filt[j, "peak"])
    if (!gene_name %in% rownames(mrna) || !peak_name %in% rownames(atac)) return(NULL)
    atac_vec <- as.numeric(atac[peak_name, ])
    if (bin) atac_vec <- as.numeric(atac_vec > 0)
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
      CalibreDiff_single(
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

  progress_message <- function(done, total, prefix = "Completed") {
    line <- paste0(
      "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
      prefix, " differential models: ", done, "/", total
    )
    if (verbose) message(line)
    if (!is.null(progress_file)) {
      cat(line, "\n", file = progress_file, append = TRUE)
    }
  }
  if (verbose) message("Fitting differential models (", ncore, " cores)...")
  progress_message(0L, n_pair, prefix = "Starting")
  batch_size <- min(progress_every, n_pair)
  batch_starts <- seq.int(1L, n_pair, by = batch_size)
  batch_ends <- pmin(batch_starts + batch_size - 1L, n_pair)
  batch_results <- vector("list", length(batch_starts))
  for (i_batch in seq_along(batch_starts)) {
    batch_indices <- seq.int(batch_starts[[i_batch]], batch_ends[[i_batch]])
    if (ncore > 1L && .Platform$OS.type != "windows") {
      batch_results[[i_batch]] <- pbmcapply::pbmclapply(
        batch_indices, fit_one_pair, mc.cores = ncore,
        mc.set.seed = FALSE, ignore.interactive = TRUE
      )
    } else {
      batch_results[[i_batch]] <- lapply(batch_indices, fit_one_pair)
    }
    progress_message(batch_ends[[i_batch]], n_pair)
  }
  res_list <- unlist(batch_results, recursive = FALSE)
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
  object@misc$Calibre <- list(
    res_df = res_df,
    n_cell = n_cell,
    celltype = if (celltype_diff_mode) {
      paste(celltype_diff, collapse = "_vs_")
    } else if (is.null(celltype)) "AllCells" else celltype,
    analysis = "differential",
    comparison_type = if (celltype_diff_mode) "celltype" else "condition",
    comparison_col = comparison_col,
    condition_filter = condition,
    celltype_diff = celltype_diff,
    condition_diff = condition_diff,
    condition_col = condition_col,
    condition_levels = condition_levels,
    donor_col = donor_col,
    donor_effect = donor_effect,
    donor_robust = donor_robust,
    atac_centering = "none",
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
  filtered_object <- CalibreFilter(
    object, stat = stat, fdrcutoff = fdrcutoff,
    correction = correction, returnObject = TRUE,
    conservative = conservative, verbose = verbose
  )
  threshold <- filtered_object@misc$Calibre$threshold
  statistic <- switch(
    stat,
    chisq = filtered_object@misc$Calibre$res_df$chi_sq,
    zval = abs(filtered_object@misc$Calibre$res_df$zval),
    pval = -log2(pmax(
      2 * stats::pnorm(-abs(filtered_object@misc$Calibre$res_df$zval)),
      .Machine$double.xmin
    ))
  )
  filtered_object@misc$Calibre$res_df$significant <-
    is.finite(statistic) & statistic >= threshold
  if (returnObject) return(filtered_object)
  filtered_object@misc$Calibre$res_df[
    filtered_object@misc$Calibre$res_df$significant, , drop = FALSE
  ]
}

# -----------------------------------------------------------
# Thresholding and significance filter
# -----------------------------------------------------------
#' Filter an Calibre Seurat result
#' @param object Seurat object returned by [Calibre()].
#' @param stat Selection statistic.
#' @param fdrcutoff Target FDR.
#' @param correction Apply the optional finite-cell correction.
#' @param returnObject Return the updated Seurat object.
#' @param conservative Use the conservative FDR numerator offset.
#' @param verbose Report the number of discoveries.
#' @return Updated Seurat object or significant-pair data frame.
#' @export
CalibreFilter <- function(
    object,
    stat = "pval",
    fdrcutoff = 0.05,
    correction = FALSE,
    returnObject = FALSE,
    conservative = FALSE,
    verbose = TRUE
) {
  object <- .calibre_sync_object(object)
  res_df <- object@misc$Calibre$res_df
  n_cell <- object@misc$Calibre$n_cell
  if (is.null(res_df) || !nrow(res_df)) {
    warning("No results found in object@misc$Calibre$res_df")
    return(if (returnObject) object else res_df)
  }
  if (stat == "chisq") {
    coef_real <- res_df$chi_sq
    coef_syn <- res_df$chi_sq_syn
  } else if (stat == "zval") {
    coef_real <- abs(res_df$zval)
    coef_syn <- abs(res_df$zval_syn)
  } else if (stat == "beta") {
    coef_real <- abs(res_df$beta)
    coef_syn <- abs(res_df$beta_syn)
  } else if (stat == "pval") {
    coef_real <- -log2(pmax(2 * stats::pnorm(-abs(res_df$zval)), .Machine$double.xmin))
    coef_syn <- -log2(pmax(2 * stats::pnorm(-abs(res_df$zval_syn)), .Machine$double.xmin))
  } else {
    stop("stat must be either 'chisq', 'zval', 'pval', or 'beta'")
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
  object@misc$Calibre$threshold <- calibrated$threshold
  object@misc$Calibre$threshold_audit <- calibrated$audit
  object@misc$Calibre$fdrcutoff <- fdrcutoff
  object@misc$Calibre$stat <- stat
  object@misc$Calibre$n_sig <- length(sig_idx)
  object <- .calibre_sync_object(object)
  if (returnObject) object else sig_df
}
