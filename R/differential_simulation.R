# Utilities for the model-based differential-association simulation.

da_simulation_methods <- function() {
  data.frame(
    method_id = c(
      "epilink_poisson", "epilink_nb", "glm_poisson", "glm_nb",
      "scmultimap"
    ),
    method = c(
      "EpiLink-Poisson", "EpiLink-NB", "Poisson-GLM", "NB-GLM",
      "scMultiMap"
    ),
    stringsAsFactors = FALSE
  )
}

#' Simulate paired RNA/ATAC data with a condition-specific association
#'
#' The DGP is Poisson or negative binomial with
#'   log(mu) = intercept + condition_effect * condition +
#'             atac_effect * atac + beta_diff * atac * condition.
#' Donors determine condition assignment but have no donor-specific effect in
#' this pilot.
simulate_da_data <- function(
    n_donors = 20L, cells_per_donor = 100L, n_donors_control = 10L,
    n_pairs = 40L, effect_fraction = 0.25, effect_size = 0.8,
    atac_probability = 0.2, intercept = 1, atac_effect = 0.5,
    condition_effect = 0.2, dgp = c("poisson", "nb", "zinb"), nb_size = 2,
    zero_inflation = 0.30, seed = 1L
) {
  dgp <- match.arg(dgp)
  n_donors <- as.integer(n_donors)
  cells_per_donor <- as.integer(cells_per_donor)
  n_donors_control <- as.integer(n_donors_control)
  n_pairs <- as.integer(n_pairs)
  if (n_donors < 2L || cells_per_donor < 2L || n_pairs < 1L) {
    stop("n_donors, cells_per_donor, or n_pairs is too small", call. = FALSE)
  }
  if (n_donors_control < 1L || n_donors_control >= n_donors) {
    stop("n_donors_control must split donors into two nonempty groups",
         call. = FALSE)
  }
  if (effect_fraction < 0 || effect_fraction > 1) {
    stop("effect_fraction must be in [0, 1]", call. = FALSE)
  }
  if (atac_probability <= 0 || atac_probability >= 1) {
    stop("atac_probability must be in (0, 1)", call. = FALSE)
  }
  if (nb_size <= 0) stop("nb_size must be positive", call. = FALSE)
  if (zero_inflation < 0 || zero_inflation >= 1) {
    stop("zero_inflation must be in [0, 1)", call. = FALSE)
  }
  with_epilink_seed(as.integer(seed), {
    donor_ids <- paste0("donor", seq_len(n_donors))
    donor <- factor(
      rep(donor_ids, each = cells_per_donor), levels = donor_ids
    )
    donor_condition <- c(
      rep("Control", n_donors_control),
      rep("Case", n_donors - n_donors_control)
    )
    condition <- factor(
      rep(donor_condition, each = cells_per_donor),
      levels = c("Control", "Case")
    )
    condition01 <- as.numeric(condition) - 1
    n_cells <- length(donor)
    cell_ids <- paste0("cell", seq_len(n_cells))

    atac <- matrix(
      stats::rbinom(n_pairs * n_cells, 1, atac_probability),
      nrow = n_pairs, ncol = n_cells
    )
    truth <- rep(FALSE, n_pairs)
    n_true <- as.integer(round(n_pairs * effect_fraction))
    if (n_true > 0L && effect_size != 0) {
      truth[sample.int(n_pairs, n_true)] <- TRUE
    }
    beta_diff <- ifelse(truth, effect_size, 0)
    condition_matrix <- matrix(
      condition01, nrow = n_pairs, ncol = n_cells, byrow = TRUE
    )
    eta <- intercept + condition_effect * condition_matrix +
      atac_effect * atac + beta_diff * atac * condition_matrix
    mean_matrix <- exp(eta)
    rna <- if (dgp == "poisson") {
      matrix(stats::rpois(length(mean_matrix), lambda = c(mean_matrix)),
             nrow = n_pairs, ncol = n_cells)
    } else {
      counts <- matrix(stats::rnbinom(
        length(mean_matrix), mu = c(mean_matrix), size = nb_size
      ), nrow = n_pairs, ncol = n_cells)
      if (dgp == "zinb") {
        structural_zero <- matrix(
          stats::rbinom(length(mean_matrix), 1L, zero_inflation) == 1L,
          nrow = n_pairs, ncol = n_cells
        )
        counts[structural_zero] <- 0L
      }
      counts
    }

    gene_ids <- paste0("gene", seq_len(n_pairs))
    peak_start <- seq_len(n_pairs) * 1000L
    peak_ids <- paste0("chr1-", peak_start, "-", peak_start + 199L)
    rownames(rna) <- gene_ids
    rownames(atac) <- peak_ids
    colnames(rna) <- colnames(atac) <- cell_ids
    metadata <- data.frame(
      newCT = "all", condition = condition, donor = donor,
      row.names = cell_ids, stringsAsFactors = FALSE
    )
    pairs <- data.frame(
      gene = gene_ids, peak = peak_ids, pair_id = seq_len(n_pairs),
      stringsAsFactors = FALSE
    )
    list(
      rna = rna, atac = atac, metadata = metadata, pairs = pairs,
      truth = truth, beta_diff = beta_diff, n_cells = n_cells,
      n_pairs = n_pairs,
      parameters = list(
        n_donors = n_donors, cells_per_donor = cells_per_donor,
        n_donors_control = n_donors_control,
        effect_fraction = effect_fraction, effect_size = effect_size,
        atac_probability = atac_probability, intercept = intercept,
        atac_effect = atac_effect, condition_effect = condition_effect,
        dgp = dgp, nb_size = nb_size, zero_inflation = zero_inflation,
        seed = seed
      )
    )
  })
}

make_da_seurat <- function(simulation) {
  if (!requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("Signac", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE)) {
    stop("Seurat, Signac, and Matrix are required", call. = FALSE)
  }
  object <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(simulation$rna, sparse = TRUE),
    meta.data = simulation$metadata, assay = "RNA"
  )
  object[["ATAC"]] <- Signac::CreateChromatinAssay(
    counts = Matrix::Matrix(simulation$atac, sparse = TRUE),
    sep = c("-", "-")
  )
  object
}

score_da_selection <- function(selected, truth) {
  selected <- unique(as.integer(selected))
  selected <- selected[is.finite(selected) & selected >= 1L &
                         selected <= length(truth)]
  true_positives <- sum(truth[selected])
  false_positives <- length(selected) - true_positives
  n_true <- sum(truth)
  data.frame(
    discoveries = length(selected), true_positives = true_positives,
    false_positives = false_positives,
    fdr = if (length(selected)) false_positives / length(selected) else 0,
    power = if (n_true) true_positives / n_true else NA_real_,
    stringsAsFactors = FALSE
  )
}

da_pair_keys <- function(x) paste(x$gene, x$peak, sep = "\r")

fit_da_glm <- function(simulation, family = c("poisson", "nb"), alpha = 0.1,
                       n_cores = 1L) {
  family <- match.arg(family)
  condition <- as.numeric(simulation$metadata$condition) - 1
  fit_one <- function(index) {
    data <- data.frame(
      exprs = as.numeric(simulation$rna[index, ]),
      condition = condition,
      atac = as.numeric(simulation$atac[index, ] > 0)
    )
    fit <- try(if (family == "poisson") {
      stats::glm(
        exprs ~ condition + atac + atac:condition,
        family = stats::poisson(), data = data
      )
    } else {
      MASS::glm.nb(exprs ~ condition + atac + atac:condition, data = data)
    }, silent = TRUE)
    if (inherits(fit, "try-error")) return(c(statistic = NA, pval = NA))
    coefficients <- summary(fit)$coefficients
    term <- intersect(c("condition:atac", "atac:condition"),
                      rownames(coefficients))
    if (length(term) != 1L) return(c(statistic = NA, pval = NA))
    c(statistic = coefficients[term, 3L], pval = coefficients[term, 4L])
  }
  indices <- seq_len(simulation$n_pairs)
  fitted <- if (n_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(indices, fit_one, mc.cores = n_cores,
                       mc.set.seed = FALSE)
  } else lapply(indices, fit_one)
  fitted <- do.call(rbind, fitted)
  results <- cbind(
    simulation$pairs,
    data.frame(
      statistic = fitted[, "statistic"], pval = fitted[, "pval"],
      stringsAsFactors = FALSE
    )
  )
  results$padj <- stats::p.adjust(results$pval, method = "BH")
  list(
    results = results,
    selected = which(is.finite(results$padj) & results$padj <= alpha),
    n_fitted = sum(is.finite(results$pval))
  )
}

# The theoretical DA branch fits the full model
#   exprs ~ condition + atac + atac:condition
# and tests the interaction (the differential-association estimand).  This is
# the differential analogue of the model-based simulation's theoretical GLM
# branch; no EpiLink residual null or finite-cell correction is applied.
fit_da_glm_theory <- function(simulation, object,
                               family = c("poisson", "nb"), alpha = 0.1,
                               n_cores = 1L) {
  family <- match.arg(family)
  fitted <- EpiLinkDiff(
    object = object,
    gene_peak_df = simulation$pairs[, c("gene", "peak"), drop = FALSE],
    celltype_col = "newCT", family = family, celltype = NULL,
    condition_col = "condition", donor_col = NULL,
    donor_effect = FALSE, donor_robust = NULL,
    covariates = NULL, stat = "zval", bin = TRUE,
    selection = "rate", returnObject = FALSE, fdrcutoff = alpha,
    ncore = n_cores, theory = TRUE, correction = FALSE,
    conservative = FALSE, verbose = FALSE, seed = 1L
  )
  if (is.null(fitted)) {
    return(list(results = NULL, selected = integer(), n_fitted = 0L))
  }
  selected_keys <- da_pair_keys(fitted[
    is.finite(fitted$padj) & fitted$padj <= alpha,
    c("gene", "peak"), drop = FALSE
  ])
  selected <- which(
    da_pair_keys(simulation$pairs) %in% selected_keys
  )
  list(
    results = fitted,
    selected = selected,
    n_fitted = sum(is.finite(fitted$pval))
  )
}

fit_da_epilink <- function(simulation, object,
                           family = c("poisson", "nb"), alpha = 0.1,
                           seed = 1L, n_cores = 1L,
                           null_generation = "s", correction = TRUE,
                           conservative = TRUE, donor_robust = FALSE) {
  family <- match.arg(family)
  fitted <- EpiLinkDiff(
    object = object,
    gene_peak_df = simulation$pairs[, c("gene", "peak")],
    celltype_col = "newCT", family = family, celltype = NULL,
    condition_col = "condition",
    donor_col = if (isTRUE(donor_robust)) "donor" else NULL,
    donor_effect = FALSE,
    donor_robust = if (isTRUE(donor_robust)) TRUE else NULL,
    covariates = NULL, null_generation = null_generation,
    stat = "zval", bin = TRUE, selection = "rate",
    returnObject = TRUE, fdrcutoff = alpha, ncore = n_cores,
    theory = FALSE, conservative = conservative, correction = correction,
    verbose = FALSE, seed = seed
  )
  if (is.null(fitted)) {
    return(list(results = NULL, selected = integer(), n_fitted = 0L))
  }
  results <- fitted@misc$EpiLink$res_df
  selected_keys <- da_pair_keys(results[
    !is.na(results$significant) & results$significant,
    c("gene", "peak"), drop = FALSE
  ])
  pair_keys <- da_pair_keys(simulation$pairs)
  list(
    results = results,
    selected = which(pair_keys %in% selected_keys),
    n_fitted = sum(is.finite(results$pval))
  )
}

fit_da_scmultimap <- function(simulation, object, alpha = 0.1, seed = 1L,
                              n_permutations = 100L, irls = TRUE,
                              permutation_cores = 1L) {
  if (!requireNamespace("scMultiMap", quietly = TRUE)) {
    stop("scMultiMap is required", call. = FALSE)
  }
  metadata <- simulation$metadata
  pairs <- simulation$pairs[, c("gene", "peak"), drop = FALSE]
  condition_levels <- levels(metadata$condition)
  all_cells <- rownames(metadata)
  cell_condition <- as.character(metadata$condition)
  control_cells <- all_cells[cell_condition == condition_levels[[1L]]]
  case_cells <- all_cells[cell_condition == condition_levels[[2L]]]
  if (!length(control_cells) || !length(case_cells)) {
    stop("both conditions must contain at least one cell", call. = FALSE)
  }

  run_subset <- function(cells) {
    suppressMessages(run_scMultiMap_compat(
      object[, cells], pairs, gene_assay = "RNA", peak_assay = "ATAC",
      bsample = NULL, irls = irls, verbose = FALSE
    ))
  }
  align <- function(result) {
    pair_key <- da_pair_keys(pairs)
    result_key <- da_pair_keys(result)
    out <- result[match(pair_key, result_key), , drop = FALSE]
    out$gene <- pairs$gene
    out$peak <- pairs$peak
    out
  }
  control_result <- align(run_subset(control_cells))
  case_result <- align(run_subset(case_cells))
  observed_diff <- case_result$covar - control_result$covar

  n_permutations <- max(1L, as.integer(n_permutations))
  permutation_cells <- with_epilink_seed(
    as.integer(seed), lapply(seq_len(n_permutations), function(i) {
      sample(all_cells, length(case_cells), replace = FALSE)
    })
  )
  fit_permutation <- function(index) {
    case_null <- permutation_cells[[index]]
    control_null <- setdiff(all_cells, case_null)
    null_control <- align(run_subset(control_null))
    null_case <- align(run_subset(case_null))
    null_case$covar - null_control$covar
  }
  null_list <- if (permutation_cores > 1L &&
                   .Platform$OS.type != "windows") {
    parallel::mclapply(
      seq_len(n_permutations), fit_permutation,
      mc.cores = permutation_cores, mc.set.seed = FALSE
    )
  } else lapply(seq_len(n_permutations), fit_permutation)
  null_diff <- do.call(cbind, null_list)
  p_upper <- rowMeans(observed_diff > null_diff, na.rm = TRUE)
  p_lower <- rowMeans(observed_diff < null_diff, na.rm = TRUE)
  pval <- pmin(p_upper, p_lower) * 2
  pval[rowSums(is.finite(null_diff)) == 0L |
         !is.finite(observed_diff)] <- NA_real_
  results <- data.frame(
    simulation$pairs,
    statistic = observed_diff, pval = pval,
    padj = stats::p.adjust(pval, method = "BH"),
    stringsAsFactors = FALSE
  )
  list(
    results = results,
    selected = which(is.finite(results$padj) & results$padj <= alpha),
    n_fitted = sum(is.finite(results$pval)),
    n_permutations = n_permutations
  )
}

run_da_benchmark <- function(simulation, alpha = 0.1, seed = 1L,
                             pair_cores = 1L, null_generation = "s",
                             correction = TRUE, conservative = TRUE,
                             donor_robust = FALSE,
                             scmultimap_permutations = 100L,
                             scmultimap_irls = TRUE) {
  object <- make_da_seurat(simulation)
  registry <- da_simulation_methods()
  method_seeds <- epilink_pair_seeds(nrow(registry), as.integer(seed))
  rows <- vector("list", nrow(registry))
  for (index in seq_len(nrow(registry))) {
    method_id <- registry$method_id[[index]]
    started <- proc.time()[["elapsed"]]
    fitted <- tryCatch({
      if (method_id == "epilink_poisson") {
        fit_da_epilink(
          simulation, object, family = "poisson", alpha = alpha,
          seed = method_seeds[[index]], n_cores = pair_cores,
          null_generation = null_generation, correction = correction,
          conservative = conservative, donor_robust = donor_robust
        )
      } else if (method_id == "epilink_nb") {
        fit_da_epilink(
          simulation, object, family = "nb", alpha = alpha,
          seed = method_seeds[[index]], n_cores = pair_cores,
          null_generation = null_generation, correction = correction,
          conservative = conservative, donor_robust = donor_robust
        )
      } else if (method_id == "glm_poisson") {
        fit_da_glm_theory(simulation, object, family = "poisson",
                          alpha = alpha, n_cores = pair_cores)
      } else if (method_id == "glm_nb") {
        fit_da_glm_theory(simulation, object, family = "nb",
                          alpha = alpha, n_cores = pair_cores)
      } else {
        fit_da_scmultimap(
          simulation, object, alpha = alpha, seed = method_seeds[[index]],
          n_permutations = scmultimap_permutations,
          irls = scmultimap_irls, permutation_cores = pair_cores
        )
      }
    }, error = function(e) structure(
      list(selected = integer(), n_fitted = 0L),
      error = conditionMessage(e)
    ))
    score <- score_da_selection(fitted$selected, simulation$truth)
    score$method <- registry$method[[index]]
    score$method_id <- method_id
    score$n_fitted <- fitted$n_fitted
    score$runtime_seconds <- proc.time()[["elapsed"]] - started
    score$error <- if (is.null(attr(fitted, "error"))) NA_character_ else
      attr(fitted, "error")
    rows[[index]] <- score
  }
  do.call(rbind, rows)
}

summarize_da_metrics <- function(metrics) {
  long <- rbind(
    data.frame(metrics[, c("dgp", "n_cells", "effect_size", "method", "method_id")],
               metric = "FDR", value = metrics$fdr),
    data.frame(metrics[, c("dgp", "n_cells", "effect_size", "method", "method_id")],
               metric = "Power", value = metrics$power)
  )
  grouping <- interaction(
    long$dgp, long$n_cells, long$effect_size,
    long$method_id, long$metric, drop = TRUE
  )
  rows <- lapply(split(long, grouping), function(x) {
    values <- x$value[is.finite(x$value)]
    data.frame(
      dgp = x$dgp[[1]], n_cells = x$n_cells[[1]],
      effect_size = x$effect_size[[1]],
      method = x$method[[1]], method_id = x$method_id[[1]],
      metric = x$metric[[1]],
      mean = if (length(values)) mean(values) else NA_real_,
      se = if (length(values) > 1L) stats::sd(values) / sqrt(length(values))
      else NA_real_,
      n = length(values), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$dgp, out$n_cells, out$effect_size,
            out$method_id, out$metric), , drop = FALSE]
}
