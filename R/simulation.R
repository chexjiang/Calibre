#' Simulate paired ATAC and RNA measurements
#'
#' @param n_cells Number of cells.
#' @param n_pairs Number of one-to-one enhancer-gene pairs.
#' @param effect_fraction Fraction of pairs with a nonzero coefficient.
#' @param effect_size Nonzero log-fold-change coefficient.
#' @param dgp RNA count distribution: Poisson, NB, or ZINB.
#' @param atac_probability Bernoulli accessibility probability.
#' @param intercept RNA log-mean intercept.
#' @param nb_size NB size (inverse overdispersion) parameter.
#' @param zero_inflation Structural-zero probability for ZINB.
#' @param seed Dataset seed.
#' @return A list containing count matrices, coefficients, and truth labels.
#' @export
simulate_calibre_data <- function(n_cells, n_pairs, effect_fraction = 0.1,
                                   effect_size = 0.5,
                                   dgp = c("poisson", "nb", "zinb"),
                                   atac_probability = 0.2, intercept = 1,
                                   nb_size = 2, zero_inflation = 0.3,
                                   seed = 1L) {
  dgp <- match.arg(dgp)
  n_cells <- as.integer(n_cells)
  n_pairs <- as.integer(n_pairs)
  if (n_cells < 2L || n_pairs < 1L) stop("n_cells and n_pairs are too small", call. = FALSE)
  if (effect_fraction < 0 || effect_fraction > 1) stop("effect_fraction must be in [0, 1]", call. = FALSE)
  if (atac_probability <= 0 || atac_probability >= 1) stop("atac_probability must be in (0, 1)", call. = FALSE)
  if (nb_size <= 0) stop("nb_size must be positive", call. = FALSE)
  if (zero_inflation < 0 || zero_inflation >= 1) stop("zero_inflation must be in [0, 1)", call. = FALSE)

  with_calibre_seed(as.integer(seed), {
    atac <- matrix(
      stats::rbinom(n_pairs * n_cells, size = 1, prob = atac_probability),
      nrow = n_pairs, ncol = n_cells
    )
    n_true <- as.integer(round(n_pairs * effect_fraction))
    truth <- rep(FALSE, n_pairs)
    if (n_true > 0L && effect_size != 0) truth[sample.int(n_pairs, n_true)] <- TRUE
    beta <- ifelse(truth, effect_size, 0)
    mu <- exp(intercept + beta * atac)
    rna <- if (dgp == "poisson") {
      matrix(stats::rpois(length(mu), lambda = c(mu)), nrow = n_pairs)
    } else {
      counts <- matrix(stats::rnbinom(length(mu), mu = c(mu), size = nb_size), nrow = n_pairs)
      if (dgp == "zinb") {
        structural_zero <- matrix(
          stats::rbinom(length(mu), size = 1, prob = zero_inflation) == 1,
          nrow = n_pairs
        )
        counts[structural_zero] <- 0
      }
      counts
    }
    rownames(rna) <- paste0("gene", seq_len(n_pairs))
    peak_index <- seq_len(n_pairs)
    peak_start <- peak_index * 1000L
    peak_chromosome <- paste0("chr", ((peak_index - 1L) %% 5L) + 1L)
    rownames(atac) <- paste0(peak_chromosome, "-", peak_start, "-", peak_start + 199L)
    colnames(rna) <- colnames(atac) <- paste0("cell", seq_len(n_cells))
    list(
      rna = rna, atac = atac, truth = truth, beta = beta,
      dgp = dgp, n_cells = n_cells, n_pairs = n_pairs,
      parameters = list(
        effect_fraction = effect_fraction, effect_size = effect_size,
        atac_probability = atac_probability, intercept = intercept,
        nb_size = nb_size, zero_inflation = zero_inflation, seed = seed
      )
    )
  })
}

#' Available simulation method specifications
#' @return A data frame defining analysis family, correction, and null generator.
#' @export
calibre_simulation_methods <- function() {
  data.frame(
    method = c(
      "Poisson-GLM", "NB-GLM",
      "Calibre-Poisson", "Calibre-NB", "scMultiMap", "Signac", "ArchR"
    ),
    method_id = c(
      "glm_poisson_bh", "glm_nb_bh",
      "calibre_poisson_residual", "calibre_nb_residual", "scmultimap",
      "signac", "archr"
    ),
    family = c("poisson", "nb", "poisson", "nb", NA, NA, NA),
    correction = c("BH", "BH", "Calibre", "Calibre", "scMultiMap", "Signac", "ArchR"),
    finite_sample_correction = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
    null_generation = c(NA, NA, "s", "s", NA, NA, NA),
    stringsAsFactors = FALSE
  )
}

make_simulation_seurat <- function(simulation) {
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Signac", quietly = TRUE)) {
    stop("Seurat and Signac are required for simulation method evaluation")
  }
  metadata <- data.frame(newCT = rep("all", simulation$n_cells),
                         row.names = colnames(simulation$rna))
  object <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(simulation$rna, sparse = TRUE),
    meta.data = metadata, assay = "RNA"
  )
  object[["ATAC"]] <- Signac::CreateChromatinAssay(
    counts = Matrix::Matrix(simulation$atac, sparse = TRUE)
  )
  object
}

fit_simulated_pairs <- function(simulation, object, method_spec, q, seed,
                                calibre_correction = TRUE,
                                calibre_conservative = TRUE,
                                archr_cores = 1L,
                                signac_pvalue_correction = TRUE) {
  pairs <- data.frame(gene = rownames(simulation$rna),
                      peak = rownames(simulation$atac), stringsAsFactors = FALSE)
  if (identical(method_spec$correction, "scMultiMap")) {
    results <- suppressMessages(run_scMultiMap_compat(
      object, pairs, gene_assay = "RNA", peak_assay = "ATAC", irls = FALSE
    ))
    valid <- is.finite(results$pval)
    selected_pairs <- results[valid & results$padj <= q, c("gene", "peak"), drop = FALSE]
    selected_key <- paste(selected_pairs$gene, selected_pairs$peak, sep = "\r")
    pair_key <- paste(pairs$gene, pairs$peak, sep = "\r")
    return(list(selected = match(selected_key, pair_key), n_fitted = sum(valid)))
  }
  if (identical(method_spec$correction, "Signac")) {
    results <- with_calibre_seed(seed, suppressWarnings(
      run_signac_candidate_benchmark(
        object, pairs, q = q,
        correct_pvalues = signac_pvalue_correction
      )
    ))
    valid <- is.finite(results$pval)
    return(list(selected = which(valid & results$padj <= q),
                n_fitted = sum(valid)))
  }
  if (identical(method_spec$correction, "ArchR")) {
    score <- run_archr_simulation(
      simulation, q = q, seed = seed, n_cores = archr_cores
    )
    selected <- attr(score, "selected")
    if (is.null(selected)) {
      # run_archr_simulation returns a scored row; this branch is retained for
      # direct method-registry use, while the standalone runner stores the row
      # with its full metadata.
      return(list(selected = integer(), n_fitted = score$n_fitted[[1]]))
    }
    return(list(selected = selected, n_fitted = score$n_fitted[[1]]))
  }
  theory <- identical(method_spec$correction, "BH")
  fitted <- suppressWarnings(Calibre(
    object = object, gene_peak_df = pairs, family = method_spec$family,
    celltype = "all", celltype_col = "newCT", covariates = NULL,
    null_generation = "s", stat = "zval", bin = TRUE, selection = "rate",
    returnObject = !theory, fdrcutoff = q, ncore = 1,
    theory = theory, correction = if (theory) FALSE else calibre_correction,
    conservative = if (theory) FALSE else calibre_conservative,
    verbose = FALSE, seed = seed
  ))
  if (is.null(fitted)) return(list(selected = integer(), n_fitted = 0L))
  if (theory) {
    results <- fitted
    selected <- results$pair_index[results$padj <= q]
  } else {
    results <- fitted@misc$Calibre$res_df
    significant <- CalibreFilter(fitted, stat = "zval", fdrcutoff = q,
                                 correction = calibre_correction,
                                 conservative = calibre_conservative,
                                 returnObject = FALSE, verbose = FALSE)
    selected <- significant$pair_index
  }
  list(selected = selected, n_fitted = nrow(results))
}

score_simulation_selection <- function(selected, truth) {
  selected <- unique(as.integer(selected))
  selected <- selected[is.finite(selected) & selected >= 1L & selected <= length(truth)]
  discoveries <- length(selected)
  true_positives <- if (discoveries) sum(truth[selected]) else 0L
  false_positives <- discoveries - true_positives
  n_true <- sum(truth)
  data.frame(
    discoveries = discoveries,
    true_positives = true_positives,
    false_positives = false_positives,
    fdr = false_positives / max(discoveries, 1L),
    power = if (n_true) true_positives / n_true else NA_real_
  )
}

#' Run selected methods on one simulated dataset
#'
#' @param simulation Output from [simulate_calibre_data()].
#' @param method_ids Method identifiers from [calibre_simulation_methods()].
#' @param q Target FDR.
#' @param seed Master inference seed.
#' @param include_scent Run the optional, very slow SCENT adapter.
#' @param scent_cores Cores passed to SCENT.
#' @param calibre_correction Apply finite-cell correction to both Calibre models.
#' @param calibre_conservative Use the plus-one finite-sample numerator for Calibre.
#' @param signac_pvalue_correction Apply Signac's two-sided p-value multiplier.
#' @return One row per method with FDR, power, fit count, and runtime.
#' @export
run_calibre_simulation <- function(simulation, method_ids = NULL, q = 0.1,
                                    seed = 1L, include_scent = FALSE,
                                    scent_cores = 1L,
                                    calibre_correction = TRUE,
                                    calibre_conservative = TRUE,
                                    archr_cores = 1L,
                                    signac_pvalue_correction = TRUE) {
  registry <- calibre_simulation_methods()
  if (is.null(method_ids)) {
    method_ids <- registry$method_id
  } else {
    method_ids <- calibre_normalize_method_ids(method_ids)
  }
  unknown <- setdiff(method_ids, registry$method_id)
  if (length(unknown)) stop("unknown simulation methods: ", paste(unknown, collapse = ", "), call. = FALSE)
  registry <- registry[match(method_ids, registry$method_id), , drop = FALSE]
  method_seeds <- calibre_pair_seeds(nrow(registry) + as.integer(include_scent), seed)
  object <- make_simulation_seurat(simulation)

  rows <- lapply(seq_len(nrow(registry)), function(index) {
    started <- proc.time()[["elapsed"]]
    fitted <- tryCatch(
      fit_simulated_pairs(simulation, object, registry[index, ], q,
                          method_seeds[index], calibre_correction,
                          calibre_conservative, archr_cores,
                          signac_pvalue_correction),
      error = function(e) structure(list(selected = integer(), n_fitted = 0L), error = conditionMessage(e))
    )
    score <- score_simulation_selection(fitted$selected, simulation$truth)
    score$method <- registry$method[index]
    score$method_id <- registry$method_id[index]
    is_calibre <- identical(registry$correction[index], "Calibre")
    score$finite_sample_correction <- is_calibre && calibre_correction
    score$plus_one <- is_calibre && calibre_conservative
    score$n_fitted <- fitted$n_fitted
    score$runtime_seconds <- proc.time()[["elapsed"]] - started
    score$error <- if (is.null(attr(fitted, "error"))) NA_character_ else attr(fitted, "error")
    score
  })
  if (include_scent) {
    rows[[length(rows) + 1L]] <- run_scent_simulation(
      simulation, q = q, seed = method_seeds[length(method_seeds)], n_cores = scent_cores
    )
  }
  do.call(rbind, rows)
}

#' Run SCENT on one simulation dataset
#'
#' This adapter is complete but deliberately excluded from default runs because
#' SCENT performs adaptive bootstrapping up to tens of thousands of resamples.
#' @param simulation Output from [simulate_calibre_data()].
#' @param q Target FDR.
#' @param seed Random seed.
#' @param n_cores Cores passed to SCENT.
#' @return A one-row simulation score.
#' @export
run_scent_simulation <- function(simulation, q = 0.1, seed = 1L, n_cores = 1L) {
  if (!requireNamespace("SCENT", quietly = TRUE)) stop("SCENT is not installed", call. = FALSE)
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is not installed", call. = FALSE)
  metadata <- data.frame(
    celltypes = rep("all", simulation$n_cells),
    covariate = rep(1, simulation$n_cells),
    cell = colnames(simulation$rna), row.names = colnames(simulation$rna)
  )
  pairs <- data.frame(
    gene = rownames(simulation$rna), peak = rownames(simulation$atac),
    stringsAsFactors = FALSE
  )
  object <- SCENT::CreateSCENTObj(
    rna = Matrix::Matrix(simulation$rna, sparse = TRUE),
    atac = Matrix::Matrix(simulation$atac, sparse = TRUE),
    meta.data = metadata, peak.info = pairs,
    covariates = "covariate", celltypes = "celltypes"
  )
  started <- proc.time()[["elapsed"]]
  object <- with_calibre_seed(seed, SCENT::SCENT_algorithm(
    object, celltype = "all", ncores = n_cores, regr = "poisson", bin = TRUE
  ))
  results <- object@SCENT.result
  if (ncol(results) == 7L) {
    names(results) <- c("gene", "peak", "beta", "se", "z", "p", "boot_p")
  }
  if (!"boot_p" %in% names(results)) stop("SCENT result has no boot_p column", call. = FALSE)
  indices <- match(results$gene, rownames(simulation$rna))
  selected <- indices[stats::p.adjust(results$boot_p, method = "BH") <= q]
  score <- score_simulation_selection(selected, simulation$truth)
  score$method <- "SCENT"
  score$method_id <- "scent"
  score$n_fitted <- nrow(results)
  score$runtime_seconds <- proc.time()[["elapsed"]] - started
  score$error <- NA_character_
  score
}
