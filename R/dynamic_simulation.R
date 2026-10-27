# Dynamic Calibre simulation utilities.
#
# The accessibility process is deliberately independent of pseudotime:
#   A_ip ~ Bernoulli(pi),  pi does not depend on T_i.
# The RNA process contains a common gene trajectory and an accessibility effect
# that may change over pseudotime:
#   log(mu_ip) = alpha + g(T_i) + A_ip { beta + delta h(T_i) }.

simulate_calibre_dynamic_data <- function(
    n_cells, n_pairs, effect_fraction = 0.2, effect_size = 0.5,
    dynamic_shape = "linear", dgp = c("poisson", "nb", "zinb"),
    atac_probability = 0.2, intercept = 1, atac_effect = 0.5,
    expression_trend = 1, nb_size = 10, zero_inflation = 0.30,
    seed = 1L) {
  dgp <- match.arg(dgp)
  dynamic_shape <- match.arg(
    dynamic_shape, c("linear", "early", "late", "reversal")
  )
  n_cells <- as.integer(n_cells)
  n_pairs <- as.integer(n_pairs)
  if (n_cells < 20L || n_pairs < 1L) {
    stop("n_cells must be at least 20 and n_pairs must be positive", call. = FALSE)
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

  with_calibre_seed(as.integer(seed), {
    pseudotime <- stats::runif(n_cells)

    # IMPORTANT: this is the only accessibility DGP. There is no
    # pseudotime term in the accessibility probability.
    atac <- matrix(
      stats::rbinom(n_pairs * n_cells, size = 1, prob = atac_probability),
      nrow = n_pairs, ncol = n_cells
    )

    shape <- switch(
      dynamic_shape,
      linear = pseudotime - 0.5,
      early = exp(-0.5 * ((pseudotime - 0.25) / 0.12)^2),
      late = exp(-0.5 * ((pseudotime - 0.75) / 0.12)^2),
      reversal = tanh((pseudotime - 0.5) / 0.12)
    )
    # The fitted model anchors the by-smooth at mean pseudotime, so center the
    # truth at the same point. Normalize the dynamic component to make
    # effect_size comparable across shapes.
    shape <- shape - approx(
      x = pseudotime[order(pseudotime)],
      y = shape[order(pseudotime)],
      xout = mean(pseudotime), ties = "ordered"
    )$y
    shape <- shape / max(abs(shape))
    expression_curve <- expression_trend * (
      0.65 * sin(2 * pi * pseudotime) +
        0.25 * cos(4 * pi * pseudotime)
    )

    n_true <- as.integer(round(n_pairs * effect_fraction))
    dynamic_truth <- rep(FALSE, n_pairs)
    if (n_true > 0L && effect_size > 0) {
      dynamic_truth[sample.int(n_pairs, n_true)] <- TRUE
    }
    dynamic_amplitude <- ifelse(dynamic_truth, effect_size, 0)
    association <- atac_effect + outer(dynamic_amplitude, shape)
    mu <- exp(
      matrix(intercept + expression_curve, nrow = n_pairs,
             ncol = n_cells, byrow = TRUE) + atac * association
    )
    rna <- if (dgp == "poisson") {
      matrix(stats::rpois(length(mu), lambda = c(mu)), nrow = n_pairs)
    } else {
      counts <- matrix(
        stats::rnbinom(length(mu), mu = c(mu), size = nb_size),
        nrow = n_pairs
      )
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
    rownames(atac) <- paste0(
      peak_chromosome, "-", peak_start, "-", peak_start + 199L
    )
    colnames(rna) <- colnames(atac) <- paste0("cell", seq_len(n_cells))
    list(
      rna = rna,
      atac = atac,
      pseudotime = pseudotime,
      truth = dynamic_truth,
      dynamic_truth = dynamic_truth,
      dynamic_amplitude = dynamic_amplitude,
      dynamic_shape = dynamic_shape,
      shape = shape,
      beta = rep(atac_effect, n_pairs),
      dgp = dgp,
      n_cells = n_cells,
      n_pairs = n_pairs,
      parameters = list(
        effect_fraction = effect_fraction,
        effect_size = effect_size,
        dynamic_shape = dynamic_shape,
        atac_probability = atac_probability,
        accessibility_independent_of_pseudotime = TRUE,
        intercept = intercept,
        atac_effect = atac_effect,
        expression_trend = expression_trend,
        nb_size = nb_size,
        zero_inflation = zero_inflation,
        seed = seed
      )
    )
  })
}

dynamic_simulation_methods <- function() {
  data.frame(
    method = c(
      "Calibre-dynamic Poisson",
      "Calibre-dynamic NB",
      "Dynamic GAM Poisson",
      "Dynamic GAM NB"
    ),
    method_id = c(
      "calibre_dynamic_poisson",
      "calibre_dynamic_nb",
      "gam_dynamic_poisson",
      "gam_dynamic_nb"
    ),
    family = c("poisson", "nb", "poisson", "nb"),
    nullstrap = c(TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

make_dynamic_simulation_seurat <- function(simulation) {
  object <- make_simulation_seurat(simulation)
  object@meta.data$pseudotime <- as.numeric(simulation$pseudotime)
  object
}

fit_dynamic_simulated_pairs <- function(
    simulation, object, method_spec, q = 0.1, seed = 1L,
    pseudotime_formula = "s(pseudotime, k = 5, bs = 'tp')",
    null_generation = "s", correction = TRUE, conservative = TRUE,
    stat = "chi_sq", ncore = 1L) {
  # Simulation pairs are already a fixed one-to-one candidate universe. Fit
  # the same pair-level model used by CalibreDynamic directly, avoiding the
  # extra Seurat/filtering layer used for real data.
  pair_seeds <- calibre_pair_seeds(simulation$n_pairs, seed)
  fit_one <- function(j) {
    df <- data.frame(
      exprs = as.numeric(simulation$rna[j, ]),
      atac = as.numeric(simulation$atac[j, ]),
      pseudotime = as.numeric(simulation$pseudotime)
    )
    result <- try(CalibreDynamic_single(
      df,
      pseudotime_formula = pseudotime_formula,
      family = method_spec$family,
      null_generation = null_generation,
      theory = !isTRUE(method_spec$nullstrap),
      seed = pair_seeds[j]
    ), silent = TRUE)
    if (inherits(result, "try-error") || is.null(result)) return(NULL)
    result$pair_index <- j
    result$gene <- rownames(simulation$rna)[j]
    result$peak <- rownames(simulation$atac)[j]
    result
  }
  results_list <- if (ncore > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(simulation$n_pairs), fit_one,
                       mc.cores = ncore, mc.set.seed = FALSE)
  } else lapply(seq_len(simulation$n_pairs), fit_one)
  results <- do.call(rbind, Filter(Negate(is.null), results_list))
  if (is.null(results) || !nrow(results)) {
    return(list(
      selected = integer(), selected_gam = integer(), n_fitted = 0L
    ))
  }
  gam_padj <- stats::p.adjust(results$pval, method = "BH")
  selected_gam <- results$pair_index[
    is.finite(gam_padj) & gam_padj <= q
  ]
  if (!isTRUE(method_spec$nullstrap)) {
    return(list(
      selected = selected_gam, selected_gam = selected_gam,
      n_fitted = nrow(results)
    ))
  }
  if (!stat %in% c("chi_sq", "zval")) {
    stop("stat must be 'chi_sq' or 'zval'", call. = FALSE)
  }
  real_stat <- if (stat == "chi_sq") {
    results$chi_sq
  } else {
    abs(results$zval)
  }
  null_stat <- if (stat == "chi_sq") {
    results$chi_sq_syn
  } else {
    abs(results$zval_syn)
  }
  correction_factor <- if (correction) {
    1 / (1 + sqrt(log2(nrow(results)) / simulation$n_cells))
  } else 1
  threshold <- binary_search(
    real_stat, null_stat, q = q * correction_factor,
    conservative = conservative
  )
  list(
    selected = results$pair_index[
      is.finite(real_stat) & real_stat >= threshold
    ],
    selected_gam = selected_gam,
    n_fitted = nrow(results)
  )
}

score_dynamic_selection <- function(selected, truth) {
  selected <- unique(as.integer(selected))
  selected <- selected[
    is.finite(selected) & selected >= 1L & selected <= length(truth)
  ]
  discoveries <- length(selected)
  true_positives <- if (discoveries) sum(truth[selected]) else 0L
  false_positives <- discoveries - true_positives
  n_true <- sum(truth)
  data.frame(
    discoveries = discoveries,
    true_positives = true_positives,
    false_positives = false_positives,
    fdr = false_positives / max(discoveries, 1L),
    power = if (n_true) true_positives / n_true else NA_real_,
    stringsAsFactors = FALSE
  )
}

summarize_dynamic_simulation <- function(metrics) {
  if (is.null(metrics) || !nrow(metrics)) return(data.frame())
  long <- rbind(
    transform(metrics, metric = "FDR", value = fdr),
    transform(metrics, metric = "Power", value = power)
  )
  group <- interaction(
    long$dgp, long$n_cells, long$dynamic_shape, long$effect_size,
    long$method_id, long$metric, drop = TRUE
  )
  rows <- lapply(split(long, group), function(x) {
    values <- x$value[is.finite(x$value)]
    data.frame(
      dgp = x$dgp[[1L]], n_cells = x$n_cells[[1L]],
      dynamic_shape = x$dynamic_shape[[1L]],
      effect_size = x$effect_size[[1L]],
      method = x$method[[1L]], method_id = x$method_id[[1L]],
      metric = x$metric[[1L]],
      mean = if (length(values)) mean(values) else NA_real_,
      se = if (length(values) > 1L) stats::sd(values) / sqrt(length(values)) else NA_real_,
      n = length(values), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
