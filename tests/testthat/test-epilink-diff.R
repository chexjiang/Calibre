make_diff_data <- function(n_donor = 8L, cells_per_donor = 40L, seed = 13L) {
  set.seed(seed)
  donor <- factor(rep(paste0("d", seq_len(n_donor)), each = cells_per_donor))
  condition <- rep(rep(0:1, each = n_donor / 2L), each = cells_per_donor)
  atac_raw <- rbinom(length(donor), 1, 0.4)
  atac <- atac_raw - ave(atac_raw, donor)
  donor_effect <- rep(rnorm(n_donor, sd = 0.2), each = cells_per_donor)
  mu <- exp(0.4 + donor_effect + 0.5 * atac + 0.8 * atac * condition)
  data.frame(
    exprs = rpois(length(mu), mu),
    atac = atac,
    condition_id = condition,
    atac_condition = atac * condition,
    individual_id = donor
  )
}

test_that("differential pair fit reports the interaction effect", {
  result <- EpiLinkDiff_single(
    make_diff_data(), family = "poisson", null_generation = "p",
    individual = TRUE, seed = 22
  )
  expect_true(all(c(
    "beta_control", "beta_diff", "beta_disease", "se_diff",
    "wald_z", "wald_p", "beta_syn", "zval_syn"
  ) %in% names(result)))
  expect_equal(result$beta, result$beta_diff)
  expect_equal(result$zval, result$wald_z)
  expect_equal(result$beta_disease, result$beta_control + result$beta_diff)
  expect_gt(result$beta_diff, 0)
})

test_that("differential null generation is reproducible and keeps a common slope", {
  data <- make_diff_data()
  first <- EpiLinkDiff_single(
    data, null_generation = "p", individual = TRUE, seed = 91
  )
  second <- EpiLinkDiff_single(
    data, null_generation = "p", individual = TRUE, seed = 91
  )
  expect_equal(first, second)
  expect_error(
    EpiLinkDiff_single(data, null_generation = "n", individual = TRUE),
    "does not preserve"
  )
})

test_that("cluster robust variance can be used without individual fixed effects", {
  data <- make_diff_data()
  result <- EpiLinkDiff_single(
    data, null_generation = "p", individual = FALSE,
    cluster = data$individual_id, seed = 92
  )
  expect_true(is.finite(result$se_diff))
  expect_true(is.finite(result$zval))
  expect_true(is.finite(result$zval_syn))
})

test_that("EpiLinkDiff returns labeled differential results", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  set.seed(42)
  n_donor <- 8L
  cells_per_donor <- 30L
  donor <- factor(rep(paste0("d", seq_len(n_donor)), each = cells_per_donor))
  condition <- factor(
    rep(rep(c("CTRL", "AD"), each = n_donor / 2L), each = cells_per_donor),
    levels = c("CTRL", "AD")
  )
  condition01 <- as.numeric(condition) - 1
  peak <- rbinom(length(donor), 1, 0.4)
  centered_peak <- peak - ave(peak, donor)
  expression <- rpois(
    length(donor), exp(0.5 + 0.4 * centered_peak + 0.7 * centered_peak * condition01)
  )
  cells <- paste0("cell", seq_along(donor))
  rna <- Matrix::Matrix(rbind(
    g1 = expression,
    g2 = rpois(length(donor), 1)
  ), sparse = TRUE)
  atac <- Matrix::Matrix(rbind(
    `chr1-100-200` = peak,
    `chr1-300-400` = rbinom(length(donor), 1, 0.3)
  ), sparse = TRUE)
  colnames(rna) <- colnames(atac) <- cells
  metadata <- data.frame(
    newCT = "test", condition = condition, donor = donor, row.names = cells
  )
  object <- Seurat::CreateSeuratObject(rna, meta.data = metadata, assay = "RNA")
  object[["ATAC"]] <- Signac::CreateChromatinAssay(atac)
  pairs <- data.frame(gene = "g1", peak = "chr1-100-200")
  expect_error(
    EpiLinkDiff(object, pairs, condition_col = "condition", donor_effect = TRUE),
    "donor_col is required"
  )

  fitted <- EpiLinkDiff(
    object, pairs, condition_col = "condition", donor_col = "donor",
    donor_effect = TRUE,
    celltype = "test", null_generation = "p", ncore = 1,
    fdrcutoff = 0.2, returnObject = TRUE, verbose = FALSE, seed = 7
  )
  result <- fitted@misc$EpiLink$res_df
  expect_identical(fitted@misc$EpiLink$analysis, "differential")
  expect_identical(fitted@misc$EpiLink$condition_levels, c("CTRL", "AD"))
  expect_identical(fitted@misc$EpiLink$donor_col, "donor")
  expect_true(fitted@misc$EpiLink$donor_effect)
  expect_true(fitted@misc$EpiLink$donor_robust)
  expect_identical(result$reference_condition, "CTRL")
  expect_identical(result$comparison_condition, "AD")
  expect_true(all(c("calibrated_q", "significant") %in% names(result)))
})
