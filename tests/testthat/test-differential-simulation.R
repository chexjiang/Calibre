test_that("DA simulator has donor-level conditions and correct truth", {
  simulation <- simulate_da_data(
    n_donors = 6, cells_per_donor = 10, n_donors_control = 3,
    n_pairs = 20, effect_fraction = 0.25, effect_size = 0.8, seed = 4
  )
  expect_equal(dim(simulation$rna), c(20, 60))
  expect_equal(dim(simulation$atac), c(20, 60))
  expect_equal(sum(simulation$truth), 5)
  donor_conditions <- tapply(
    as.character(simulation$metadata$condition),
    simulation$metadata$donor, function(x) length(unique(x))
  )
  expect_true(all(donor_conditions == 1L))
  expect_equal(
    as.integer(table(simulation$metadata$condition)), c(30L, 30L)
  )
})

test_that("DA global null contains no true pairs", {
  simulation <- simulate_da_data(
    n_donors = 4, cells_per_donor = 6, n_donors_control = 2,
    n_pairs = 10, effect_fraction = 0.5, effect_size = 0, seed = 5
  )
  expect_false(any(simulation$truth))
  expect_true(all(simulation$beta_diff == 0))
})

test_that("DA scoring computes FDR and power", {
  score <- score_da_selection(c(1, 3, 4), c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(score$true_positives, 2)
  expect_equal(score$false_positives, 1)
  expect_equal(score$fdr, 1 / 3)
  expect_equal(score$power, 1)
})
