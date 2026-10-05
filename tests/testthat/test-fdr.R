test_that("exhaustive FDR calibration handles a non-monotone empirical ratio", {
  real <- c(10, 9, 8, 7)
  null <- c(9.5, 7.5)
  fit <- binary_search(real, null, q = 0.4, return_audit = TRUE)

  expect_equal(fit$audit$threshold, c(7, 8, 9, 10))
  expect_equal(fit$audit$fdp_hat, c(2.2 / 4, 1.2 / 3, 1.2 / 2, 0.2))
  expect_equal(fit$threshold, 8)
})

test_that("calibration returns no discoveries when no threshold is valid", {
  fit <- binary_search(c(5, 4), c(6, 5), q = 0.1,
                       conservative = TRUE, return_audit = TRUE)
  expect_equal(fit$threshold, Inf)
})

test_that("invalid FDR inputs fail clearly", {
  expect_error(binary_search("x", 1, q = 0.1))
})

test_that("fast threshold counts exactly match exhaustive scans", {
  set.seed(20260818)
  real <- c(sample(c(-1, 0, 0.5, 1, 2), 200, replace = TRUE), NA, Inf)
  synthetic <- c(sample(c(-2, 0, 0.5, 1.5, 3), 180, replace = TRUE), NA, -Inf)
  fit <- binary_search(real, synthetic, q = 0.2, conservative = TRUE,
                       return_audit = TRUE)
  real_finite <- real[is.finite(real)]
  synthetic_finite <- synthetic[is.finite(synthetic)]

  expect_equal(
    fit$audit$n_real,
    vapply(fit$audit$threshold, function(x) sum(real_finite >= x), integer(1))
  )
  expect_equal(
    fit$audit$n_null,
    vapply(fit$audit$threshold, function(x) sum(synthetic_finite >= x), integer(1))
  )
})
