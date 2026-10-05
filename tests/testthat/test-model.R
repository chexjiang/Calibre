make_pair_data <- function(n = 300L, seed = 1L) {
  set.seed(seed)
  atac <- rbinom(n, 1, 0.3)
  batch <- factor(sample(c("a", "b"), n, replace = TRUE))
  depth <- runif(n, 0.8, 1.2)
  mu <- exp(0.2 + 0.7 * atac + 0.25 * (batch == "b") + log(depth))
  data.frame(exprs = rpois(n, mu), atac = atac, batch = batch, log_depth = log(depth))
}

test_that("pair fit is reproducible for a fixed seed", {
  data <- make_pair_data()
  first <- EpiLink_single(data, c("batch", "log_depth"),
                          null_generation = "s", seed = 101)
  second <- EpiLink_single(data, c("batch", "log_depth"),
                           null_generation = "s", seed = 101)
  expect_equal(first, second)
  expect_true(all(c("beta", "zval", "beta_syn", "zval_syn") %in% names(first)))
})

test_that("different pair seeds produce different synthetic statistics", {
  data <- make_pair_data()
  first <- EpiLink_single(data, c("batch", "log_depth"),
                          null_generation = "s", seed = 101)
  second <- EpiLink_single(data, c("batch", "log_depth"),
                           null_generation = "s", seed = 102)
  expect_false(isTRUE(all.equal(first$zval_syn, second$zval_syn)))
})

test_that("theory mode returns only the observed model", {
  result <- EpiLink_single(make_pair_data(), c("batch", "log_depth"),
                           theory = TRUE)
  expect_true(all(c("beta", "zval", "pval", "disp") %in% names(result)))
  expect_false("zval_syn" %in% names(result))
})
