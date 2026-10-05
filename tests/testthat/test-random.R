test_that("pair seeds are reproducible and pair-specific", {
  first <- epilink_pair_seeds(100, 42)
  second <- epilink_pair_seeds(100, 42)
  expect_identical(first, second)
  expect_length(unique(first), 100)
  expect_false(identical(first, epilink_pair_seeds(100, 43)))
})

test_that("local seeding does not alter the caller RNG state", {
  set.seed(91)
  before <- .Random.seed
  value <- with_epilink_seed(12, runif(3))
  expect_length(value, 3)
  expect_identical(.Random.seed, before)
})

