test_that("CRISPR validation requires peak overlap and matching gene", {
  universe <- data.frame(
    gene = c("GENE1", "GENE2", "GENE1", "GENE3"),
    peak = c("chr1-100-120", "chr1-100-120", "chr1-300-320", "chr1-500-520")
  )
  positive <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(105, 305), c(110, 310)),
    gene = c("GENE1", "GENE2")
  )
  negative <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(505, 510), gene = "GENE3"
  )

  result <- evaluate_binary_link_validation(
    universe = universe,
    selected = universe[c(1, 2), ],
    positive_ranges = positive,
    negative_ranges = negative,
    method = "test",
    cutoff = 0.1
  )

  expect_equal(result$TP, 1)
  expect_equal(result$FP, 0)
  expect_equal(result$FN, 0)
  expect_equal(result$TN, 1)
  expect_equal(result$n_validation_positive, 1)
  expect_equal(result$n_validation_negative, 1)
})
