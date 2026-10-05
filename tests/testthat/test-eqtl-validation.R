test_that("ranked eQTL curve uses signed descending scores", {
  scored <- data.frame(
    gene = c("G1", "G2", "G3", "G4"),
    peak = paste0("chr1-", c(10, 20, 30, 40), "-", c(11, 21, 31, 41)),
    score = c(4, 3, -10, 1)
  )
  units <- data.table::data.table(
    variant_id = paste0("v", 1:4), peak = scored$peak, gene = scored$gene,
    eqtl_gene = c("G1", "X", "G3", "X"), pip = 0.9,
    positive = c(TRUE, FALSE, TRUE, FALSE)
  )
  curve <- ranked_curve(scored, units, c(0.25, 1))
  expect_equal(curve$TP[[1]], 1)
  expect_equal(curve$score_threshold[[1]], 4)
  expect_equal(curve$recall[[2]], 1)
})

test_that("SCENT odds ratio and gene bootstrap are reproducible", {
  counts <- data.frame(
    gene = c("A", "B"),
    causal_annotated = c(4, 2), causal_not_annotated = c(1, 2),
    common_annotated = c(10, 10), common_not_annotated = c(90, 90)
  )
  observed <- odds_ratio_from_counts(counts)
  expect_gt(observed$odds_ratio_ha, 1)
  one <- bootstrap_scent_odds(counts, replicates = 20, seed = 7)
  two <- bootstrap_scent_odds(counts, replicates = 20, seed = 7)
  expect_equal(one, two)
})

test_that("peak-eQTL units require the correct target gene", {
  universe <- data.frame(
    gene = c("RIGHT", "WRONG"), peak = rep("chr1-100-120", 2)
  )
  eqtl <- data.frame(
    chr = "chr1", start = 110L, end = 110L, variant_id = "v1",
    gene = "RIGHT", pip = 0.9
  )
  units <- make_ranked_validation_units(universe, eqtl, 0.5)
  expect_equal(sum(units$positive), 1)
  expect_true(units[gene == "RIGHT"]$positive)
  expect_false(units[gene == "WRONG"]$positive)
})
