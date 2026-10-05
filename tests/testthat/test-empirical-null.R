test_that("empirical ATAC permutation is deterministic and shared by all peaks", {
  atac <- matrix(seq_len(24), nrow = 4,
                 dimnames = list(paste0("peak", 1:4), paste0("cell", 1:6)))
  first <- permute_empirical_atac(atac, seed = 17)
  second <- permute_empirical_atac(atac, seed = 17)
  expect_identical(first$counts, second$counts)
  expect_identical(first$permutation, second$permutation)
  expect_identical(colnames(first$counts), colnames(atac))
  expect_setequal(first$permutation$source_cell, colnames(atac))
  expect_identical(
    unname(first$counts[, 1]),
    unname(atac[, first$permutation$source_index[1]])
  )
})

test_that("empirical result alignment distinguishes raw and operational calls", {
  pairs <- data.frame(
    gene = c("g1", "g2"), peak = c("p1", "p2"), pair_id = 1:2
  )
  results <- data.frame(
    gene = c("g2", "g1"), peak = c("p2", "p1"),
    pval = c(0.01, 0.04), padj = c(0.02, 0.08)
  )
  aligned <- align_empirical_result(
    results, pairs, "test", "Test", alpha = 0.05,
    operational_keys = empirical_pair_key(results[1, ])
  )
  expect_equal(aligned$pval, c(0.04, 0.01))
  expect_true(all(aligned$raw_rejected))
  expect_equal(aligned$operational_rejected, c(FALSE, TRUE))
})

test_that("1000 empirical-null replicates map to ten stable folders", {
  expect_equal(empirical_null_folder(c(1)), "01")
  expect_equal(empirical_null_folder(c(100)), "01")
  expect_equal(empirical_null_folder(c(101)), "02")
  expect_equal(empirical_null_folder(c(1000)), "10")
  expect_match(
    empirical_null_result_path("out", 1000, "signac"),
    "out/10/replicate_1000_signac[.]tsv[.]gz$"
  )
})

test_that("Poisson and NB method jobs include their paired GLM", {
  expect_identical(
    empirical_null_method_group("epilink_poisson"),
    c("epilink_poisson", "glm_poisson")
  )
  expect_identical(
    empirical_null_method_group("glm_nb"),
    c("epilink_nb", "glm_nb")
  )
  expect_identical(empirical_null_method_group("signac"), "signac")
})
