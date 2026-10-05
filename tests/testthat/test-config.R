test_that("configured historical inputs exist", {
  skip_if_not_installed("yaml")
  config <- yaml::read_yaml(file.path("..", "..", "config", "datasets.yml"))
  for (dataset in c("pbmc", "bmmc")) {
    paths <- unlist(config[[dataset]][c("rna", "atac", "metadata", "candidate_pairs")])
    expect_true(all(file.exists(file.path("..", "..", paths))), info = dataset)
  }
  expect_equal(config$pbmc$output_dir, "results/pbmc")
  expect_equal(config$bmmc$output_dir, "results/bmmc")
  # These dataset-specific references reproduce the incoming CRISPR count
  # tables (PBMC/result_saori used hg19; BMMC/result used hg38).
  expect_match(config$pbmc$crispr_positive, "hg38_positive[.]tsv$")
  expect_match(config$pbmc$crispr_negative, "hg38_negative[.]tsv$")
  expect_equal(config$pbmc$crispr_reference_build, "hg38")
  expect_match(config$bmmc$crispr_positive, "hg38_positive[.]tsv$")
  expect_match(config$bmmc$crispr_negative, "hg38_negative[.]tsv$")
  expect_equal(config$bmmc$crispr_reference_build, "hg38")
  expect_equal(config$defaults$statistic, "zval")
  expect_true(config$defaults$correction)
  expect_true(config$defaults$conservative)
})

test_that("simulation profiles cover Poisson, NB, and ZINB", {
  skip_if_not_installed("yaml")
  config <- yaml::read_yaml(file.path("..", "..", "config", "simulation.yml"))
  expect_equal(unlist(config$defaults$data_generating_models), c("poisson", "nb", "zinb"))
  expect_false(config$defaults$include_scent)
  expect_true(config$defaults$epilink_correction)
  expect_true(config$defaults$epilink_conservative)
  expect_true(all(c("diagnostic", "power_pilot", "trend", "sample_size", "quick", "full") %in% names(config$profiles)))
  expect_equal(unlist(config$profiles$diagnostic$effect_sizes), c(0, 0.2))
  expect_equal(config$profiles$power_pilot$effect_fraction, 0.2)
  expect_equal(length(unlist(config$profiles$trend$effect_sizes)), 10)
  expect_equal(unlist(config$profiles$sample_size$n_cells), c(100, 200, 300, 400))
})
