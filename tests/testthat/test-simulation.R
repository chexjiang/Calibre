test_that("all simulation DGPs are deterministic and produce valid counts", {
  for (dgp in c("poisson", "nb", "zinb")) {
    first <- simulate_epilink_data(80, 40, effect_size = 0.6, dgp = dgp, seed = 19)
    second <- simulate_epilink_data(80, 40, effect_size = 0.6, dgp = dgp, seed = 19)
    expect_identical(first$rna, second$rna)
    expect_identical(first$atac, second$atac)
    expect_true(all(first$rna >= 0 & first$rna == round(first$rna)))
    expect_true(all(first$atac %in% c(0, 1)))
    expect_equal(sum(first$truth), 4)
  }
})

test_that("ZINB adds structural zeros relative to the matched NB draw", {
  nb <- simulate_epilink_data(500, 100, dgp = "nb", zero_inflation = 0.5, seed = 31)
  zinb <- simulate_epilink_data(500, 100, dgp = "zinb", zero_inflation = 0.5, seed = 31)
  expect_gt(mean(zinb$rna == 0), mean(nb$rna == 0) + 0.2)
})

test_that("zero effect is an all-null simulation", {
  simulation <- simulate_epilink_data(50, 20, effect_fraction = 0.5,
                                       effect_size = 0, seed = 32)
  expect_false(any(simulation$truth))
  expect_true(all(simulation$beta == 0))
})

test_that("small simulation method comparison returns auditable metrics", {
  simulation <- simulate_epilink_data(100, 30, effect_fraction = 0.2,
                                       effect_size = 0.8, dgp = "poisson", seed = 7)
  result <- run_epilink_simulation(
    simulation,
    method_ids = c("glm_poisson_bh", "epilink_poisson_residual", "scmultimap"),
    q = 0.1, seed = 9
  )
  expect_equal(result$method_id,
               c("glm_poisson_bh", "epilink_poisson_residual", "scmultimap"))
  expect_true(all(result$fdr >= 0 & result$fdr <= 1))
  expect_true(all(result$power >= 0 & result$power <= 1))
  expect_true(all(result$n_fitted == 30))
  expect_true(all(is.na(result$error)))
  expect_true(result$finite_sample_correction[result$method_id == "epilink_poisson_residual"])
  expect_true(result$plus_one[result$method_id == "epilink_poisson_residual"])
  expect_false(any(result$plus_one[result$method_id != "epilink_poisson_residual"]))
})

test_that("Signac candidate benchmark is deterministic and complete", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  simulation <- simulate_epilink_data(
    100, 30, effect_fraction = 0.2, effect_size = 0.8, seed = 41
  )
  first <- run_epilink_simulation(simulation, method_ids = "signac", seed = 42)
  second <- run_epilink_simulation(simulation, method_ids = "signac", seed = 42)
  expect_equal(first[, c("discoveries", "fdr", "power", "n_fitted")],
               second[, c("discoveries", "fdr", "power", "n_fitted")])
  expect_equal(first$n_fitted, 30)
  expect_true(is.na(first$error))
})

test_that("simulation method registry has unique identifiers", {
  registry <- epilink_simulation_methods()
  expect_length(unique(registry$method_id), nrow(registry))
  expect_true(all(stats::na.omit(registry$family) %in% c("poisson", "nb")))
  expect_equal(registry$method_id,
               c("glm_poisson_bh", "glm_nb_bh",
                 "epilink_poisson_residual", "epilink_nb_residual",
                 "scmultimap", "signac"))
  expect_true(all(registry$null_generation[registry$correction == "EpiLink"] == "s"))
  expect_true(all(registry$finite_sample_correction[registry$correction == "EpiLink"]))
  expect_equal(registry$correction[registry$method_id == "scmultimap"], "scMultiMap")
  expect_equal(registry$correction[registry$method_id == "signac"], "Signac")
})
