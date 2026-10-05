test_that("small Seurat analysis is reproducible across worker counts", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  set.seed(5)
  n <- 120L
  atac_signal <- rbinom(n, 1, 0.35)
  rna <- Matrix::Matrix(rbind(
    g1 = rpois(n, exp(0.5 + 0.8 * atac_signal)),
    g2 = rpois(n, 2)
  ), sparse = TRUE)
  atac <- Matrix::Matrix(rbind(
    `chr1-100-200` = atac_signal,
    `chr1-300-400` = rbinom(n, 1, 0.4)
  ), sparse = TRUE)
  colnames(rna) <- colnames(atac) <- paste0("cell", seq_len(n))
  metadata <- data.frame(newCT = rep("test", n), row.names = colnames(rna))
  object <- Seurat::CreateSeuratObject(rna, meta.data = metadata, assay = "RNA")
  object[["ATAC"]] <- Signac::CreateChromatinAssay(atac)
  pairs <- data.frame(
    gene = c("g1", "g2"),
    peak = c("chr1-100-200", "chr1-300-400")
  )

  serial <- EpiLink(object, pairs, celltype = "test", fdrcutoff = 0.2,
                    correction = TRUE, seed = 44, ncore = 1,
                    returnObject = TRUE, verbose = FALSE)
  parallel <- EpiLink(object, pairs, celltype = "test", fdrcutoff = 0.2,
                      correction = TRUE, seed = 44, ncore = 2,
                      returnObject = TRUE, verbose = FALSE)
  expect_equal(serial@misc$EpiLink$res_df, parallel@misc$EpiLink$res_df)
  expect_equal(serial@misc$EpiLink$threshold, parallel@misc$EpiLink$threshold)
  expect_equal(serial@misc$EpiLink$res_df$pair_seed, epilink_pair_seeds(2, 44))
  expect_true(serial@misc$EpiLink$correction)
})
