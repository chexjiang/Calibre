test_that("candidate filtering uses cell-specific matrices", {
  rna <- Matrix::Matrix(matrix(c(1, 0, 1, 0, 0, 0), nrow = 2, byrow = TRUE), sparse = TRUE)
  atac <- Matrix::Matrix(matrix(c(1, 1, 0, 0, 0, 0), nrow = 2, byrow = TRUE), sparse = TRUE)
  rownames(rna) <- c("g1", "g2")
  rownames(atac) <- c("chr1-100-200", "chr1-300-400")
  pairs <- data.frame(gene = c("g1", "g2"),
                      peak = c("chr1-100-200", "chr1-300-400"))
  colnames(rna) <- colnames(atac) <- paste0("cell", 1:3)
  object <- Seurat::CreateSeuratObject(rna, assay = "RNA")
  object[["ATAC"]] <- Signac::CreateChromatinAssay(atac)
  kept <- FilterGenePeakPairs(object, pairs, expression.rate = 0.1, peak.rate = 0.1)
  expect_equal(kept$gene, "g1")
  expect_equal(kept$peak, "chr1-100-200")
})
