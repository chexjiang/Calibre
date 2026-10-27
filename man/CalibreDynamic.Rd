% Generated from R/Calibre.R; keep in sync with the roxygen block.
\name{CalibreDynamic}
\alias{CalibreDynamic}
\alias{EpiLinkDynamic}
\title{Run Calibre-dynamic on a Seurat multiome object}
\usage{
CalibreDynamic(
  object,
  gene_peak_df,
  celltype_col = "newCT",
  family = "poisson",
  celltype = NULL,
  pseudotime_col = "pseudotime",
  pseudotime_formula = "s(pseudotime, k = 5, bs = 'tp')",
  covariates = NULL,
  null_generation = "s",
  stat = "zval",
  bin = TRUE,
  selection = "rate",
  returnObject = FALSE,
  fdrcutoff = 0.05,
  ncore = 10,
  theory = FALSE,
  conservative = FALSE,
  correction = FALSE,
  verbose = FALSE,
  seed = 123
)
}
\arguments{
\item{object}{Seurat object with RNA and ATAC assays.}
\item{gene_peak_df}{Candidate pairs with \code{gene} and \code{peak} columns.}
\item{celltype_col}{Cell-type metadata column.}
\item{family}{Either \code{"poisson"} or \code{"nb"}.}
\item{celltype}{Optional cell type; \code{NULL} uses all cells.}
\item{pseudotime_col}{Pseudotime metadata column.}
\item{pseudotime_formula}{An \code{s(pseudotime, ...)} smooth term, for example \code{"s(pseudotime, k = 5, bs = 'tp')"}.}
\item{covariates}{Optional cell-level model terms.}
\item{null_generation}{Null generator. Formal analyses use standardized Pearson residual resampling, \code{"s"}.}
\item{stat}{Selection statistic, either \code{"zval"} or \code{"pval"}.}
\item{bin}{Whether to binarize ATAC counts.}
\item{selection}{Candidate-pair filtering strategy.}
\item{returnObject}{Whether to return the Seurat object instead of significant rows.}
\item{fdrcutoff}{Target FDR.}
\item{ncore}{Number of forked workers.}
\item{theory}{Return observed GAM results without null generation.}
\item{conservative}{Use the conservative FDR numerator offset.}
\item{correction}{Apply the finite-cell correction in \code{CalibreFilter}.}
\item{verbose}{Report progress.}
\item{seed}{Master seed used to generate pair-specific seeds.}
}
\value{A Seurat object, significant-pair data frame, or theory result table.}
\description{
Fits a Poisson or negative-binomial GAM with an enhancer-by-pseudotime smooth
and calibrates the dynamic-link threshold using a standardized-residual nullstrap.
}
