% Generated from R/Calibre.R; keep in sync with the roxygen block.
\name{CalibreDynamic}
\alias{CalibreDynamic}
\alias{EpiLinkDynamic}
\title{Run Calibre-dynamic on a Seurat multiome object}
\usage{
CalibreDynamic(
  object,
  gene_peak_df,
  pseudotime_col = "pseudotime",
  pseudotime_formula = "quadratic",
  rescue_pseudotime_formula = "s(pseudotime, k = 5, bs = 'tp')",
  covariates = NULL,
  fdrcutoff = 0.1,
  ncore = 10,
  stat = "chisq",
  family = "poisson",
  null_generation = "s",
  bin = TRUE,
  selection = "rate",
  returnObject = FALSE,
  theory = FALSE,
  conservative = FALSE,
  correction = FALSE,
  verbose = FALSE,
  seed = 123,
  pair_seeds = NULL,
  rescue_dispersion = 1000,
  progress_file = NULL,
  progress_every = 10000L,
  checkpoint_dir = NULL,
  return_dynamic = FALSE,
  dynamic_grid = 100L
)
}
\arguments{
\item{object}{A Seurat object containing \code{RNA} and \code{ATAC} assays,
already subsetted to the pre-specified lineage. All cells in the input object
are used.}
\item{gene_peak_df}{A data frame of candidate gene--CRE pairs.}
\item{pseudotime_col}{Metadata column containing finite numeric pseudotime
values. The default is \code{"pseudotime"}.}
\item{pseudotime_formula}{Model for the dynamic effect. The default is
\code{"quadratic"}.}
\item{rescue_pseudotime_formula}{Formula used for automatic dispersion rescue.
The default is \code{"s(pseudotime, k = 5, bs = 'tp')"}.}
\item{covariates}{Optional cell-level covariates included in the GAM. The
default is \code{NULL}.}
\item{fdrcutoff}{Target FDR level. The default is \code{0.1}.}
\item{ncore}{Number of parallel workers on supported platforms. The default is
\code{10}.}
\item{stat}{Test statistic used to rank dynamic associations for FDR
thresholding; the default is \code{"chisq"}.}
\item{family}{GLM family. The default is \code{"poisson"}, with \code{"nb"} as
the negative binomial alternative.}
\item{null_generation}{Synthetic-null data generation scheme. The default,
\code{"s"}, uses standardized Pearson-residual resampling.}
\item{bin}{Whether to binarize ATAC counts before model fitting. The default is
\code{TRUE}.}
\item{selection}{Candidate-pair filtering strategy. The default is
\code{"rate"}.}
\item{returnObject}{Whether to return the updated Seurat object. The default is
\code{FALSE}, which returns the significant CRE--gene data frame.}
\item{theory}{Return observed GAM results without null generation. The default is \code{FALSE}.}
\item{conservative}{Use the conservative FDR numerator offset. The default is \code{FALSE}.}
\item{correction}{Apply the finite-cell correction in \code{CalibreFilter}. The default is \code{FALSE}.}
\item{verbose}{Report progress. The default is \code{FALSE}.}
\item{seed}{Master seed used to generate pair-specific seeds. The default is \code{123}.}
\item{pair_seeds}{Optional pair-specific seeds aligned to the input rows of \code{gene_peak_df}. The default is \code{NULL}.}
\item{rescue_dispersion}{For quadratic fits, refit pairs above this dispersion using \code{rescue_pseudotime_formula}. The default is \code{1000}; \code{NULL} disables automatic rescue.}
\item{progress_file}{Optional path to append parent-process progress logs. The default is \code{NULL}.}
\item{progress_every}{Number of pairs per progress batch. The default is \code{10000L}.}
\item{checkpoint_dir}{Optional directory for per-batch checkpoint files. The default is \code{NULL}.}
\item{return_dynamic}{Return compact dynamic parameters that can be evaluated at any pseudotime. No pair-by-cell table is created. The default is \code{FALSE}.}
\item{dynamic_grid}{Deprecated compatibility argument; ignored. The default is \code{100L}.}
}
\value{A Seurat object, significant-pair data frame, theory result table, or—when
\code{return_dynamic = TRUE}—a result containing compact dynamic
parameters.}
\description{
Fits a Poisson or negative-binomial GAM with an enhancer-by-pseudotime smooth
and calibrates the dynamic-link threshold using a standardized-residual nullstrap.
The input object is assumed to be pre-subsetted to one lineage, and all input
cells are used.
}
