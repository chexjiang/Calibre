% Generated from R/Calibre.R; keep in sync with the roxygen block.
\name{CalibreDynamicPredict}
\alias{CalibreDynamicPredict}
\title{Evaluate saved Calibre-dynamic spline parameters}
\usage{
CalibreDynamicPredict(dynamic_parameters, pseudotime)
}
\arguments{
\item{dynamic_parameters}{One element of \code{dynamic_parameters} returned by \code{CalibreDynamic()}.}
\item{pseudotime}{Finite numeric pseudotime values at which to evaluate the model.}
}
\value{A data frame containing the dynamic component, total effect, optional
standard error and confidence interval, and fold change.}
\description{Reconstructs the dynamic accessibility effect from the saved
mgcv spline basis and coefficients without refitting the model.}
