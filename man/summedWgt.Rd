% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R-PRES.r
\name{summedWgt}
\alias{summedWgt}
\title{Calculate summed model weights for a set of variables for a specific parameter.}
\usage{
summedWgt(covnames, param, aic.tab)
}
\arguments{
\item{covnames}{a vector of string values containing the names of the variables of interest.}

\item{param}{a string with the name of the parameter of interest, e.g., \code{param="psi"}.}

\item{aic.tab}{an \code{aic.tab} object that includes the set of models from which the
summed model weights are to be calculated.}
}
\value{
returns a data frame with the summed weights and evidence ratios.
}
\description{
Calculates the summed model weights (and evidence ratios) for a set of variables
  \code{covnames}, for a parameter \code{param} and a set of models stored in an
  \code{aic.tab} object.
}
\seealso{
\code{\link{createAicTable}}
}
\author{
Darryl MacKenzie.
}
\keyword{htest}
