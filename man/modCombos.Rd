% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R-PRES.r
\name{modCombos}
\alias{modCombos}
\title{Create a list of formulae for all possible first-order combinations for a set of variables.}
\usage{
modCombos(param, covs)
}
\arguments{
\item{param}{a string with the name of the parameter of interest, e.g., \code{param="psi"}.}

\item{covs}{a vector of string values containing the names of the variables of interest.}
}
\value{
returns a list of formulae of each combination of variables.
}
\description{
Creates a list of formulae for the parameter \code{param} which contains all possible
first-order (i.e., no interactions) combinations for the set of variables \code{cov}.
}
\seealso{
\code{\link{createAicTable}}, \code{\link{predict.occMod}}
}
\author{
Darryl MacKenzie.
}
\keyword{models}
