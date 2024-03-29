% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R-PRES.r
\name{summary.aicTable}
\alias{summary.aicTable}
\title{Summary of  \code{aic.table} object}
\usage{
\method{summary}{aicTable}(object, ...)
}
\arguments{
\item{object}{an object of class \code{aicTable} which is the result of the function
\code{\link{createAicTable}}}

\item{...}{possible additional arguments}
}
\value{
Produces a summary of the model comparison using AIC as a data frame that
includes the model name, relative difference in AIC, model weight, number of parameters,
twice the negative log-likelihood and any warnings.
}
\description{
Provides a brief summary of an object of class \code{aicTable}, which is the result of
the function \code{\link{createAicTable}}.
}
\author{
Darryl MacKenzie
}
\keyword{htest}
