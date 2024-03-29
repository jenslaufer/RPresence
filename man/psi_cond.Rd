% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R-PRES.r
\name{psi_cond}
\alias{psi_cond}
\title{get conditional occupancy estimates from single-season ("so") model}
\usage{
psi_cond(m0, paodata)
}
\arguments{
\item{m0}{model object (result from running "occMod").}

\item{paodata}{pao data object}
}
\value{
numeric matrix with column of estimates and column of std. errors
}
\description{
\code{psicond} Computes conditional (on detection history) occupancy estimates
  from single-season ("so") model.  Uses "delta" method, with numerical gradient.
}
