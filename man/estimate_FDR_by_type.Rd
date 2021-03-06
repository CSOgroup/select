% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_FDR.R
\name{estimate_FDR_by_type}
\alias{estimate_FDR_by_type}
\title{Interaction type aware FDR estimation}
\usage{
estimate_FDR_by_type(
  results,
  real,
  random,
  var = "int_type",
  FDR.thresh = 0.1,
  invert.scores = FALSE
)
}
\arguments{
\item{temp}{list of m matrices (n x n)}
}
\value{
n lists of n vectors of m elements
}
\description{
Interaction type aware FDR estimation
}
