% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_random.R
\name{gen.random.am}
\alias{gen.random.am}
\title{Generate random alteration matrix. Enables parallel calculation}
\usage{
gen.random.am(
  M,
  N = 1000,
  permuteFun = birewireBlocEvents,
  n.cores = 1,
  folder = "./",
  r.seed = 100,
  ...
)
}
\arguments{
\item{am}{Input alteration matrix}

\item{params}{parameter list}
}
\value{
A collection of random Alteration matrices
}
\description{
Generate random alteration matrix. Enables parallel calculation
}
