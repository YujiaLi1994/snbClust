% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sim.Independent.R
\name{Sim.Independent}
\alias{Sim.Independent}
\title{Simulation 1}
\usage{
Sim.Independent(
  ngenes = 150,
  eff_a = 1,
  percent_DE = 1,
  sim_disp,
  empirical_dist
)
}
\arguments{
\item{ngenes}{total number of genes}

\item{eff_a}{effect size}

\item{percent_DE}{percentage of DE genes}

\item{sim_disp}{gene dispersion from breast cancer data}

\item{empirical_dist}{average gene expression from breast cancer data}
}
\value{
data with independent genes for simulation 1 and 2
}
\description{
Simulation 1
}
\references{
Rahman T, Li Y, Ma T, et al. A sparse negative binomial mixture model for clustering RNA-seq count data[J]. arXiv preprint arXiv:1912.02399, 2019.
}
