% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snbClust.fusedlasso.R
\name{em_glm_Yujia}
\alias{em_glm_Yujia}
\title{sparse negative binomial mixture model for clustering RNA-seq count data with fused lasso penalty.}
\usage{
em_glm_Yujia(
  data,
  phi = 1,
  lib,
  c_center = NULL,
  k = 3,
  lambda,
  penalty.option = "MCP",
  tune = 1,
  ita = 3,
  max_iter_IRLS = 20,
  max_iter_ADMM = 50,
  max_iter_EM = 10,
  threshold = 0.005,
  BIC.ita = 1
)
}
\arguments{
\item{data}{data matrix, where rows represent genes and columns represents samples.}

\item{phi}{dispersion parameter for each gene.}

\item{lib}{library size normalization factor for each sample.}

\item{k}{number of clusters.}

\item{lambda}{penalty parameter needs to be tuned or specified.}

\item{penalty.option}{penalty type, the default is MCP.}

\item{tune}{a parameter controlling the ADMM algorithm}

\item{ita}{a parameter in MCP penalty}

\item{center}{Initial centers for each clusters. We suggest using
other clustering algorithms such as sparse K-means to speed up.}
}
\value{
A list of two components:
\itemize{
\item{result: }{a list containing mean matrix (beta),
likelihood, penalized log likelihood and cluster assignment probability (z)}
\item{BIC: }{ Extended BIC of the result}

}
}
\description{
sparse negative binomial mixture model for clustering RNA-seq count data with fused lasso penalty
}
\references{
Rahman T, Li Y, Ma T, et al. A sparse negative binomial mixture model for clustering RNA-seq count data[J]. arXiv preprint arXiv:1912.02399, 2019.
}
