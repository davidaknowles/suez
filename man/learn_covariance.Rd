% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/learn_covariance.R
\name{learn_covariance}
\alias{learn_covariance}
\title{Function to learn covariance structure across samples using latent factors}
\usage{
learn_covariance(
  input,
  anno,
  kinship_matrix = NULL,
  P = 20,
  same_ind_initial_weight = 0.01,
  kinship_matrix_initial_weight = 0.01,
  ...
)
}
\arguments{
\item{input}{[genes x samples] Normalized (e.g. quantile normalized) matrix of log2 gene expression}

\item{anno}{[samples x 2] data.frame with `individual` and `condition` columns. If a kinship matrix is provided the row and column names must correspond to the values used in the `individual` column.}

\item{kinship_matrix}{Representing relatedness between individuals. Optional.}

\item{P}{Number of latent factors to model. suez can prune out unnecessary factors (in principle).}

\item{same_ind_initial_weight}{Initialization.}

\item{kinship_matrix_initial_weight}{Initialization.}

\item{...}{passed to rstan::optimizing}
}
\description{
Function to learn covariance structure across samples using latent factors
}
