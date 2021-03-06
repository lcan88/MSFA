% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bayesian.R
\name{sp_eigen}
\alias{sp_eigen}
\title{Normalized eigenvalues of a symmetric matrix.}
\usage{
sp_eigen(SigPhi)
}
\arguments{
\item{SigPhi}{Symmetric matrix, of size \eqn{p \times p}{p x p}.}
}
\value{
The normalized eigenvalues of \code{SigPhi}, a vector of length \eqn{p}{p}.
}
\description{
This is a rather simple function, performing the eigenvalue computation for a
covariance matrix. Used for the selection of common latent factor dimensions.
}
