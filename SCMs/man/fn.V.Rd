% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn.V.R
\name{fn.V}
\alias{fn.V}
\title{Compute Loss for Synthetic Control Weights}
\usage{
fn.V(
  variables.v = stop("variables.v missing"),
  X0.scaled = stop("X0.scaled missing"),
  X1.scaled = stop("X1.scaled missing"),
  Z0 = stop("Z0 missing"),
  Z1 = stop("Z1 missing"),
  margin.ipop = 5e-04,
  sigf.ipop = 5,
  bound.ipop = 10
)
}
\arguments{
\item{variables.v}{Numeric vector. Variables for the diagonal of the V matrix.}

\item{X0.scaled}{Matrix. Scaled feature matrix for control units.}

\item{X1.scaled}{Matrix. Scaled feature matrix for treated unit.}

\item{Z0}{Matrix. Outcome matrix for control units.}

\item{Z1}{Matrix. Outcome matrix for treated unit.}

\item{margin.ipop}{Numeric. Margin parameter for interior point optimization. Default is 0.0005.}

\item{sigf.ipop}{Numeric. Significant digits for ipop. Default is 5.}

\item{bound.ipop}{Numeric. Bound for ipop. Default is 10.}
}
\value{
Numeric. The computed loss value.
}
\description{
This function computes the loss for synthetic control weights given feature matrices and outcome vectors.
}
\examples{
# Example usage (replace with actual example when available)
# loss <- fn.V(variables.v = c(0.5, 0.5), X0.scaled = matrix(...), X1.scaled = matrix(...),
#              Z0 = matrix(...), Z1 = matrix(...))
}
