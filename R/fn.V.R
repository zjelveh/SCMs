#' @title Compute Loss for Synthetic Control Weights
#' @description This function computes the loss for synthetic control weights given feature matrices and outcome vectors.
#'
#' @param variables.v Numeric vector. Variables for the diagonal of the V matrix.
#' @param X0.scaled Matrix. Scaled feature matrix for control units.
#' @param X1.scaled Matrix. Scaled feature matrix for treated unit.
#' @param Z0 Matrix. Outcome matrix for control units.
#' @param Z1 Matrix. Outcome matrix for treated unit.
#' @param margin.ipop Numeric. Margin parameter for interior point optimization. Default is 0.0005.
#' @param sigf.ipop Numeric. Significant digits for ipop. Default is 5.
#' @param bound.ipop Numeric. Bound for ipop. Default is 10.
#'
#' @return Numeric. The computed loss value.
#'
#' @examples
#' # Example usage (replace with actual example when available)
#' # loss <- fn.V(variables.v = c(0.5, 0.5), X0.scaled = matrix(...), X1.scaled = matrix(...),
#' #              Z0 = matrix(...), Z1 = matrix(...))
fn.V <- function(
    variables.v = stop("variables.v missing"),
    X0.scaled = stop("X0.scaled missing"),
    X1.scaled = stop("X1.scaled missing"),
    Z0 = stop("Z0 missing"),
    Z1 = stop("Z1 missing"),
    margin.ipop = 0.0005,
    sigf.ipop = 5,
    bound.ipop = 10
) {
  library(kernlab)
  
  # Create diagonal V matrix from variables
  V <- diag(x = as.numeric(abs(variables.v) / sum(abs(variables.v))),
            nrow = length(variables.v),
            ncol = length(variables.v))
  
  # Set up quadratic programming problem
  H <- t(X0.scaled) %*% V %*% (X0.scaled)
  a <- X1.scaled
  c <- -1 * c(t(a) %*% V %*% (X0.scaled))
  A <- t(rep(1, length(c)))
  b <- 1
  l <- rep(0, length(c))
  u <- rep(1, length(c))
  r <- 0
  
  # Solve quadratic programming problem using interior point optimization
  res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r, 
              bound = bound.ipop,
              margin = margin.ipop, 
              maxiter = 1000, 
              sigf = sigf.ipop)
  
  # Extract solution weights
  solution.w <- as.matrix(primal(res))
  
  # Compute loss for features
  loss.w <- as.numeric(t(X1.scaled - X0.scaled %*% solution.w) %*%
                         (V) %*% (X1.scaled - X0.scaled %*% solution.w))
  
  # Compute loss for outcomes
  loss.v <- as.numeric(t(Z1 - Z0 %*% solution.w) %*%
                         (Z1 - Z0 %*% solution.w))
  loss.v <- loss.v / nrow(Z0)
  
  return(invisible(loss.v))
}