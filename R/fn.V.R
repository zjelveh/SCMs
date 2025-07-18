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
  # kernlab functions imported via NAMESPACE
  
  # Validate variables.v before creating diagonal matrix
  if (is.null(variables.v) || length(variables.v) <= 0) {
    stop("variables.v is NULL or has length 0")
  }
  
  if (any(is.na(variables.v))) {
    stop("variables.v contains NA values")
  }
  
  if (any(!is.finite(variables.v))) {
    stop("variables.v contains non-finite values")
  }
  
  var_length <- length(variables.v)
  if (var_length <= 0 || is.na(var_length) || !is.finite(var_length)) {
    stop(paste("Invalid length of variables.v:", var_length))
  }
  
  # Create diagonal V matrix from variables
  V <- diag(x = as.numeric(abs(variables.v) / sum(abs(variables.v))),
            nrow = var_length,
            ncol = var_length)
  
  # Validate input matrices
  if (is.null(X0.scaled) || !is.matrix(X0.scaled)) {
    stop("X0.scaled is NULL or not a matrix")
  }
  
  if (is.null(X1.scaled) || !is.matrix(X1.scaled)) {
    stop("X1.scaled is NULL or not a matrix")
  }
  
  if (nrow(X0.scaled) != var_length) {
    stop(paste("Dimension mismatch: X0.scaled has", nrow(X0.scaled), "rows but variables.v has length", var_length))
  }
  
  if (nrow(X1.scaled) != var_length) {
    stop(paste("Dimension mismatch: X1.scaled has", nrow(X1.scaled), "rows but variables.v has length", var_length))
  }
  
  # Set up quadratic programming problem
  H <- t(X0.scaled) %*% V %*% (X0.scaled)
  a <- X1.scaled
  c <- -1 * c(t(a) %*% V %*% (X0.scaled))
  A <- t(rep(1, length(c)))
  b <- 1
  l <- rep(0, length(c))
  u <- rep(1, length(c))
  r <- 0
  
  # Check if clarabel should be used for V optimization
  use_clarabel <- getOption("SCMs.prefer_clarabel", FALSE)
  
  if (use_clarabel && requireNamespace("clarabel", quietly = TRUE)) {
    # Use clarabel for QP solving
    tryCatch({
      solution.w <- solve_qp_clarabel(c, H, A, b, l, u)
      attr(solution.w, "v_solver_used") <- "clarabel"
    }, error = function(e) {
      # If clarabel fails, fall back to kernlab::ipop but log the fallback
      warning(paste0("Clarabel failed for V optimization, falling back to kernlab::ipop: ", e$message))
      res <- kernlab::ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r, 
                           bound = bound.ipop,
                           margin = margin.ipop, 
                           maxiter = 1000, 
                           sigf = sigf.ipop)
      
      # Extract solution weights
      solution.w <- as.matrix(kernlab::primal(res))
      attr(solution.w, "v_solver_used") <- "kernlab::ipop"
      attr(solution.w, "v_fallback_reason") <- e$message
    })
  } else {
    # Use kernlab::ipop for QP solving
    res <- kernlab::ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r, 
                         bound = bound.ipop,
                         margin = margin.ipop, 
                         maxiter = 1000, 
                         sigf = sigf.ipop)
    
    # Extract solution weights
    solution.w <- as.matrix(kernlab::primal(res))
    attr(solution.w, "v_solver_used") <- "kernlab::ipop"
  }
  
  # Compute loss for features
  loss.w <- as.numeric(t(X1.scaled - X0.scaled %*% solution.w) %*%
                         (V) %*% (X1.scaled - X0.scaled %*% solution.w))
  
  # Compute loss for outcomes
  loss.v <- as.numeric(t(Z1 - Z0 %*% solution.w) %*%
                         (Z1 - Z0 %*% solution.w))
  loss.v <- loss.v / nrow(Z0)
  
  return(invisible(loss.v))
}

#' @title Solve QP using Clarabel for V optimization
#' @description Solve quadratic programming problem using clarabel solver
#' @param c Linear objective coefficients
#' @param H Quadratic objective matrix (Hessian)
#' @param A Equality constraint matrix
#' @param b Equality constraint RHS
#' @param l Lower bounds
#' @param u Upper bounds
#' @return Solution vector as matrix
solve_qp_clarabel <- function(c, H, A, b, l, u) {
  # Convert kernlab::ipop problem to clarabel format
  # ipop: minimize c'x + 0.5*x'Hx subject to Ax = b, l <= x <= u
  # clarabel: minimize 0.5*x'Px + q'x subject to Ax + s = b, s in cones
  
  n_vars <- length(c)
  
  # Clarabel QP formulation
  P <- H  # Quadratic term (ipop uses 0.5*x'Hx, clarabel uses 0.5*x'Px)
  q <- c  # Linear term
  
  # Constraints: Ax = b, l <= x <= u
  # Reformulate as: 
  # 1. Ax = b (zero cone)
  # 2. -x + s_l = -l, s_l >= 0 (lower bounds)
  # 3. x + s_u = u, s_u >= 0 (upper bounds)
  
  # Build constraint matrix
  A_eq <- A  # Equality constraints Ax = b
  A_lower <- -diag(n_vars)  # Lower bound constraints -x + s_l = -l
  A_upper <- diag(n_vars)   # Upper bound constraints x + s_u = u
  
  A_constraints <- rbind(A_eq, A_lower, A_upper)
  b_constraints <- c(b, -l, u)
  
  # Cone specification
  n_eq <- nrow(A_eq)
  n_lower <- length(l)
  n_upper <- length(u)
  
  # Call clarabel solver
  control_settings <- clarabel::clarabel_control(
    verbose = FALSE,
    max_iter = 1000L,
    tol_gap_abs = 1e-8,
    tol_gap_rel = 1e-8
  )
  
  result <- clarabel::clarabel(
    A = A_constraints,
    b = b_constraints,
    q = q,
    P = P,
    cones = list(z = n_eq, l = n_lower + n_upper),  # Zero cones + linear cones
    control = control_settings
  )
  
  # Check clarabel status
  if (result$status != 2) {
    status_desc <- names(clarabel::solver_status_descriptions())[result$status]
    stop(paste0("Clarabel V optimization failed! Status: ", result$status, 
                " (", status_desc, "). Problem may be infeasible or ill-conditioned."))
  }
  
  # Extract solution (first n_vars components)
  solution <- result$x[1:n_vars]
  
  # Return as matrix to match ipop format
  return(as.matrix(solution))
}