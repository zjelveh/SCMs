# Auxiliary function that creates the constraints to be passed to the optimization problem
w.constr.OBJ <- function(w.constr, A, Z, V, J, KM, M) {
  # Default method to estimate weights as in Abadie et al. (2010)
  if (is.null(w.constr)) {
    w.constr <- list(lb = 0,
                     p = "L1",
                     dir = "==",
                     Q = 1,
                     name = "simplex")
    
  } else if (w.constr[["name"]] == "simplex") {
    
    if (!"Q" %in% names(w.constr)) {
      Q <- 1
    } else {
      Q <- w.constr[["Q"]]
    }
    
    w.constr <- list(lb = 0,
                     p = "L1",
                     dir = "==",
                     Q = Q,
                     name = 'simplex')
    
  } else if (w.constr[["name"]] == "ols") {
    w.constr <- list(lb = -Inf,
                     dir = "NULL",
                     p = "no norm",
                     name = 'ols')
    
  } else if (w.constr[["name"]] == "lasso") {
    
    if (!"Q" %in% names(w.constr)) {
      w.constr[["Q"]] <- shrinkage.EST("lasso", A, as.matrix(Z), V, J, KM)$Q
    }
    
    w.constr <- list(lb = -Inf,
                     p = "L1",
                     dir = "<=",
                     Q = w.constr[["Q"]],
                     name = 'lasso')
    
  } else if (w.constr[["name"]] == "ridge") {
    
    if (!("Q" %in% names(w.constr))) {
      
      feature.id <- unlist(purrr::map(stringr::str_split(rownames(Z), "\\."), 2))
      
      Qfeat <- c()
      for (feat in unique(feature.id)) {
        Af <- A[feature.id == feat,, drop = FALSE]
        Zf <- Z[feature.id == feat,, drop = FALSE]
        Vf <- V[feature.id == feat, feature.id == feat, drop = FALSE]
        if (nrow(Af) >= 2) {
          #5) {
          QQ <- tryCatch({
            aux <- shrinkage.EST("ridge", Af, as.matrix(Zf), Vf, J, KM)
            Q <- aux$Q
          }, warning = function(warn) { },
          error = function(err) { }, finally = { })
          Qfeat <- c(Qfeat, QQ)
        }
      }
      
      if (is.null(Qfeat)) {
        aux <- shrinkage.EST("ridge", A, as.matrix(Z), V, J, KM)
        Qfeat <- aux$Q
      }
      # Filter out infinite values and use a reasonable default if needed
      Qfeat_finite <- Qfeat[is.finite(Qfeat)]
      if (length(Qfeat_finite) > 0) {
        w.constr[["Q"]] <- max(min(Qfeat_finite), .5)
      } else {
        w.constr[["Q"]] <- 1.0  # Default Q value
      }
      if (exists("aux") && is.finite(aux$lambda)) {
        w.constr[["lambda"]] <- aux$lambda
      } else {
        w.constr[["lambda"]] <- 0.1  # Default lambda value
      }
      
    }
    
    w.constr <- list(lb = -Inf,
                     p = "L2",
                     dir = "<=",
                     Q = w.constr[["Q"]],
                     name = "ridge",
                     lambda = w.constr[["lambda"]])
    
  } else if (w.constr[["name"]] == "L1-L2") {
    
    if (!("Q2" %in% names(w.constr))) {
      
      feature.id <- unlist(purrr::map(stringr::str_split(rownames(Z), "\\."), 2))
      Qfeat <- c()
      for (feat in unique(feature.id)) {
        Af <- A[feature.id == feat,, drop = FALSE]
        Zf <- Z[feature.id == feat,, drop = FALSE]
        Vf <- V[feature.id == feat, feature.id == feat, drop = FALSE]
        
        if (nrow(Af) >= 5) {
          QQ <- tryCatch({
            aux <- shrinkage.EST("ridge", Af, as.matrix(Zf), Vf, J, KM)
            Q2 <- aux$Q
          }, warning = function(warn) { },
          error = function(err) { }, finally = { })
          Qfeat <- c(Qfeat, QQ)
        }
      }
      
      if (is.null(Qfeat)) {
        aux <- shrinkage.EST("ridge", A, as.matrix(Z), V, J, KM)
        Qfeat <- aux$Q
      }
      # Filter out infinite values and use a reasonable default if needed
      Qfeat_finite <- Qfeat[is.finite(Qfeat)]
      if (length(Qfeat_finite) > 0) {
        w.constr[["Q2"]] <- max(min(Qfeat_finite), .5)
      } else {
        w.constr[["Q2"]] <- 1.0  # Default Q2 value
      }
      if (exists("aux") && is.finite(aux$lambda)) {
        w.constr[["lambda"]] <- aux$lambda
      } else {
        w.constr[["lambda"]] <- 0.1  # Default lambda value
      }
    }
    
    w.constr <- list(lb = 0,
                     p = "L1-L2",
                     dir = "==/<=",
                     Q = 1,
                     Q2 = w.constr[["Q2"]],
                     name = "L1-L2",
                     lambda = w.constr[["lambda"]])
    
  } else if (w.constr[["name"]] == "pensynth") {
    
    # Validate lambda (CV should have been done at higher level)
    if (!is.numeric(w.constr[["lambda"]]) || w.constr[["lambda"]] <= 0) {
      stop("Pensynth lambda must be a positive number")
    }
    
    w.constr <- list(lb = 0,
                     p = "pensynth", 
                     dir = "simplex",
                     Q = 1,
                     lambda = w.constr[["lambda"]],
                     name = "pensynth")
    
  } else {
    
    # if constraint is entirely user specified just check everything is fine
    if (!(all(c('p', 'dir', 'Q', 'lb') %in% names(w.constr)))) {
      stop("If 'name' is not specified, w.constr should be a list whose elements
            must be named 'p','dir','Q','lb'.")
    }
    
    if (!(w.constr[["p"]] %in% c("no norm", "L1", "L2", "L1-L2"))) {
      stop("Specify either p = 'no norm' (no constraint on the norm of weights),
          p = 'L1' (L1-norm), p = 'L2' (L2-norm)")
    } else if (w.constr[["p"]] == "no norm") {
      w.constr[["dir"]] <- "NULL"
    }
    
    if (!(w.constr[["dir"]] %in% c("<=", "==", "==/<=", "NULL"))) {
      stop("Specify either dir = '<=' (inequality constraint on the norm of the weights)
            or dir = '==' (equality constraint on the norm of the weights) or dir = 'NULL'
            in case you don't want to specify a constraint on the norm of the weights.")
    }
    
    if (!(w.constr[["lb"]] == 0 || w.constr[["lb"]] == -Inf)) {
      stop("Specify either lb = 0 or lb = -Inf.")
    }
    
    w.constr[["lb"]] <- w.constr[["lb"]]
    w.constr[["name"]] <- "user provided"
  }
  return(w.constr)
}

# Estimate optimal lambda for pensynth using K-fold cross-validation
estimate_optimal_lambda_cv <- function(A, Z, V, Y.donors, Y.pre, nlambda = 20) {
  # A is the treated unit aggregated features (vector)
  # Z is the donor aggregated features matrix (features x donors)  
  # V is the feature weighting matrix (features x features)
  # Y.donors is the donor outcome time series (T_pre x J)
  # Y.pre is the treated unit outcome time series (T_pre x 1)
  
  # Ensure inputs are matrices
  if (!is.matrix(Y.donors)) Y.donors <- as.matrix(Y.donors)
  if (!is.matrix(Y.pre)) Y.pre <- as.matrix(Y.pre)
  
  J <- ncol(Z)  # Number of donors
  T_pre <- nrow(Y.donors)  # Number of pre-treatment periods
  
  if (T_pre < 2) {
    stop("Need at least 2 pre-treatment periods for hold-out validation")
  }
  
  # Pre-compute V-weighted matrices for feature matching
  v_weights <- diag(V)
  sqrt_v <- sqrt(v_weights)
  A_weighted <- sqrt_v * A
  Z_weighted <- sqrt_v * Z
  
  # Calculate components for QP (following original pensynth)
  X0VX0 <- crossprod(Z_weighted)  # X0v'X0v (J x J)
  X1VX0 <- crossprod(A_weighted, Z_weighted)  # X1v'X0v (1 x J)
  Delta <- colSums(sweep(Z_weighted, 1, A_weighted, "-")^2)  # ||X1v - X0v_j||^2 for each donor j
  
  # Create lambda sequence following original pensynth approach
  lambda_seq <- lambda_sequence_pensynth(X1VX0, Delta, nlambda)
  
  # Hold-out validation: use outcome time series Z1, Z0
  # Split time series: use first T_pre-1 periods for training, last period for validation
  train_periods <- 1:(T_pre - 1)
  val_period <- T_pre
  
  Z1_train <- Y.pre[train_periods, , drop = FALSE]  # Training: treated outcomes
  Z0_train <- Y.donors[train_periods, , drop = FALSE]  # Training: donor outcomes
  Z1_val <- Y.pre[val_period, , drop = FALSE]  # Validation: treated outcome
  Z0_val <- Y.donors[val_period, , drop = FALSE]  # Validation: donor outcomes
  
  cv_errors <- numeric(nlambda)
  
  # For each lambda, fit pensynth model and evaluate on hold-out period
  for (i in 1:nlambda) {
    tryCatch({
      weights <- fit_pensynth_single_lambda(A, Z, V, lambda_seq[i], A_weighted, sqrt_v)
      
      # Predict held-out treated outcome using fitted weights
      pred_val <- Z0_val %*% weights
      
      # Calculate MSE on held-out period
      cv_errors[i] <- as.numeric((Z1_val - pred_val)^2)
      
    }, error = function(e) {
      cv_errors[i] <<- Inf
    })
  }
  
  # Return lambda with minimum CV error
  optimal_idx <- which.min(cv_errors)
  
  if (cv_errors[optimal_idx] == Inf) {
    stop("Cross-validation failed for all lambda values")
  }
  
  return(lambda_seq[optimal_idx])
}

# Lambda sequence calculation following original pensynth
lambda_sequence_pensynth <- function(X1VX0, Delta, nlambda) {
  lmin <- 1e-11
  lmax <- sum(abs(as.numeric(X1VX0)/Delta))
  return(exp(seq(log(lmin), log(lmax), length.out = nlambda)))
}

# Fit pensynth for a single lambda value
fit_pensynth_single_lambda <- function(A, Z, V, lambda, A_weighted = NULL, sqrt_v = NULL, use_clarabel = getOption("SCMs.prefer_clarabel", FALSE)) {
  J <- ncol(Z)
  
  # Use clarabel if available and requested
  if (use_clarabel && requireNamespace("clarabel", quietly = TRUE)) {
    return(fit_pensynth_single_lambda_clarabel(A, Z, V, lambda, A_weighted, sqrt_v))
  }
  
  # Use pre-computed weighted matrices if provided, otherwise compute them
  if (is.null(A_weighted) || is.null(sqrt_v)) {
    v_weights <- diag(V)
    sqrt_v <- sqrt(v_weights)
    A_weighted <- sqrt_v * A
  }
  
  # Calculate Z_weighted for current donor subset
  Z_weighted <- sqrt_v * Z  # Broadcasting: (K x 1) * (K x J) = (K x J)
  
  # Vectorized computation of ||X1v - X0v_j||^2 for all donors j
  # A_weighted is length K, Z_weighted is (K x J)
  # Use sweep() to subtract A_weighted from each column of Z_weighted, then compute column-wise squared norms
  Delta <- colSums(sweep(Z_weighted, 1, A_weighted, "-")^2)
  # Set up CVXR optimization
  x <- CVXR::Variable(J)
  
  # Combined objective: standard SCM fit + pensynth penalty
  base_objective <- CVXR::quad_form(A - Z %*% x, V)
  penalty_term <- lambda * sum(Delta * x)
  objective <- CVXR::Minimize(base_objective + penalty_term)
  
  # Standard simplex constraints
  constraints <- list(CVXR::sum_entries(x) == 1, x >= 0)
  # Solve optimization problem
  prob <- CVXR::Problem(objective, constraints)
  sol <- CVXR::solve(prob, solver = "ECOS", num_iter = 100000L, verbose = FALSE)

  if (is.null(sol$getValue) || !(sol$status %in% c("optimal", "optimal_inaccurate"))) {
    stop("Pensynth optimization failed")
  }
  
  weights <- sol$getValue(x)
  return(as.numeric(weights))
}

# Fast clarabel implementation for single lambda pensynth
fit_pensynth_single_lambda_clarabel <- function(A, Z, V, lambda, A_weighted = NULL, sqrt_v = NULL) {
  J <- ncol(Z)
  
  # Use pre-computed weighted matrices if provided, otherwise compute them
  if (is.null(A_weighted) || is.null(sqrt_v)) {
    v_weights <- diag(V)
    sqrt_v <- sqrt(v_weights)
    A_weighted <- sqrt_v * A
  }
  
  # Calculate Z_weighted for current donor subset
  Z_weighted <- sqrt_v * Z  # Broadcasting: (K x 1) * (K x J) = (K x J)
  
  # Vectorized computation of ||X1v - X0v_j||^2 for all donors j
  Delta <- colSums(sweep(Z_weighted, 1, A_weighted, "-")^2)
  
  # Set up QP formulation following original pensynth implementation
  # P = X0VX0 = crossprod(X0v), q = -X1VX0 + lambda*Delta
  P <- crossprod(Z_weighted)  # X0v'X0v (J x J)
  q <- -crossprod(A_weighted, Z_weighted) + lambda * Delta  # -X1v'X0v + lambda*Delta (length J)
  
  # Add small regularization to ensure P is positive definite
  P <- P + 1e-8 * diag(J)
  
  # Set up constraints: sum(x) = 1, x >= 0
  # For clarabel: Ax + s = b, where s ∈ K
  # We need: [1 1 ... 1] * x + s1 = 1 (equality, s1 = 0)
  #          [-I] * x + s2 = 0 (inequality, s2 >= 0, so -x + s2 = 0 => x = s2 >= 0)
  
  A_constraints <- rbind(
    matrix(1, nrow = 1, ncol = J),    # Sum constraint: sum(x) = 1
    -diag(J)                          # Non-negativity: -x + s = 0, s >= 0
  )
  b_constraints <- c(1, rep(0, J))    # RHS: [1, 0, 0, ..., 0]
  
  # Call clarabel solver with optimized settings for CV speed
  control_settings <- clarabel::clarabel_control(
    verbose = FALSE,
    max_iter = 5000L,  # Reduced for CV speed
    tol_gap_abs = 1e-6,  # Slightly relaxed for speed
    tol_gap_rel = 1e-6
  )
  
  result <- clarabel::clarabel(
    A = A_constraints,
    b = b_constraints,
    q = as.vector(q),
    P = P,
    cones = list(z = 1L, l = J),  # 1 zero cone (equality) + J linear cones (non-negativity)
    control = control_settings
  )

  # Check clarabel status - accept only solved status (numeric code 2 = "Solved")
  if (result$status != 2) {
    status_desc <- names(clarabel::solver_status_descriptions())[result$status]
    if (result$status == 3) {  # PrimalInfeasible
      stop(paste0("Pensynth optimization problem is primal infeasible! ", 
                  "The constraints cannot be satisfied. This may indicate: ",
                  "(1) lambda value is too large, (2) data contains invalid values, ",
                  "or (3) problem setup is incorrect. Status: ", result$status, " (", status_desc, ")"))
    } else {
      stop(paste0("Clarabel pensynth optimization failed! Status: ", result$status, 
                  " (", status_desc, "). Check problem setup, lambda value, or disable clarabel with options(SCMs.prefer_clarabel = FALSE)"))
    }
  }
  
  return(as.numeric(result$x))
}

shrinkage.EST <- function(method, A, Z, V, J, KM) {
  
  lambd <- NULL
  if (method == "lasso") Q <- 1
  
  if (method == "ridge") {
    wls <- lm(A ~ Z - 1, weights = diag(V))
    sig.wls <- sigma(wls)
    lambd <- sig.wls ^ 2 * (J + KM) / sum(wls$coef ^ 2, na.rm = TRUE) # rule of thumb for lambda (Hoerl et al, 1975)
    Q <- sqrt(sum(wls$coef ^ 2, na.rm = TRUE)) / (1 + lambd) # convert lambda into Q
    
    # reduce dimensionality of the problem if more params than obs
    if (is.nan(Q) || (nrow(Z) <= ncol(Z) + 10)) {
      lasso.cols <- b.est(A = A, Z = Z, J = J, KM = KM, V = V, CVXR.solver = "OSQP",
                          w.constr = list(name = "lasso", dir = "<=", lb = 0, p = "L1", Q = 1))
      active.cols <- abs(lasso.cols) > 1e-8
      if (sum(active.cols) >= (max(nrow(A) - 10, 2))) {
        active.cols <- rank(-abs(lasso.cols)) <= max(nrow(A) - 10, 2)
      }
      Z.sel <- Z[, active.cols, drop = FALSE]
      wls <- lm(A ~ Z.sel - 1, weights = diag(V))
      sig.wls <- sigma(wls)
      lambd <- sig.wls ^ 2 * (ncol(Z.sel) + KM) / sum(wls$coef ^ 2, na.rm = TRUE) # rule of thumb for lambda (Hoerl et al, 1975)
      Q <- sqrt(sum(wls$coef ^ 2, na.rm = TRUE)) / (1 + lambd) # convert lambda into Q
    }
  }
  
  return(list(Q = Q, lambda = lambd))
}

# Auxiliary function that solves the (un)constrained problem to estimate b
# depending on the desired method
b.est <- function(A, Z, C = NULL, J, KM, w.constr, V, CVXR.solver = "ECOS") {
  
  # Check for global clarabel preference for all constraint types
  use_clarabel <- getOption("SCMs.prefer_clarabel", FALSE)
  
  if (use_clarabel && requireNamespace("clarabel", quietly = TRUE)) {
    # Use clarabel for all constraint types (including pensynth)
    if (w.constr[["p"]] == "pensynth") {
      result <- b.est.clarabel(A, Z, C, J, KM, w.constr, V)
      attr(result, "solver_used") <- "clarabel"
      attr(result, "constraint_type") <- "pensynth"
      return(result)
    } else {
      tryCatch({
        result <- b.est.clarabel.unified(A, Z, C, J, KM, w.constr, V)
        attr(result, "solver_used") <- "clarabel"
        attr(result, "constraint_type") <- w.constr[["p"]]
        return(result)
      }, error = function(e) {
        # If clarabel fails, fall back to CVXR but log the fallback
        warning(paste0("Clarabel failed for constraint '", w.constr[["p"]], "', falling back to CVXR: ", e$message))
        result <- b.est.cvxr.fallback(A, Z, C, J, KM, w.constr, V, CVXR.solver)
        attr(result, "solver_used") <- "CVXR"
        attr(result, "constraint_type") <- w.constr[["p"]]
        attr(result, "fallback_reason") <- e$message
        return(result)
      })
    }
  }
  
  # Fallback to CVXR/specialized functions
  if (w.constr[["p"]] == "pensynth") {
    result <- b.est.cvxr.pensynth(A, Z, C, J, KM, w.constr, V, CVXR.solver)
    attr(result, "solver_used") <- "CVXR"
    attr(result, "constraint_type") <- "pensynth"
    return(result)
  }
  

  # Use CVXR solver (default path)
  result <- b.est.cvxr.fallback(A, Z, C, J, KM, w.constr, V, CVXR.solver)
  attr(result, "solver_used") <- "CVXR"
  attr(result, "constraint_type") <- w.constr[["p"]]
  return(result)
}

#' @title CVXR fallback function for W optimization
#' @description Solve W optimization using CVXR (original implementation)
#' @param A Treated unit features matrix
#' @param Z Donor features matrix
#' @param J Number of donor units
#' @param KM Number of additional variables
#' @param w.constr Constraint specification
#' @param V Feature weighting matrix
#' @param CVXR.solver CVXR solver to use
#' @return Numeric vector of weights
b.est.cvxr.fallback <- function(A, Z, C = NULL, J, KM, w.constr, V, CVXR.solver) {
  dire <- w.constr[["dir"]]
  lb <- w.constr[["lb"]]
  p <- w.constr[["p"]]
  QQ <- w.constr[["Q"]]
  

  if (p == "L1-L2") {
    Q2 <- w.constr[["Q2"]]
  } else {
    Q2 <- NULL
  }
  
  # Determine problem dimensions
  if (!is.null(C)) {
    # Include constant terms - optimize over both w (J+KM weights) and γ (constant coefficients)
    n_const <- ncol(C)
    x <- CVXR::Variable(J + KM + n_const)
    # Split variables: w = x[1:(J+KM)], γ = x[(J+KM+1):(J+KM+n_const)]
    w_vars <- x[1:(J+KM)]
    gamma_vars <- x[(J+KM+1):(J+KM+n_const)]
    objective <- CVXR::Minimize(CVXR::quad_form(A - Z %*% w_vars - C %*% gamma_vars, V))
  } else {
    # No constant terms - original formulation
    x <- CVXR::Variable(J + KM)
    objective <- CVXR::Minimize(CVXR::quad_form(A - Z %*% x, V))
  }
  
  if (p == "no norm") {
    # least squares
    constraints <- list()
    
  } else if (p == "L1") {
    if (dire == "==") {
      # simplex
      constraints <- list(CVXR::sum_entries(x[1:J]) == QQ, x[1:J] >= lb)
    } else if (dire == "<=") {
      # lasso
      #constraints <- list(CVXR::norm1(x[1:J]) <= QQ, x[1:J] >= lb)
      constraints <- list(CVXR::norm1(x[1:J]) <= QQ)
    }
    
  } else if (p == "L2") {
    # ridge
    if (dire == "==") {
      constraints <- list(CVXR::sum_squares(x[1:J]) == CVXR::power(QQ, 2))
    } else if (dire == "<=") {
      constraints <- list(CVXR::sum_squares(x[1:J]) <= CVXR::power(QQ, 2))
    }
    
  } else if (p == "L1-L2") {
    constraints <- list(CVXR::sum_entries(x[1:J]) == QQ,
                        CVXR::power(CVXR::cvxr_norm(x[1:J], 2), 2) <= CVXR::power(Q2, 2),
                        x[1:J] >= lb)
                        
  } else if (p == "pensynth") {
    # This should not be reached due to early return in b.est()
    stop("Pensynth constraint should be handled by specialized functions")
  }
  
  # num_iter required in large lasso/L1-L2/ridge problems is often >10k, which is the default in OSQP/ECOS
  prob <- CVXR::Problem(objective, constraints)
  sol <- CVXR::solve(prob, solver = CVXR.solver, num_iter = 100000L, verbose = FALSE)
  
  
  if (is.null(sol$getValue)) {
    stop("Solver failed - getValue method not available")
  }
  
  b <- sol$getValue(x)
  alert <- !(sol$status %in% c("optimal", "optimal_inaccurate"))
  
  if (alert == TRUE) {
    stop(paste0("Estimation algorithm not converged! The algorithm returned the value:",
                sol$status, ". To check to what errors it corresponds go to
               'https://cvxr.rbind.io/cvxr_examples/cvxr_gentle-intro/'. Typically, this occurs
                because the problem is badly-scaled. If so, scaling the data fixes the issue. Another
                fix could be changing the algorithm via the option 'solver'. Check your available options
                using CVXR::installed_solvers() and consult
                'https://cvxr.rbind.io/cvxr_examples/cvxr_using-other-solvers/'"),
         call. = FALSE)
  }
  
  b <- b[, 1, drop = TRUE]
  
  if (!is.null(C)) {
    # When C is present, return both donor weights and constant coefficients
    n_const <- ncol(C)
    donor_weights <- b[1:(J+KM)]
    constant_coeffs <- b[(J+KM+1):(J+KM+n_const)]
    
    names(donor_weights) <- colnames(Z)
    names(constant_coeffs) <- colnames(C)
    
    # Return list with both components
    return(list(w = donor_weights, gamma = constant_coeffs))
  } else {
    # Original case - return only donor weights
    names(b) <- colnames(Z)
    return(b)
  }
}

# CVXR implementation for pensynth (current)
b.est.cvxr.pensynth <- function(A, Z, C = NULL, J, KM, w.constr, V, CVXR.solver) {
  lambda <- w.constr[["lambda"]]
  
  # Extract V weights for penalty calculation  
  v_weights <- diag(V)
  
  # Calculate pensynth penalty components
  # Z is (K x J), v_weights is length K, A is length K
  sqrt_v <- sqrt(v_weights)
  Z_weighted <- sqrt_v * Z  # Broadcasting: (K x 1) * (K x J) = (K x J)
  A_weighted <- sqrt_v * A  # Broadcasting: (K x 1) * (K x 1) = (K x 1)
  
  # Vectorized computation of ||X1v - X0v_j||^2 for all donors j
  # A_weighted is length K, Z_weighted is (K x J)
  # Use sweep() to subtract A_weighted from each column of Z_weighted, then compute column-wise squared norms
  Delta <- colSums(sweep(Z_weighted, 1, A_weighted, "-")^2)
  
  # Set up CVXR optimization
  x <- CVXR::Variable(J + KM)
  
  # Combined objective: standard SCM fit + pensynth penalty
  base_objective <- CVXR::quad_form(A - Z %*% x, V)
  penalty_term <- lambda * sum(Delta * x[1:J])
  objective <- CVXR::Minimize(base_objective + penalty_term)
  
  # Standard simplex constraints (sum to 1, non-negative)
  constraints <- list(CVXR::sum_entries(x[1:J]) == 1, x[1:J] >= 0)
  
  # Solve optimization problem
  prob <- CVXR::Problem(objective, constraints)
  sol <- CVXR::solve(prob, solver = CVXR.solver, num_iter = 100000L, verbose = FALSE)
  
  if (is.null(sol$getValue)) {
    stop("Solver failed - getValue method not available")
  }
  
  if (!(sol$status %in% c("optimal", "optimal_inaccurate"))) {
    stop(paste0("Pensynth estimation not converged! Status: ", sol$status))
  }
  
  b <- sol$getValue(x)
  b <- b[, 1, drop = TRUE]
  names(b) <- colnames(Z)
  
  return(b)
}

# Clarabel implementation for pensynth (future)
b.est.clarabel <- function(A, Z, C = NULL, J, KM, w.constr, V) {
  # High-performance clarabel implementation for pensynth
  
  # Only handle pensynth constraint for now
  if (w.constr[["name"]] != "pensynth") {
    warning("Clarabel implementation currently only supports pensynth constraint")
    return(b.est.cvxr(A, Z, J, KM, w.constr, V, "ECOS"))
  }
  
  lambda <- w.constr[["lambda"]]
  
  # Extract V weights for penalty calculation  
  v_weights <- diag(V)
  sqrt_v <- sqrt(v_weights)
  Z_weighted <- sqrt_v * Z  # Broadcasting: (K x 1) * (K x J) = (K x J)
  A_weighted <- sqrt_v * A  # Broadcasting: (K x 1) * (K x 1) = (K x 1)
  
  # Vectorized computation of ||X1v - X0v_j||^2 for all donors j
  Delta <- colSums(sweep(Z_weighted, 1, A_weighted, "-")^2)
  
  # Set up QP formulation following original pensynth implementation
  # P = X0VX0 = crossprod(X0v), q = -X1VX0 + lambda*Delta
  P <- crossprod(Z_weighted)  # X0v'X0v (J x J)
  q <- -crossprod(A_weighted, Z_weighted) + lambda * Delta  # -X1v'X0v + lambda*Delta (length J)
  
  # Add small regularization to ensure P is positive definite
  P <- P + 1e-8 * diag(J)
  
  # Set up constraints: sum(x) = 1, x >= 0
  # For clarabel: Ax + s = b, where s ∈ K
  A_constraints <- rbind(
    matrix(1, nrow = 1, ncol = J),    # Sum constraint: sum(x) = 1
    -diag(J)                          # Non-negativity: -x + s = 0, s >= 0
  )
  b_constraints <- c(1, rep(0, J))    # RHS: [1, 0, 0, ..., 0]
  
  # Call clarabel solver with high-precision settings
  control_settings <- clarabel::clarabel_control(
    verbose = FALSE,
    max_iter = 1000L,
    tol_gap_abs = 1e-8,
    tol_gap_rel = 1e-8
  )
  
  result <- clarabel::clarabel(
    A = A_constraints,
    b = b_constraints,
    q = as.vector(q),
    P = P,
    cones = list(z = 1L, l = J),  # 1 zero cone (equality) + J linear cones (non-negativity)
    control = control_settings
  )
  
  # Check clarabel status - fail hard if not solved (numeric code 2 = "Solved")
  if (result$status != 2) {
    status_desc <- names(clarabel::solver_status_descriptions())[result$status]
    if (result$status == 3) {  # PrimalInfeasible
      stop(paste0("Pensynth optimization problem is primal infeasible! ", 
                  "The constraints cannot be satisfied. This may indicate: ",
                  "(1) lambda value is too large, (2) data contains invalid values, ",
                  "or (3) problem setup is incorrect. Status: ", result$status, " (", status_desc, ")"))
    } else {
      stop(paste0("Clarabel pensynth estimation failed! Status: ", result$status, 
                  " (", status_desc, "). Check problem setup, lambda value, or disable clarabel with options(SCMs.prefer_clarabel = FALSE)"))
    }
  }
  
  weights <- result$x
  names(weights) <- colnames(Z)
  
  return(weights)
}

#' @title Clarabel W Optimization - OLS (Unconstrained)
#' @description Solve unconstrained quadratic programming problem using clarabel
#' @param A Treated unit features matrix
#' @param Z Donor features matrix  
#' @param J Number of donor units
#' @param KM Number of additional variables
#' @param w.constr Constraint specification
#' @param V Feature weighting matrix
#' @return Numeric vector of weights
b.est.clarabel.ols <- function(A, Z, C = NULL, J, KM, w.constr, V) {
  # QP formulation: minimize 0.5 * x'Px + q'x
  if (!is.null(C)) {
    # With constant terms: minimize ||A - Z*w - C*γ||_V^2
    # Combined matrix [Z C], variables [w; γ]
    ZC <- cbind(Z, C)
    P <- 2 * t(ZC) %*% V %*% ZC
    q <- -2 * as.vector(t(ZC) %*% V %*% A)
    n_vars <- J + KM + ncol(C)
  } else {
    # Original case: minimize ||A - Z*w||_V^2
    P <- 2 * t(Z) %*% V %*% Z
    q <- -2 * as.vector(t(Z) %*% V %*% A)
    n_vars <- J + KM
  }
  
  # For unconstrained QP, clarabel might not accept empty constraint matrices
  # Let's add a dummy constraint that doesn't restrict anything: 0*x = 0
  A_constraints <- matrix(0, nrow = 1, ncol = n_vars)  # All zeros
  b_constraints <- 0  # RHS = 0
  
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
    cones = list(z = 1L),  # One zero cone for the dummy constraint
    control = control_settings
  )
  
  # Check clarabel status - accept only solved status (numeric code 2 = "Solved")
  if (result$status != 2) {
    status_desc <- names(clarabel::solver_status_descriptions())[result$status]
    stop(paste0("Clarabel OLS optimization failed! Status: ", result$status, " (", status_desc, ")"))
  }
  
  solution <- as.numeric(result$x)
  
  if (!is.null(C)) {
    # When C is present, return both donor weights and constant coefficients
    n_const <- ncol(C)
    donor_weights <- solution[1:(J+KM)]
    constant_coeffs <- solution[(J+KM+1):(J+KM+n_const)]
    
    names(donor_weights) <- colnames(Z)
    names(constant_coeffs) <- colnames(C)
    
    # Return list with both components
    return(list(w = donor_weights, gamma = constant_coeffs))
  } else {
    # Original case - return only donor weights
    names(solution) <- colnames(Z)
    return(solution)
  }
}

#' @title Clarabel W Optimization - Simplex
#' @description Solve simplex-constrained QP using clarabel (sum to 1, non-negative)
#' @param A Treated unit features matrix
#' @param Z Donor features matrix
#' @param J Number of donor units
#' @param KM Number of additional variables
#' @param w.constr Constraint specification
#' @param V Feature weighting matrix
#' @return Numeric vector of weights
b.est.clarabel.simplex <- function(A, Z, C = NULL, J, KM, w.constr, V) {
  # QP formulation: minimize 0.5 * x'Px + q'x
  if (!is.null(C)) {
    # With constant terms: minimize ||A - Z*w - C*γ||_V^2
    # Combined matrix [Z C], variables [w; γ]
    ZC <- cbind(Z, C)
    P <- 2 * t(ZC) %*% V %*% ZC
    q <- -2 * as.vector(t(ZC) %*% V %*% A)
    n_vars <- J + KM + ncol(C)
  } else {
    # Original case: minimize ||A - Z*w||_V^2
    P <- 2 * t(Z) %*% V %*% Z
    q <- -2 * as.vector(t(Z) %*% V %*% A)
    n_vars <- J + KM
  }
  
  # Get constraint parameters
  QQ <- w.constr[["Q"]]  # Should be 1 for simplex
  lb <- w.constr[["lb"]]  # Lower bound, typically 0
  
  # Constraint matrices for simplex: sum(x[1:J]) = QQ, x[1:J] >= lb
  # A_constraints * x + s = b_constraints, s in cones
  # Note: constraints only apply to donor weights (first J variables), not constant coefficients
  
  # Equality constraint: sum(x[1:J]) = QQ
  A_eq <- matrix(0, nrow = 1, ncol = n_vars)
  A_eq[1, 1:J] <- 1
  
  # Inequality constraints: x[1:J] >= lb, reformulated as -x[1:J] + s = -lb, s >= 0
  A_ineq <- matrix(0, nrow = J, ncol = n_vars)
  A_ineq[1:J, 1:J] <- -diag(J)
  
  # Combine constraints
  A_constraints <- rbind(A_eq, A_ineq)
  b_constraints <- c(QQ, rep(-lb, J))
  
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
    cones = list(z = 1L, l = J),  # 1 zero cone (equality) + J linear cones (non-negativity)
    control = control_settings
  )
  
  # Check clarabel status - accept only solved status (numeric code 2 = "Solved")
  if (result$status != 2) {
    status_desc <- names(clarabel::solver_status_descriptions())[result$status]
    stop(paste0("Clarabel simplex optimization failed! Status: ", result$status, " (", status_desc, ")"))
  }
  
  solution <- as.numeric(result$x)
  
  if (!is.null(C)) {
    # When C is present, return both donor weights and constant coefficients
    n_const <- ncol(C)
    donor_weights <- solution[1:(J+KM)]
    constant_coeffs <- solution[(J+KM+1):(J+KM+n_const)]
    
    names(donor_weights) <- colnames(Z)
    names(constant_coeffs) <- colnames(C)
    
    # Return list with both components
    return(list(w = donor_weights, gamma = constant_coeffs))
  } else {
    # Original case - return only donor weights
    names(solution) <- colnames(Z)
    return(solution)
  }
}

#' @title Unified Clarabel W Optimization Dispatcher
#' @description Dispatch to appropriate clarabel constraint-specific function
#' @param A Treated unit features matrix
#' @param Z Donor features matrix
#' @param J Number of donor units
#' @param KM Number of additional variables
#' @param w.constr Constraint specification
#' @param V Feature weighting matrix
#' @return Numeric vector of weights
b.est.clarabel.unified <- function(A, Z, C = NULL, J, KM, w.constr, V) {
  # Extract constraint parameters
  p <- w.constr[["p"]]
  dire <- w.constr[["dir"]]
  
  # Dispatch to appropriate constraint-specific function
  if (p == "no norm") {
    return(b.est.clarabel.ols(A, Z, C, J, KM, w.constr, V))
    
  } else if (p == "L1" && dire == "==") {
    return(b.est.clarabel.simplex(A, Z, C, J, KM, w.constr, V))
    
  } else if (p == "L1" && dire == "<=") {
    # Lasso - not implemented yet
    stop("Clarabel lasso constraint not yet implemented")
    
  } else if (p == "L2") {
    # Ridge constraint - temporarily not implemented
    stop("Ridge constraint with clarabel not yet fully implemented. Use CVXR by setting options(SCMs.prefer_clarabel = FALSE)")
    
  } else if (p == "L1-L2") {
    # L1-L2 - not implemented yet
    stop("Clarabel L1-L2 constraint not yet implemented")
    
  } else if (p == "pensynth") {
    # Should not reach here - pensynth handled separately
    stop("Pensynth constraint should be handled by specialized functions")
    
  } else {
    stop(paste0("Unknown constraint type: p=", p, ", dir=", dire))
  }
}


V.prep <- function(type, B, T0.features, I) {
  if (type == "separate") {
    # Default (separate fit)
    V <- diag(dim(B)[1])
    
  } else if (type == "pooled") {
    
    dim.V <- unlist(lapply(T0.features, function(x) sum(unlist(x))))
    max.dim <- max(dim.V)
    ones <- matrix(1, nrow = I, ncol = 1)
    eye <- diag(1, nrow = max.dim, ncol = max.dim)
    V <- kronecker(ones %*% t(ones), eye) # structure if T0*M was balanced across treated unit
    
    sel <- matrix(TRUE, nrow = nrow(V), ncol = ncol(V))
    for (i in seq_len(I)) {
      # trim V according to length of pre-treatment period
      if (dim.V[i] < max.dim) {
        shift <- (i - 1) * max.dim
        sel[(shift + dim.V[i] + 1):(shift + max.dim),] <- FALSE
        sel[, (shift + dim.V[i] + 1):(shift + max.dim)] <- FALSE
      }
    }
    row.trim <- rowSums(sel) != 0
    col.trim <- colSums(sel) != 0
    
    V <- V[row.trim, col.trim] / I ^ 2
  }
  
  rownames(V) <- rownames(B)
  colnames(V) <- rownames(B)
  
  return(V)
}



mat2list <- function(mat, cols = TRUE) {
  # select rows
  
  names <- strsplit(rownames(mat), "\\.")
  rnames <- unlist(lapply(names, "[[", 1))
  tr.units <- unique(rnames)
  
  # select columns
  matlist <- list()
  if (cols == TRUE) {
    if (ncol(mat) > 1) {
      
      names <- strsplit(colnames(mat), "\\.")
      
      cnames <- unlist(lapply(names, "[[", 1))
      
      for (tr in tr.units) {
        matlist[[tr]] <- mat[rnames == tr, cnames == tr, drop = FALSE]
      }
    } else if (ncol(mat) == 1) {
      
      for (tr in tr.units) {
        matlist[[tr]] <- mat[rnames == tr, 1, drop = FALSE]
      }
    } else {
      
      for (tr in tr.units) {
        matlist[[tr]] <- mat[rnames == tr, 0, drop = FALSE]
      }
    }
  } else if (cols == FALSE) {
    for (tr in tr.units) {
      matlist[[tr]] <- mat[rnames == tr,, drop = FALSE]
    }
  }
  
  return(matlist)
}


most_similar <- function(
    dataset,
    outcome,
    col_name_unit_name,
    name_treated_unit,
    treated_period,
    min_period, 
    col_name_period){
  
  dataset = data.table::as.data.table(dataset)
  
  dataset_pre_treatment = dataset[between(get(col_name_period), min_period, (treated_period-1))]
  dataset_pre_treatment_control = dataset_pre_treatment[get(col_name_unit_name)!=name_treated_unit]
  dataset_pre_treatment_treated = dataset_pre_treatment[get(col_name_unit_name)==name_treated_unit] 
  
  ybar_treated = mean(dataset_pre_treatment_treated[[outcome]])
  ybar_controls = dataset_pre_treatment_control[, .(ybar_control=mean(get(outcome))), by=col_name_unit_name]
  ybar_controls[, diff_with_treated := abs(ybar_treated - ybar_control)]
  ybar_controls = ybar_controls[order(diff_with_treated)]
  control_units = ybar_controls[ 1 : round(nrow(ybar_controls) / 2) ][[col_name_unit_name]]
  dataset = dataset[get(col_name_unit_name)%in%c(name_treated_unit, control_units)]
  return(dataset)
}
