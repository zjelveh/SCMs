#' @title Synthetic Control Estimation for Single Treated Unit
#' @description This function estimates synthetic control weights for a single treated unit using various
#' constraint specifications and optimization methods. It implements the core synthetic control methodology
#' with flexible constraint options and feature weighting schemes.
#'
#' @param data A list of class "scdata" containing the prepared data from \code{scdata()}.
#' @param w.constr List specifying constraints on the synthetic control weights. Options include:
#'   \itemize{
#'     \item \code{list(name = "simplex")} - Standard simplex constraint (weights sum to 1, non-negative)
#'     \item \code{list(name = "lasso", Q = 0.1)} - L1 penalty with regularization parameter Q
#'     \item \code{list(name = "ridge", Q = 0.1)} - L2 penalty with regularization parameter Q  
#'     \item \code{list(name = "ols")} - Ordinary least squares (no constraints)
#'     \item \code{list(name = "L1-L2", Q1 = 0.1, Q2 = 0.1)} - Combined L1 and L2 penalties
#'   }
#'   Default is NULL (no constraints).
#' @param feature_weights Character. Method for weighting pre-treatment periods/features:
#'   \itemize{
#'     \item \code{"uniform"} - Equal weights for all periods (default)
#'     \item \code{"optimize"} - Optimize feature weights using nested optimization
#'   }
#' @param V Character or matrix specifying the feature weighting approach:
#'   \itemize{
#'     \item \code{"separate"} - Optimize V matrix separately (default)
#'     \item \code{"nested"} - Nested optimization of V and W matrices
#'     \item Custom matrix - User-provided V matrix
#'   }
#' @param V.mat Matrix. Custom V (feature weighting) matrix. If provided, overrides V parameter.
#' @param solver Character. CVXR solver to use ("ECOS", "OSQP", "GUROBI", etc.). 
#'   Must be installed on system. Default is "ECOS".
#' @param save.data Logical. Whether to save input data in output object (for debugging). Default is NULL.
#'
#' @details
#' The function implements the synthetic control method by solving the optimization problem:
#' \deqn{\min_W ||X_1 - X_0 W||_V^2 + penalty(W)}
#' 
#' where \eqn{X_1} is the treated unit's pre-treatment characteristics, \eqn{X_0} is the control units' 
#' characteristics matrix, \eqn{W} are the synthetic control weights, and \eqn{V} is the feature weighting matrix.
#' 
#' Different constraint types impose different penalty functions:
#' \itemize{
#'   \item Simplex: \eqn{\sum W_i = 1, W_i \geq 0}
#'   \item Lasso: \eqn{Q \sum |W_i|}  
#'   \item Ridge: \eqn{Q \sum W_i^2}
#'   \item OLS: No constraints
#' }
#'
#' @return A list of class "scest" containing:
#' \itemize{
#'   \item \code{w.weights} - Vector of synthetic control weights
#'   \item \code{Y.pre.fit} - Fitted pre-treatment outcomes
#'   \item \code{Y.post.fit} - Fitted post-treatment outcomes  
#'   \item \code{loss.v} - Pre-treatment fit loss
#'   \item \code{specs} - Estimation specifications
#'   \item \code{data} - Original data (if save.data = TRUE)
#'   \item \code{V.mat} - Feature weighting matrix used
#'   \item Additional estimation metadata
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic simplex-constrained synthetic control
#' result <- scest(data = scm_data, 
#'                 w.constr = list(name = "simplex"))
#'
#' # Lasso-regularized synthetic control  
#' result_lasso <- scest(data = scm_data,
#'                       w.constr = list(name = "lasso", Q = 0.1))
#'
#' # With optimized feature weights
#' result_opt <- scest(data = scm_data,
#'                     w.constr = list(name = "simplex"),
#'                     feature_weights = "optimize")
#'                     
#' # Using custom V matrix
#' custom_V <- diag(nrow(scm_data$Y.pre))
#' result_custom <- scest(data = scm_data,
#'                        w.constr = list(name = "simplex"),
#'                        V.mat = custom_V)
#' }
scest <- function(data,
                  w.constr  = NULL,
                  feature_weights = 'uniform',
                  V         = "separate",
                  V.mat     = NULL,
                  solver    = "ECOS",
                  save.data = NULL) {

  # Check if input data is of correct class
  if (!methods::is(data, "scdata")) {
    stop("data should be the object returned by running scdata!")
  }

  # Check and process weight constraints
  if (!is.null(w.constr)) {
    if (!is.list(w.constr)) {
      stop("w.constr should be a list!")
    }
    if (!"name" %in% names(w.constr)) {
      w.constr[["name"]] <- "NULL"
    } else if (!w.constr[["name"]] %in% c("simplex","lasso","ridge","ols","L1-L2")) {
      stop("If 'name' is specified in w.constr, it should be 'simplex', 'lasso', 'ridge', 'ols', or 'L1-L2'.")
    }
  }

  # Check if specified solver is available
  if (!(solver %in% CVXR::installed_solvers())) {
    stop(paste0("The specified solver - ", solver," - is not available on your machine!"))
  }

  # Extract data matrices and specifications
  A <- data$A
  B <- data$B
  P <- data$P
  Z <- copy(B)
  X0sd <- data$X0_sds
  Y.donors <- data$Y.donors
  outcome.var <- data$specs$outcome.var
  Z0 <- Y.donors
  Z1 <- data$Y.pre

  class.type <- 'scpi_data'
  V.type <- V

  # Extract data specifications
  J <- data$specs$J
  KM <- data$specs$KM
  I <- 1
  KMI <- KM
  Jtot <- J
  M <- data$specs$M
  T0.features <- data$specs$T0.features
  out.in.features <- data$specs$out.in.features

  # Prepare weighting matrix V
  if (is.character(V) == FALSE) {
    stop("The object V should be a string! If you want to input a matrix use the option V.mat!")
  }


  if (!is.null(V.mat)) {
    if (!is.matrix(V.mat)) {
      stop("The object V.mat should be a matrix!")
    }

    if (nrow(V.mat) != nrow(B) || ncol(V.mat) != nrow(B)) {
      stop(paste0("V.mat should be a ", nrow(B), "x", nrow(B)," matrix, but currently it is
                a ", nrow(V.mat), "x", ncol(V.mat), " matrix!"))
    }
  } else if(feature_weights=='optimize'){
    # Optimize feature weights

    require(optimx)
    SV1 = rep(1/nrow(B), nrow(B))

    Margin.ipop = 0.05
    Sigf.ipop = 5
    Bound.ipop = 10
    all.methods = FALSE
    optimxmethod = c("Nelder-Mead","BFGS")


    rgV.optim.1 <- optimx(par=SV1, fn=fn.V,
                          gr=NULL,
                          hess=NULL,
                          method=optimxmethod,
                          hessian=FALSE,
                          itnmax = NULL,# c(10,10),
                          control=list(kkt=FALSE,
                                       starttests=FALSE,
                                       dowarn=FALSE,
                                       all.methods=all.methods),
                          X0.scaled = B,
                          X1.scaled = A,
                          Z0 = Z0,
                          Z1 = Z1,
                          margin.ipop = Margin.ipop,
                          sigf.ipop = Sigf.ipop,
                          bound.ipop = Bound.ipop
    )

    # get minimum
    # if(verbose==TRUE){print(rgV.optim.1)}
    rgV.optim <- collect.optimx(rgV.optim.1, "min")

    solution.v   <- abs(rgV.optim$par)/sum(abs(rgV.optim$par))
    V.mat = diag(solution.v)
  } else{
    if(feature_weights=='uniform'){
      V.type = 'separate'
    }
    V.mat <- V.prep(type = V.type, B, T0.features, I)
  }
  V <- V.mat

 # Handle NA values in V
  if (mean(is.na(diag(V)) > 0)) {
    V <- diag(rep(1/length(nrow(B)), (nrow(B))))
  }
  rownames(V) <- rownames(Z)
  colnames(V) <- rownames(V)

  # Estimate synthetic control
  w.constr <- w.constr.OBJ(w.constr, A, Z, V, J, KM, M)
  if (w.constr[["name"]] == "lasso") solver <- "OSQP"
  b <- b.est(A = A, Z = Z, J = J, KM = KM, w.constr = w.constr, V = V, CVXR.solver = solver)

  # Process results
  if (KMI == 0) {
    w <- b
    r <- NULL
  } else {
    w <- b[1:Jtot]
    r <- b[(Jtot + 1):length(b)]
  }

  # Calculate fitted values and residuals
  A.hat <- Z %*% b
  res <- A - A.hat
  Z_scaled <- Z * X0sd
  A.hat <- Z_scaled %*% b

  # Prepare pre-treatment fit
  # ALWAYS use time-varying predictions regardless of matching features
  # This ensures inference gets the time series it expects
  
  # Ensure w is a proper column matrix for matrix multiplication
  if (is.vector(w)) {
    w <- matrix(w, ncol = 1)
  }
  fit.pre <- Y.donors %*% w

  # Post-treatment prediction
  fit.post <- P %*% b

  # Prepare estimation results
  est.results <- list(b = b, w = w, r = r, Y.pre.fit = fit.pre, Y.post.fit = fit.post,
                      A.hat = A.hat, res = res, V = V, w.constr = w.constr)

  # Prepare data frame to return
  df <- list(A = data$A, B = data$B, P = data$P, Z = Z, specs = data$specs,
             Y.pre = data$Y.pre, Y.post = data$Y.post,
             Y.pre.agg = data$Y.pre, Y.post.agg = data$Y.post.agg)

  # Prepare return object
  to.return <- list(data = df, est.results = est.results)
  class(to.return) <- 'scest'
  to.return$data$specs$class.type <- 'scpi_scest'

  return(to.return)
}
