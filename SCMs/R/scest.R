#' @title Synthetic Control Estimation for Single Treated Unit
#' @description This function performs synthetic control estimation for a single treated unit.
#'
#' @param data A list of class "scdata" containing the prepared data for synthetic control estimation.
#' @param w.constr A list specifying constraints on the weights. Default is NULL.
#' @param feature_weights Character. Method for assigning weights to features. Can be "uniform" or "optimize". Default is "uniform".
#' @param V Character or matrix. Specifies the weighting matrix. Default is "separate".
#' @param V.mat Matrix. Custom weighting matrix. Default is NULL.
#' @param solver Character. The solver to use for optimization. Default is "ECOS".
#' @param save.data Not used in the current implementation.
#'
#' @return A list of class "scest" containing the estimation results and data.
#'
#' @import methods
#' @import CVXR
#' @import optimx
#'
#' @export
#'
#' @examples
#' # Example usage (replace with actual example when available)
#' # result <- scest(data = my_scdata, w.constr = list(name = "simplex"), feature_weights = "optimize")
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
  if (out.in.features == TRUE) {
    fit.pre <- A.hat[1:T0.features[outcome.var], , drop = FALSE]
    names <- strsplit(rownames(fit.pre), "\\.")
    rownames(fit.pre) <- unlist(lapply(names, function(x) paste(x[1], x[3], sep=".")))
  } else {
    fit.pre <- Y.donors %*% w
  }

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
