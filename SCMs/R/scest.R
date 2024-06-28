scest <- function(data,
                  w.constr  = NULL,
                  feature_weights = 'uniform',
                  V         = "separate",
                  V.mat     = NULL,
                  solver    = "ECOS",
                  plot      = FALSE,
                  plot.name = NULL,
                  plot.path = NULL,
                  save.data = NULL) {

  ##########################
  if ((methods::is(data, "scdata") || methods::is(data, "scdataMulti")) == FALSE) {
    stop("data should be the object returned by running scdata or scdata_multi!")
  }

  if (is.null(w.constr) == FALSE) { # The user has specified W

    if (is.list(w.constr) == FALSE) {
      stop("w.constr should be a list!")
    }

    if (!"name" %in% names(w.constr)) {
      w.constr[["name"]] <- "NULL"
    } else {
      if (!w.constr[["name"]] %in% c("simplex","lasso","ridge","ols","L1-L2")) {
        stop("If 'name' is specified in w.constr, then it should be chosen among
             'simplex', 'lasso', 'ridge', 'ols', 'L1-L2'.")
      }
    }
  }

  if (!(solver %in% CVXR::installed_solvers())) {
    stop(paste0("The specified solver - ", solver," - is not available on your machine! Run
                CVXR::installed_solvers() to see the list of available options on your machine."))
  }

  # Data matrices
  A <- data$A
  B <- data$B
  P <- data$P
  Z <- copy(B)
  X0sd <- data$X0_sds
  Y.donors <- data$Y.donors
  outcome.var <- data$specs$outcome.var
  Z0 <- Y.donors
  Z1 <- data$Y.pre

  if (class(data)[1] == 'scdata') {
    class.type <- 'scpi_data'
  } else if (class(data)[1] == 'scdataMulti') {
    class.type <- 'scpi_data_multi'
  }

  V.type <- V

  # Data specs
  if (class.type == 'scpi_data') {
    J               <- data$specs$J
    KM              <- data$specs$KM
    I               <- 1
    KMI             <- KM
    Jtot            <- J
    M               <- data$specs$M

  } else if (class.type == 'scpi_data_multi') {
    J               <- data$specs$J # total number of donors
    Jtot            <- sum(unlist(J))
    KM              <- data$specs$KM # total number of covs used for adjustment per treated unit
    KMI             <- data$specs$KMI # total number of covariates used for adjustment
    I               <- data$specs$I  # number of treated units
    T0.M            <- lapply(data$specs$T0.features, sum) # observations per treated unit
    M               <- data$specs$M
  }

  T0.features     <- data$specs$T0.features
  out.in.features <- data$specs$out.in.features

  ##########################
  ## Set up the estimation problem

  # Create weighting matrix
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

  if(mean(is.na(diag(V))>0)){
    V <- diag(rep(1/length(nrow(B)), (nrow(B))))
  }
  rownames(V) <- rownames(Z)
  colnames(V) <- rownames(V)
  #############################

  # Create lists of matrices
  # Estimate SC
  if (class.type == 'scpi_data') { # single treated unit

    w.constr <- w.constr.OBJ(w.constr, A, Z, V, J, KM, M)
    if (w.constr[["name"]] == "lasso") solver <- "OSQP"

    b <- b.est(A = A, Z = Z, J = J, KM = KM, w.constr = w.constr, V = V, CVXR.solver = solver)

  } else if (class.type == 'scpi_data_multi') { # multiple treated units

    A.list <- mat2list(A)
    B.list <- mat2list(B)
    V.list <- mat2list(V)
    w.constr.list <- list()

    w.store <- c()
    r.store <- c()
    j.lb <- 1
    j.ub <- J[[1]]

    for (i in seq_len(I)) {
      if (i > 1){
        j.lb <- j.ub + 1
        j.ub <- j.lb + J[[i]] - 1
      }

      A.i <- A.list[[i]]
      Z.i <- cbind(B.list[[i]])
      V.i <- V.list[[i]]

      w.constr.list[[data$specs$treated.units[i]]] <- w.constr.OBJ(w.constr, A.i, Z.i, V.i, J[[i]], KM[[i]], M[[i]])
      if (w.constr.list[[i]]["name"] == "lasso") solver <- "OSQP"

      if (V.type == "separate") {
        res <- b.est(A = A.i, Z = Z.i, J = J[[i]], KM = KM[[i]], w.constr = w.constr.list[[i]], V = V.i, CVXR.solver = solver)
        w.store <- c(w.store, res[1:J[[i]]])
        if (KM[[i]] > 0) r.store <- c(r.store, res[(J[[i]]+1):length(res)])
      }
    }

    if (V.type != "separate" ) {
      if (w.constr.list[[i]]["name"] == "lasso") solver <- "OSQP"
      b <- b.est.multi(A = A, Z = Z, J = J, KMI = KMI, I = I,
                       w.constr = w.constr.list, V = V, CVXR.solver = solver)
      b <- b[,1,drop=TRUE]

    } else if (V.type == "separate") {
      b <- c(w.store, r.store)
    }
    w.constr <- w.constr.list
  }

  ##########################
  ## Create useful objects
  # Store point estimates
  if (KMI == 0) {
    w <- b                       # Estimated weights
    r <- NULL                    # Loading factors of additional covariates
  } else {
    w <- b[1:Jtot]
    r <- b[(Jtot + 1):length(b)]
  }

  # Fitted values and residuals
  A.hat    <- Z %*% b                            # Pre-treatment fit of synthetic unit
  res      <- A  -  A.hat                        # Estimation residual u

  # Reflate A.hat

  Z_scaled = Z * X0sd
  A.hat    <- Z_scaled %*% b                            # Pre-treatment fit of synthetic unit

  # Pre-treatment fit of outcome of interest
  if (class.type == 'scpi_data'){
    if (out.in.features == TRUE) {
      fit.pre  <- A.hat[1:T0.features[outcome.var], , drop = FALSE]
      names    <- strsplit(rownames(fit.pre), "\\.")
      rownames(fit.pre) <- unlist(lapply(names, function(x) paste(x[1], x[3],sep=".")))
    } else {
      fit.pre  <- Y.donors %*% w
    }
  } else if(class.type == 'scpi_data_multi') {

    i.lb <- 1
    fit.pre <- c()
    Yd.list <- mat2list(data$Y.donors)
    w.list <- mat2list(as.matrix(w), cols=FALSE)

    for (i in seq_len(I)) {
      i.ub <- i.lb + T0.features[[i]][outcome.var] - 1
      if (out.in.features[[i]] == TRUE) {
        fit.pre.i  <- A.hat[i.lb:i.ub, , drop = FALSE]
        names    <- strsplit(rownames(fit.pre.i), "\\.")
        rownames(fit.pre.i) <- unlist(lapply(names, function(x) paste(x[1],x[3],sep=".")))
      } else {
        fit.pre.i <- Yd.list[[i]] %*% w.list[[i]]
      }
      i.lb <- i.lb + sum(unlist(T0.features[[i]]), na.rm = TRUE)
      fit.pre <- rbind(fit.pre, fit.pre.i)
    }
  }


  # Post-treatment prediction of outcome of interest
  fit.post <- P %*% b

  # Name V
  rownames(V) = rownames(B)
  colnames(V) = rownames(V)

  est.results <- list(b = b,
                      w = w,
                      r = r,
                      Y.pre.fit = fit.pre,
                      Y.post.fit = fit.post,
                      A.hat = A.hat,
                      res = res,
                      V = V,
                      w.constr = w.constr)

  if (class.type == 'scpi_data') {
    df   <- list(A = data$A,
                 B = data$B,
                 P = data$P,
                 Z = Z,
                 specs  = data$specs,
                 Y.pre  = data$Y.pre,
                 Y.post = data$Y.post,
                 Y.pre.agg = data$Y.pre,
                 Y.post.agg = data$Y.post.agg)

  } else if (class.type == 'scpi_data_multi') {

    # shortcut to avoid "no visible binding for global variable 'X' when checking the package
    Treatment <- NULL

    # Y.pre and Y.post might require some extra work
    # if the predictand of interest is aggregate (either over units or over time)
    # The next function process the data in the same way it's done in scplotMulti
    treated.units   <- data$specs$treated.units
    sparse.matrices <- data$specs$sparse.matrices
    anticipation    <- data$specs$anticipation
    period.post     <- data$specs$period.post
    units.est       <- data$specs$units.est
    effect          <- data$specs$effect
    Y.df            <- data$Y.df
    Y.pre.fit       <- fit.pre
    Y.post.fit      <- fit.post

    # create to plot object
    Yprocessed <- outcomeGet(Y.pre.fit=Y.pre.fit, Y.post.fit=Y.post.fit, Y.df=Y.df,
                             units.est=units.est, treated.units=treated.units, plot.type=effect,
                             anticipation=anticipation, period.post=period.post,
                             sparse.matrices=sparse.matrices)

    Ydf.pre <- subset(Yprocessed$toplot, Treatment == 0)
    Ydf.post <- subset(Yprocessed$toplot, Treatment == 1)
    Y.pre.agg  <- as.matrix(Ydf.pre[["Actual"]])
    Y.post.agg <- as.matrix(Ydf.post[["Actual"]])
    names.pre <- paste(Ydf.pre$ID, Ydf.pre$Time, sep=".")
    names.post <- paste(Ydf.post$ID, Ydf.post$Time, sep=".")
    rownames(Y.pre.agg) <- names.pre
    rownames(Y.post.agg) <- names.post
    colnames(Y.pre.agg) <- data$specs$outcome.var
    colnames(Y.post.agg) <- data$specs$outcome.var

    df   <- list(A = data$A,
                 B = data$B,
                 P = data$P,
                 P.diff = data$P.diff,
                 Y.df = data$Y.df,
                 Y.pre = data$Y.pre,
                 Y.post = data$Y.post,
                 Y.pre.agg = Y.pre.agg,
                 Y.post.agg = Y.post.agg,
                 Z = Z,
                 specs  = data$specs)
  }

  ##########################
  ## Return to the user
  to.return <- list(data = df, est.results = est.results)
  class(to.return) <- 'scest'
  if (class.type == 'scpi_data') {
    to.return$data$specs$class.type <- 'scpi_scest'
  } else if (class.type == 'scpi_data_multi') {
    to.return$data$specs$class.type <- 'scpi_scest_multi'
  }

  ##################################################
  ## Plot
  if (plot == TRUE) {
    if (is.null(plot.name) == FALSE) {
      fig.name <- plot.name
    } else {
      fig.name <- "scest_default_plot"
    }

    if (is.null(plot.path) == FALSE) {
      fig.path <- plot.path
    } else {
      fig.path <- getwd()
    }

    if (class.type == 'scpi_data'){

      scplot(result = to.return, fig.path = fig.path,
             fig.name = fig.name, fig.format = "png", save.data = save.data)

    } else if (class.type == "scpi_data_multi"){

      scplotMulti(result = to.return)

    }

  }

  return(to.return)
}

