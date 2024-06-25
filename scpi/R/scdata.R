scdata <- function(df,
                   id.var,
                   time.var,
                   outcome.var,
                   period.pre,
                   period.post,
                   unit.tr,
                   unit.co,
                   anticipation = 0,
                   constant = FALSE,
                   verbose = TRUE,
                   covagg = list(),
                   cointegrated.data = FALSE) {

  ############################################################################
  ############################################################################
  ### Error Checking

  #  Safe copy
  data <- df

  # Store variable names and variable class
  var.names   <- names(data)    # Var names in dataframe
  var.class   <- sapply(data, class)      # Var types in dataframe

  # Check main input is a dataframe
  if (is.data.frame(data) == FALSE) {
    stop("Data input should be a dataframe object!")
  }

  # Convert from tibble to proper dataframe object
  if (tibble::is_tibble(data) == TRUE) {
    data <- as.data.frame(data)
  }

  # Check inputs are string
  if (is.character(id.var) == FALSE) {
    stop("You should specify the name of id.var as a character! (eg. id.var = 'ID')")
  }

  if (is.character(outcome.var) == FALSE) {
    stop("You should specify the name of outcome.var as a character! (eg. outcome.var = 'outcome')")
  }

  if (is.character(time.var) == FALSE) {
    stop("You should specify the name of time.var as a character! (eg. time.var = 'time')")
  }


  # Check inputs are in dataframe
  if (!(id.var %in% var.names)) {
    stop("ID variable (id.var) not found in the input dataframe!")
  }

  if (!(time.var %in% var.names)) {
    stop("Time variable (time.var) not found in the input dataframe!")
  }

  if (!(outcome.var %in% var.names)) {
    stop("Outcome variable (outcome.var) not found in the input dataframe!")
  }

  time.var.class <- var.class[var.names == time.var]
  if (!time.var.class %in% c("numeric", "integer", "Date")) {
    stop("Time variable (time.var) must be either numeric or character!")
  }

  if (!var.class[var.names == outcome.var] %in% c("numeric", "integer")) {
    stop("Outcome variable (outcome.var) must be numeric!")
  }

  features = c()
  for(item in names(covagg)){
    for(feat in covagg[[item]]){
      features = c(features, feat)
    }
  }
  features = unique(features)



  if (is.null(features) ==  TRUE) {
    features <- outcome.var
  }

  if (outcome.var %in% features) {
    out.in.features <- TRUE
  } else {
    out.in.features <- FALSE
  }


  if ((length(features) == 1) && (constant == TRUE)) {
    stop("When specifying just one feature you either specify constant == TRUE")
  }

  period.pre  <- sort(period.pre, decreasing = FALSE)
  period.post <- sort(period.post, decreasing = FALSE)

  # replace space and period with underscore to avoid issues below
  unit.tr = gsub(' |\\.', '_', unit.tr)
  unit.co = gsub(' |\\.', '_', unit.co)

  if(is.character(data[, id.var])){
    data[, id.var] = gsub(' +|\\.+', '_', data[, id.var] )
  }


  # Create ID and time variables
  if (is.numeric(data[, id.var])) {
    data[id.var] <- as.character(data[id.var])
  }

  id          <- as.matrix(unique(data[id.var]))        # ID of units
  time        <- unique(data[,time.var])      # Time periods
  time        <- time[time %in% c(period.pre, period.post)]

  # Check that specified units are in dataframe
  if (!(unit.tr %in% id)) {
    stop("There is no treated unit with the specified ID (unit.tr) in the specified ID variable (id.var)!")
  }

  if (!all(unit.co %in% id)) {
    co.not.found <- unit.co[!(unit.co %in% id)]

    stop(paste(c("The following control unit(s) are not in the input dataframe:",
                 co.not.found), collapse = " "))
  }

  if (unit.tr %in% unit.co) {
    stop("The treated unit is also contained in the donor pool!")
  }

  if (length(unit.co) < 2) {
    stop("Please provide at least two control units!")
  }

  # Check specified time periods are in dataframe
  if (!all(period.pre %in% time)) {
    pre.not.found <- period.pre[!(period.pre %in% time)]
    stop(paste(c("The following pre-treatment period(s) are not in the input dataframe:",
                 pre.not.found), collapse = " "))
  }

  if (!all(period.post %in% time)) {
    post.not.found <- period.post[!(period.post %in% time)]
    stop(paste(c("The following post-treatment period(s) are not in the input dataframe:",
                 post.not.found), collapse = " "))
  }

  if (any(period.pre %in% period.post)) {
    stop("There is an overlap between the pre-treatment period and post-treatment period!")
  }

  # Consider eventual anticipation effect
  if (!is.numeric(anticipation)) {
    stop("The object 'anticipation' has to be an integer!")
  }

  if (anticipation > 0) {
    t0 <- length(period.pre); d <- anticipation
    period.post <- c(period.pre[(t0-d+1):t0], period.post)
    period.pre  <- period.pre[1:(t0-d)]
  }

  # Order outcome variable as first feature and handle covariates for adjustment accordingly if needed
  if (out.in.features == TRUE) {
    features <- c(features[features == outcome.var], features[features != outcome.var])
  }

  # Rename time and ID variables
  var.names[var.names == time.var] <- "Time"
  var.names[var.names == id.var]   <- "ID"
  names(data) <- var.names


  ############################################################################
  ############################################################################
  ### Data preparation

  # Make the panel balanced
  data.bal <- as.data.frame(tidyr::complete(data, .data[["ID"]],.data[["Time"]]))


  # Identify rows corresponding to treatment unit
  rows.tr.pre <- which(data.bal[, "ID"] %in% c(unit.tr) &
                         data.bal[, "Time"] %in% period.pre)

  # Identify rows corresponding to control units
  rows.co.pre <- which(data.bal[, "ID"] %in% c(unit.co) &
                         data.bal[, "Time"] %in% period.pre)

  # Identify rows corresponding to treatment unit
  rows.tr.post <- which(data.bal[, "ID"] %in% c(unit.tr) &
                          data.bal[, "Time"] %in% period.post)

  # Identify rows corresponding to control units
  rows.co.post <- which(data.bal[, "ID"] %in% c(unit.co) &
                          data.bal[, "Time"] %in% period.post)

  min_period = min(period.pre)
  max_period = max(period.pre)

  A_list = list()

  data.bal.tr = data.bal[rows.tr.pre, ]
  for(cov_agg in names(covagg)){
    if(cov_agg=='every_other_period'){
      feature_names = covagg[[cov_agg]]
      feature_names = c(outcome.var, feature_names[feature_names!=outcome.var])
      n_features = length(feature_names)

      years = (period.pre[period.pre%%2==0])
      idx2 = length(A_list) + 1
      A <- as.matrix(c(as.matrix(data.bal.tr[data.bal.tr[, "Time"] %in% years,
                                             feature_names]))) # Stack features
      colnames(A) <- unit.tr
      rownames(A) <- paste(rep(unit.tr, n_features * length(years)),
                           rep(feature_names, each=length(years)),
                           rep(years, n_features) , sep='.')
      A_list[[idx2]] = data.frame(A)
    }

    if(cov_agg=='every_period'){
      feature_names = covagg[[cov_agg]]
      feature_names = c(outcome.var, feature_names[feature_names!=outcome.var])
      n_features = length(feature_names)

      years = (period.pre)
      idx2 = length(A_list) + 1
      A <- as.matrix(c(as.matrix(data.bal.tr[data.bal.tr[, "Time"] %in% years,
                                             feature_names]))) # Stack features
      colnames(A) <- unit.tr
      rownames(A) <- paste(rep(unit.tr, n_features * length(years)),
                           rep(feature_names, each=length(years)),
                           rep(years, n_features) , sep='.')
      A_list[[idx2]] = data.frame(A)
    }

    if(cov_agg=='mean'){
      feature_names = covagg[[cov_agg]]
      feature_names = c(outcome.var, feature_names[feature_names!=outcome.var])
      n_features = length(feature_names)

      years = (period.pre)
      idx2 = length(A_list) + 1

      A <- as.matrix(c(as.matrix(colMeans(data.bal.tr[data.bal.tr[, "Time"] %in% years,
                                                      feature_names])))) # Stack features
      colnames(A) <- unit.tr
      rownames(A) <- paste(rep(unit.tr, n_features),
                           feature_names,
                           rep(paste0(min_period, '.', max_period) , n_features) , sep='.')
      A_list[[idx2]] = data.frame(A)
    }


    if(cov_agg=='last_period'){
      feature_names = covagg[[cov_agg]]
      feature_names = c(outcome.var, feature_names[feature_names!=outcome.var])
      n_features = length(feature_names)

      years = max_period

      idx2 = length(A_list) + 1
      A <- as.matrix(c(as.matrix(data.bal.tr[data.bal.tr[, "Time"] %in% years,
                                             feature_names]))) # Stack features
      colnames(A) <- unit.tr
      rownames(A) <- paste(rep(unit.tr, n_features),
                           rep(feature_names, each=1),
                           rep(max_period, n_features) , sep='.')
      A_list[[idx2]] = data.frame(A)
    }


  }

  ### Estimation Data
  # Actual Pre-treatment Series
  Y.pre <- as.matrix(data.bal[rows.tr.pre, outcome.var])
  rownames(Y.pre) <- paste(unit.tr,as.character(data.bal[rows.tr.pre, "Time"]),sep=".")
  colnames(Y.pre) <- outcome.var

  # Create B
  sel <- data.bal[rows.co.pre, ] # Select rows and columns

  B_list = list()

  data.bal.co = data.bal[rows.co.pre, ]

  feature.vec = c()
  time.vec = c()

  for(cov_agg in names(covagg)){
    if(cov_agg=='every_other_period'){
      feature_names = covagg[[cov_agg]]
      feature_names = c(outcome.var, feature_names[feature_names!=outcome.var])
      n_features = length(feature_names)

      years = (period.pre[period.pre%%2==0])
      temp = data.bal.co[data.bal.co[, "Time"] %in% years,
                         c('ID', 'Time', feature_names)]
      temp = data.table::melt(data.table(temp), id.vars=c('ID', 'Time'))

      B = data.frame(data.table::dcast(temp, variable+Time~ID))

      rownames(B) = paste(rep(unit.tr, nrow(B)), B$variable, B$Time, sep='.')

      idxes = which(colnames(B) %in% c('variable', 'Time'))

      feature.vec = c(feature.vec, as.character(B$variable))

      time.vec = c(time.vec, as.numeric(B$Time))

      B = B[, -idxes]
      colnames(B) = paste(unit.tr, colnames(B), sep='.')
      idx2 = length(B_list) + 1
      B_list[[idx2]] = B
    }

    if(cov_agg=='every_period'){
      feature_names = covagg[[cov_agg]]
      feature_names = c(outcome.var, feature_names[feature_names!=outcome.var])
      n_features = length(feature_names)

      years = period.pre
      temp = data.bal.co[data.bal.co[, "Time"] %in% years,
                         c('ID', 'Time', feature_names)]

      temp = data.table::melt(data.table(temp), id.vars=c('ID', 'Time'))
      B = data.frame(data.table::dcast(temp, variable+Time~ID))



      rownames(B) = paste(rep(unit.tr, nrow(B)),
                          B$variable,
                          B$Time, sep='.')

      idxes = which(colnames(B) %in% c('variable', 'Time'))
      feature.vec = c(feature.vec, as.character(B$variable))
      time.vec = c(time.vec, as.numeric(B$Time))

      B = B[, -idxes]
      colnames(B) = paste(unit.tr, colnames(B), sep='.')
      idx2 = length(B_list) + 1
      B_list[[idx2]] = B
    }


    if(cov_agg=='mean'){
      feature_names = covagg[[cov_agg]]
      feature_names = c(outcome.var, feature_names[feature_names!=outcome.var])
      n_features = length(feature_names)

      years = period.pre
      temp = data.bal.co[data.bal.co[, "Time"] %in% years,
                         c('ID', 'Time', feature_names)]

      temp = data.table::melt(data.table(temp), id.vars=c('ID', 'Time'))
      temp = temp[, .(value=mean(value)), by=c('ID', 'variable')]
      B = data.frame(data.table::dcast(temp, variable~ID))


      rownames(B) = paste(rep(unit.tr, nrow(B)),
                          B$variable,
                          rep(paste0(min_period, '.', max_period), nrow(B)), sep='.')

      idxes = which(colnames(B) %in% c('variable', 'Time'))
      feature.vec = c(feature.vec, as.character(B$variable))
      time.vec = c(time.vec, as.numeric(B$Time))

      B = B[, -idxes]
      colnames(B) = paste(unit.tr, colnames(B), sep='.')
      idx2 = length(B_list) + 1
      B_list[[idx2]] = B


    }


    if(cov_agg=='last_period'){
      feature_names = covagg[[cov_agg]]
      feature_names = c(outcome.var, feature_names[feature_names!=outcome.var])
      n_features = length(feature_names)

      years = max_period
      temp = data.bal.co[data.bal.co[, "Time"] %in% years,
                         c('ID', 'Time', feature_names)]

      temp = data.table::melt(data.table(temp), id.vars=c('ID', 'Time'))
      B = data.frame(data.table::dcast(temp, variable+Time~ID))



      rownames(B) = paste(rep(unit.tr, nrow(B)),
                          B$variable,
                          B$Time, sep='.')

      idxes = which(colnames(B) %in% c('variable', 'Time'))
      feature.vec = c(feature.vec, as.character(B$variable))
      time.vec = c(time.vec, as.numeric(B$Time))

      B = B[, -idxes]
      colnames(B) = paste(unit.tr, colnames(B), sep='.')
      idx2 = length(B_list) + 1
      B_list[[idx2]] = B
    }
  }


  A = dplyr::bind_rows(A_list)
  rnA = gsub('\\.\\.\\.\\d+$', '', rownames(A))
  if(sum(duplicated(rnA))>0){
    A = A[-which(duplicated(rnA)), , drop=F]
  }
  rownames(A) = rnA[!duplicated(rnA)]
  A = as.matrix(A)

  B = dplyr::bind_rows(B_list)
  rnB = gsub('\\.\\.\\.\\d+$', '', rownames(B))

  if(sum(duplicated(rnB))>0){
    B = B[-which(duplicated(rnB)),, drop=F]
  }

  rownames(B) = rnB[!duplicated(rnB)]

  colnames(B) = gsub(' ', '_', colnames(B))

  B = as.matrix(B)
  feature.vec = feature.vec[feature.vec!='Time']
  # feature.vec = feature.vec[!duplicated(feature.vec)]


  ## Create matrix with pre-period outcomes for the donors
  sel <- data.bal[rows.co.pre, c(outcome.var, "ID", "Time")] # Select rows and columns
  aux <- stats::reshape(sel,
                        direction = "wide",
                        idvar     = "Time",
                        timevar   = "ID")

  # if "Time" is not in numeric format the matrix becomes a matrix of characters,
  # this is why we do everything in one step
  Y.donors <- as.matrix(aux[, colnames(aux) != "Time", drop=FALSE])
  Y.names     <- stringr::str_remove(colnames(Y.donors), outcome.var)
  Y.names     <- stringr::str_remove(Y.names,".")
  colnames(Y.donors) <- paste(unit.tr, Y.names, sep = ".")
  rownames(Y.donors) <- paste(unit.tr, as.character(aux[,'Time']), sep = ".")
  colnames(Y.donors) = gsub(' ', '_', colnames(Y.donors))
  # print(colnames(B))
  # print(colnames(Y.donors))
  Y.donors    <- Y.donors[ , colnames(B)]  # Re-order according to B

  ## Create matrix with post-period outcomes for the donors
  sel <- data.bal[rows.co.post, c(outcome.var, "ID", "Time")] # Select rows and columns
  aux.post <- stats::reshape(sel,
                             direction = "wide",
                             idvar     = "Time",
                             timevar   = "ID")

  # if "Time" is not in numeric format the matrix becomes a matrix of characters,
  # this is why we do everything in one step
  Y.donors.post <- as.matrix(aux.post[, colnames(aux) != "Time", drop=FALSE])

  Y.names.post     <- stringr::str_remove(colnames(Y.donors.post), outcome.var)
  Y.names.post    <- stringr::str_remove(Y.names.post,".")


  colnames(Y.donors.post) <- paste(unit.tr, Y.names.post, sep = ".")
  rownames(Y.donors.post) <- paste(unit.tr, as.character(aux.post[,'Time']), sep = ".")
  colnames(Y.donors.post) = gsub(' ', '_', colnames(Y.donors.post))

  Y.donors.post    <- Y.donors.post[ , colnames(B)]  # Re-order according to B

  ############################################################################
  ##############################################################################
  ### Prediction Data

  ## Actual post-treatment series
  Y.post <- as.matrix(data.bal[rows.tr.post, outcome.var])
  rownames(Y.post) <- paste(unit.tr, as.character(data.bal[rows.tr.post, "Time"]), sep = ".")
  colnames(Y.post) <- outcome.var

  ## Prediction Matrix
  # Select series of donors
  aux    <- data.bal[rows.co.post, c(outcome.var, "ID", "Time")] # Select rows and columns

  # make the df wide so that countries are one next to the other
  aux <- stats::reshape(aux,
                        direction = "wide",
                        idvar     = "Time",
                        timevar   = "ID")
  P <- as.matrix(aux[, names(aux) != "Time"])
  rownames(P) <- paste(unit.tr, as.character(aux[,'Time']), sep = ".")
  P.names     <- stringr::str_remove(colnames(P), outcome.var)
  colnames(P) <- paste(unit.tr, P.names, sep = ".")
  colnames(P) = gsub('\\.\\.', '.', colnames(P))
  colnames(P) = gsub(' ', '_', colnames(P))

  P <- P[, colnames(B), drop = F]  # Re-order as the matrix B

  # If the outcome variable is within the specified features then we need to
  # augment P with the corresponding (eventual) covariates used for adjustment.
  # If instead the outcome variable is not within the specified features
  # P is composed by the outcome variable of the donors only

  if (out.in.features == TRUE) {
    # Check that global constant is required by the user
    if (constant == TRUE) {
      P <- cbind(P, rep(1, length(rows.tr.post)))
      colnames(P) <- c(colnames(P[, 1:(dim(P)[2] - 1), drop = FALSE]), paste(unit.tr,"0.constant", sep = "."))
    }

    # Add covariates used for adjustment in outcome variable equation (if present)


  }

  T1 <- length(period.post)


  ############################################################################
  ############################################################################
  # Proceed cleaning missing data in the pre-treatment period
  # Check if there are annoying donors with ALL missing values in the pre-treatment period

  if(!is.null(dim(B))){
    empty.cols <- colSums(is.na(B)) == nrow(B)
  } else{
    empty.cols = 0
  }

  if (sum(empty.cols) > 0) {
    names <- strsplit(colnames(B[,empty.cols]), "\\.")
    removed.cols <- unlist(lapply(names, "[[", 2))
    warn.text <- paste(c("The following donors have no observations in the pre-treatment period, hence they have been removed!",
                         removed.cols), collapse = " ")
    if (verbose == TRUE) {
      warning(warn.text)
    }
    dropped.co <- unit.co %in% removed.cols
    unit.co.eff <- unit.co[!dropped.co]

    B <- B[ , !dropped.co, drop = TRUE]

  } else {
    unit.co.eff <- unit.co
  }

  A_outcome = A[grepl(outcome, rownames(A)), , drop=F]

  A_nonoutcome = A[!(grepl(outcome, rownames(A))), ,drop=F]

  A_outcome = A_outcome[order(rownames(A_outcome)), , drop=F]
  if(nrow(A_nonoutcome)> 0){
    A_nonoutcome = A_nonoutcome[order(rownames(A_nonoutcome)),, drop=F ]
    A = rbind(A_outcome, A_nonoutcome)
  } else{
    A = A_outcome
  }


  B_outcome = B[grepl(outcome, rownames(B)), , drop=F]
  B_nonoutcome = B[!(grepl(outcome, rownames(B))), ,drop=F]
  B_outcome = B_outcome[order(rownames(B_outcome)), , drop=F]
  if(nrow(B_nonoutcome)>0){
    B_nonoutcome = B_nonoutcome[order(rownames(B_nonoutcome)),, drop=F ]
    B = rbind(B_outcome, B_nonoutcome)
  } else{
    B = B_outcome
  }

  X <- cbind(A, B)

  rownames(X) <- paste(unit.tr, feature.vec, as.character(time.vec), sep = ".")

  select      <- rowSums(is.na(X)) == 0
  X.na        <- X[select, , drop = FALSE]
  if (nrow(X.na) == 0) {
    stop("Current specification has too many missing values and no observations are left!")
  }

  j1   <- dim(as.matrix(A))[2]  # Columns of A
  j2   <- j1 + dim(B)[2]        # Columns of B
  j3   <- j2

  A.na           <- X.na[, 1:j1, drop = FALSE]
  B.na           <- X.na[, (j1+1):j2, drop = FALSE]

  feature.na.vec <- feature.vec[select]
  time.na.vec    <- time.vec[select]

  # Store effective number of observations per feature
  xx <- as.data.frame(table(feature.na.vec))
  T0.features        <- c(xx$Freq)

  names(T0.features) <- xx$feature.na.vec

  T0.features        <- T0.features[match(names(T0.features), features)]

  ############################################################################
  ############################################################################
  # Store objects

  # Size of donor pool
  J  <- length(unit.co.eff)

  # Total number of covariates used for adjustment
  KM <- 0

  # Number of features
  M <- length(features)

  # Vector containing number of covariates used for adjustment in each equation
  K  <- rep(0, M)
  names(K) <- features

  # Augment P with zeros for other equations' cov adj
  if (sum(K[-1]) > 0 & out.in.features == TRUE) {
    zeros        <- matrix(0, nrow = nrow(P), ncol = sum(K[-1]))
    P  <- cbind(P, zeros)
  }

  if (constant == TRUE) {
    K <- K + 1
  }

  if (out.in.features == FALSE) {
    KM <- 0
    nK <- names(K)
    K <- rep(0, length(K))
    names(K) <- nK
    colnames(P) <- colnames(B.na)
  } else {
    colnames(P) <- c(colnames(B.na))
  }

  specs <- list(J = J,
                K = K,
                KM = KM,
                M = M,
                I = 1,
                cointegrated.data = cointegrated.data,
                period.pre = period.pre,
                period.post = period.post,
                T0.features = T0.features,
                T1.outcome = T1,
                outcome.var = outcome.var,
                features = features,
                constant = constant,
                out.in.features = out.in.features,
                treated.units = unit.tr,
                donors.units = unit.co.eff,
                effect = "unit-time",
                sparse.matrices = FALSE,
                units.est = unit.tr)

  df.sc <- list(A = A.na,
                B = B.na,
                P = P,
                P.diff = NULL,
                Y.pre = Y.pre,
                Y.post = Y.post,
                Y.donors = Y.donors,
                Y.donors.post = Y.donors.post,
                specs = specs)

  class(df.sc) <- 'scdata'

  return(df.sc = df.sc)


}
