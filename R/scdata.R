#' @title Create Synthetic Control Data Structure
#' @description This function prepares a data structure for use with Synthetic Control Methods (SCM) 
#' by formatting and organizing the input data. It validates inputs, structures the data appropriately, 
#' and creates matrices needed for synthetic control estimation.
#'
#' @param df Data frame containing the panel data in long format.
#' @param id.var Character. Name of the column containing unit identifiers (e.g., state names, country codes).
#' @param time.var Character. Name of the column containing time periods (e.g., year, quarter).
#' @param outcome.var Character. Name of the column containing the outcome variable of interest.
#' @param period.pre Numeric vector. Pre-treatment time periods used for fitting the synthetic control.
#' @param period.post Numeric vector. Post-treatment time periods for evaluating treatment effects.
#' @param unit.tr Character. Name or identifier of the treated unit.
#' @param unit.co Character vector. Names or identifiers of the potential control units (donor pool).
#' @param anticipation Numeric. Number of pre-treatment periods to exclude due to anticipation effects. Default is 0.
#' @param constant Logical. Whether to include a global constant (intercept) term in synthetic control estimation. 
#'   When TRUE, adds a constant feature that allows the synthetic control to have a level shift relative to donors.
#'   Requires multiple features (more than one covariate specification). The constant term is estimated alongside 
#'   donor weights and can be positive or negative. Default is FALSE.
#' @param verbose Logical. Whether to print diagnostic messages and warnings. Default is TRUE.
#' @param covagg List of covariate-aggregation operations.
#'   Each entry must be a list with:
#'   \itemize{
#'     \item \code{var}: variable name, or \code{"outcome_var"} for \code{outcome.var}
#'     \item \code{select_periods}: optional period selector
#'     \item \code{partition_periods}: optional grouping of selected periods
#'     \item \code{compute}: optional aggregation function/function name (default \code{"mean"})
#'   }
#'   If \code{covagg} is empty, \code{scdata()} defaults to outcome per-period matching
#'   over all pre-treatment periods.
#' @param cointegrated.data Logical. Whether to treat the data as cointegrated (affects some computations). Default is FALSE.
#'
#' @details
#' This function performs several key operations:
#' \itemize{
#'   \item Validates that all required variables exist in the data
#'   \item Checks for balanced panel structure
#'   \item Creates outcome matrices for treated and control units
#'   \item Handles missing data appropriately
#'   \item Prepares data structures required by estimation functions
#' }
#'
#' The function expects panel data in long format where each row represents one unit-time observation.
#'
#' @return A list of class 'scdata' containing:
#' \itemize{
#'   \item \code{Y.pre}: Matrix of pre-treatment outcomes (controls × time periods)
#'   \item \code{Y.post}: Matrix of post-treatment outcomes
#'   \item \code{Y.donors}: Full outcome matrix for donor units
#'   \item \code{Y.treated}: Outcome vector for treated unit
#'   \item \code{specs}: List containing specification details
#'   \item Additional metadata for estimation and plotting
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with state-level data
#' scm_data <- scdata(df = state_panel, 
#'                    id.var = "state", 
#'                    time.var = "year", 
#'                    outcome.var = "gdp_per_capita",
#'                    period.pre = 1980:1999, 
#'                    period.post = 2000:2010,
#'                    unit.tr = "California", 
#'                    unit.co = c("New York", "Texas", "Florida"))
#'                    
#' # With anticipation effects
#' scm_data <- scdata(df = policy_data,
#'                    id.var = "country",
#'                    time.var = "year", 
#'                    outcome.var = "unemployment_rate",
#'                    period.pre = 1990:2004,
#'                    period.post = 2007:2015,
#'                    unit.tr = "Germany",
#'                    unit.co = c("France", "Italy", "Spain"),
#'                    anticipation = 2)  # Exclude 2 pre-treatment years
#'                    
#' # Outcome per-period + covariate means
#' scm_data_ops <- scdata(df = germany_data,
#'                        id.var = "country",
#'                        time.var = "year",
#'                        outcome.var = "gdp",
#'                        period.pre = 1960:1990,
#'                        period.post = 1991:2003,
#'                        unit.tr = "West Germany",
#'                        unit.co = c("USA", "UK", "France"),
#'                        covagg = list(
#'                          list(
#'                            var = "outcome_var",
#'                            partition_periods = list(type = "by_period")
#'                          ),
#'                          list(var = "investment", compute = "mean"),
#'                          list(var = "trade", compute = "mean")
#'                        ))
#'                          
#' # With constant term (requires multiple features)
#' scm_data_const <- scdata(df = germany_data,
#'                          id.var = "country",
#'                          time.var = "year", 
#'                          outcome.var = "gdp",
#'                          period.pre = 1960:1990,
#'                          period.post = 1991:2003,
#'                          unit.tr = "West Germany",
#'                          unit.co = c("USA", "UK", "France"),
#'                          constant = TRUE,  # Enable constant term
#'                          covagg = list(
#'                            list(var = "outcome_var", partition_periods = list(type = "by_period")),
#'                            list(var = "investment", compute = "mean"),
#'                            list(
#'                              var = "trade",
#'                              select_periods = list(type = "explicit", periods = c(1985, 1990))
#'                            )
#'                          ))
#' }


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
  
  # Convert from tibble or data.table to proper dataframe object
  if (tibble::is_tibble(data) == TRUE) {
    data <- as.data.frame(data)
  }
  # Convert from data.table to data.frame to avoid column access issues
  if (data.table::is.data.table(data)) {
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
  
  # Resolve covagg period selection for scdata
  resolve_covagg_periods <- function(cov_spec, period.pre) {
    if (is.null(cov_spec$periods)) {
      stop("Internal covagg specification is missing periods")
    }
    periods <- sort(unique(cov_spec$periods))
    if (!is.numeric(periods) || length(periods) == 0) {
      stop("Internal covagg specification has invalid periods")
    }
    missing_periods <- setdiff(periods, period.pre)
    if (length(missing_periods) > 0) {
      stop("Internal covagg specification contains periods outside period.pre: ",
           paste(missing_periods, collapse = ", "))
    }
    periods
  }
  
  # Resolve aggregation function for scdata
  resolve_covagg_aggfun <- function(cov_spec) {
    aggfun <- NULL
    if (!is.null(cov_spec$agg_fun)) {
      aggfun <- cov_spec$agg_fun
    } else if (!is.null(cov_spec$aggfunc)) {
      aggfun <- cov_spec$aggfunc
    }
    
    if (is.null(aggfun)) {
      return(mean)
    }
    
    if (is.character(aggfun)) {
      aggfun <- match.fun(aggfun)
    }
    
    if (!is.function(aggfun)) {
      stop("covagg agg_fun/aggfunc must be a function or a function name")
    }
    
    return(aggfun)
  }
  
  # Apply aggregation function with optional na.rm
  apply_aggfun <- function(aggfun, values) {
    aggfun_formals <- formals(aggfun)
    if (!is.null(aggfun_formals) && "na.rm" %in% names(aggfun_formals)) {
      return(aggfun(values, na.rm = TRUE))
    }
    return(aggfun(values))
  }

  if (length(covagg) == 0) {
    covagg <- list(
      list(
        var = "outcome_var",
        select_periods = list(type = "all_pre"),
        partition_periods = list(type = "by_period"),
        compute = "mean"
      )
    )
  }

  covagg <- compile_covagg_signature(
    covagg = covagg,
    outcome_var = outcome.var,
    period_pre = period.pre,
    require_named = FALSE,
    context = "scdata covagg"
  )
  
  ############################################################################
  ############################################################################
  ############################################################################
  
  # Handle covariate aggregation with nested specifications
  if (length(covagg) > 0) {
    # Extract features from specifications
    features = c()
    for (spec_name in names(covagg)) {
      spec <- covagg[[spec_name]]
      # Skip the label
      if (spec_name == "label") {
        next
      }
      
      # The spec itself contains the var element
      if (!is.null(spec$var)) {
        features = c(features, spec$var)
      }
    }
    features = unique(features)
  } else {
    features <- c()
  }
  
  if (length(features) == 0) {
    features <- outcome.var
  }

  if (outcome.var %in% features) {
    out.in.features <- TRUE
  } else {
    out.in.features <- FALSE
  }


  period.pre  <- sort(period.pre, decreasing = FALSE)
  period.post <- sort(period.post, decreasing = FALSE)
  # replace space and period with underscore to avoid issues below
  unit.tr = gsub(' +|\\.+|-+', '_', unit.tr)
  unit.co = gsub(' +|\\.+|-+', '_', unit.co)

  if(is.character(data[[id.var]])){
    data[[id.var]] = gsub(' +|\\.+|-+', '_', data[[id.var]] )
  }
  # Create ID and time variables
  if (is.numeric(data[[id.var]])) {
    data[[id.var]] <- as.character(data[[id.var]])
  }

  # Check that specified units are in dataframe
  if (!(unit.tr %in% data[[id.var]])) {
    stop("There is no treated unit with the specified ID (unit.tr) in the specified ID variable (id.var)!")
  }

  if (!all(unit.co %in% data[[id.var]])) {
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
  if (!all(period.pre %in% data[[time.var]])) {
    pre.not.found <- period.pre[!(period.pre %in% data[[time.var]])]
    stop(paste(c("The following pre-treatment period(s) are not in the input dataframe:",
                 pre.not.found), collapse = " "))
  }
  
  if (!all(period.post %in% data[[time.var]])) {
    post.not.found <- period.post[!(period.post %in% data[[time.var]])]
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

  # Clean column names for data.table compatibility
  clean_names <- gsub('\\.|-+| +' , '_', names(data))
  names(data) <- clean_names
  
  # Create ID and time variables after renaming
  id          <- as.matrix(unique(data[["ID"]]))        # ID of units
  time        <- unique(data[["Time"]])      # Time periods
  time        <- time[time %in% c(period.pre, period.post)]
  
  data.bal <- as.data.frame(tidyr::complete(data, ID, Time))
  
  
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
    
  # Process covariate specifications (nested format)
  if (length(covagg) > 0) {
    
    for (cov_name in names(covagg)) {
      cov_spec <- covagg[[cov_name]]
      
      # Process each feature within this specification
      # for (cov_name in names(spec)) {
      # Skip the label
      
      if (cov_name == "label") {
        next
      }
      
      # Skip non-list elements or missing var
      if (!is.list(cov_spec) || is.null(cov_spec$var)) {
        next
      }
      
      var_name <- cov_spec$var
      var_name = gsub('\\.+| +|-+', '_', var_name)

      
      # Check if variable exists in data
      available_cols <- colnames(data.bal.tr)
      
      if (!(var_name %in% available_cols)) {
        stop("Covariate variable '", var_name, "' not found in input data while constructing treated-unit features.")
      }
      
      # Determine periods to use
      periods_to_use <- resolve_covagg_periods(cov_spec, period.pre)
      aggfun <- resolve_covagg_aggfun(cov_spec)

      # Filter data for the specified periods
      period_data <- data.bal.tr[data.bal.tr[, "Time"] %in% periods_to_use, c(outcome.var, var_name)]

      if (nrow(period_data) > 0) {
        idx2 = length(A_list) + 1

        # Apply aggregation function
        agg_value <- apply_aggfun(aggfun, period_data[[var_name]])

        A <- as.matrix(agg_value)
        colnames(A) <- unit.tr
        rownames(A) <- paste(unit.tr, cov_name, sep='.')
        A_list[[idx2]] = data.frame(A)
      }
      # }
    }
  }
  
  # Create B
  sel <- data.bal[rows.co.pre, ] # Select rows and columns
  B_list = list()
  
  data.bal.co = data.bal[rows.co.pre, ]
  
  feature.vec = c()
  time.vec = c()
  
  # Process covariate specifications (nested format)
  
  # Process covariate specifications for B matrix (nested format)
  if (length(covagg) > 0) {
    for (cov_name in names(covagg)) {
      cov_spec <- covagg[[cov_name]]
      # Skip the label
      if(cov_name=='label') {
        next
      }
      
      
      # Skip non-list elements or missing var
      if (!is.list(cov_spec) || is.null(cov_spec$var)) {
        next
      }
      
      
      var_name <- cov_spec$var
      var_name = gsub('\\.', '_', var_name)

      
      # Check if variable exists in data for B matrix
      available_cols_b <- colnames(data.bal.co)
      if (!(var_name %in% available_cols_b)) {
        stop("Covariate variable '", var_name, "' not found in input data while constructing donor features.")
      }
      # Determine periods to use
      periods_to_use <- resolve_covagg_periods(cov_spec, period.pre)
      aggfun <- resolve_covagg_aggfun(cov_spec)
      temp = data.bal.co[data.bal.co[, "Time"] %in% periods_to_use, c('ID', 'Time', var_name)]

      if (nrow(temp) > 0) {
        # Apply aggregation by ID
        agg_data <- aggregate(
          temp[[var_name]],
          by = list(ID = temp$ID),
          FUN = function(x) apply_aggfun(aggfun, x)
        )

        # Convert to matrix format
        B <- data.frame(t(agg_data$x))

        # IMPORTANT: Use the exact same transformation that Y.donors will use
        # This matches the str_remove + gsub pattern in Y.donors creation
        temp_names <- paste("gdpcap", agg_data$ID, sep = ".")
        y_names_equivalent <- stringr::str_remove(temp_names, "gdpcap")
        y_names_equivalent <- stringr::str_remove(y_names_equivalent, "\\.")
        final_clean_names <- gsub('[^A-Za-z0-9_]', '_', y_names_equivalent)
        final_clean_names <- gsub('_{2,}', '_', final_clean_names)
        # Remove trailing underscores to match Y.donors
        final_clean_names <- gsub('_+$', '', final_clean_names)

        clean_unit_tr <- gsub('[^A-Za-z0-9_]', '_', unit.tr)
        clean_unit_tr <- gsub('_{2,}', '_', clean_unit_tr)

        # B matrix columns should match Y.donors exactly
        # Only set column names if dimensions match
        if (length(final_clean_names) == ncol(B)) {
          colnames(B) <- final_clean_names
        } else {
          warning(paste("Column name dimension mismatch in aggregation: B has", ncol(B), "columns but", length(final_clean_names), "names"))
        }
        rownames(B) <- paste(clean_unit_tr, cov_name, sep='.')

        idx2 = length(B_list) + 1
        B_list[[idx2]] = B
      }
      #}
    }
  }
  
  feature.vec = feature.vec[feature.vec!='Time']
  # feature.vec = feature.vec[!duplicated(feature.vec)]
  
  A = dplyr::bind_rows(A_list)
  
  B = dplyr::bind_rows(B_list)
  
  
  A = as.matrix(A)
  
  B = as.matrix(B)
  
  
  ## Create matrix with pre-period outcomes for the donors
  sel <- data.bal[rows.co.pre, c(outcome.var, "ID", "Time")] # Select rows and columns
  aux <- stats::reshape(sel,
                        direction = "wide",
                        idvar     = "Time",
                        timevar   = "ID")
  
  colnames(aux)     <- stringr::str_remove(colnames(aux), outcome.var)
  colnames(aux)     <- stringr::str_remove(colnames(aux),"\\.")
  
  Y.donors <- as.matrix(aux[, colnames(aux) != "Time", drop=FALSE])
  
  Y.pre <- as.matrix(data.bal[rows.tr.pre, outcome.var])
  rownames(Y.pre) <- paste(unit.tr,as.character(data.bal[rows.tr.pre, "Time"]), sep=".")
  colnames(Y.pre) <- outcome.var
  
  # if "Time" is not in numeric format the matrix becomes a matrix of characters,
  # this is why we do everything in one step
  
  # # Use underscores consistently - clean Y.names to match B matrix format
  # clean_y_names <- gsub('[^A-Za-z0-9_]', '_', Y.names)
  # clean_y_names <- gsub('_{2,}', '_', clean_y_names)
  # # Remove trailing underscores
  # clean_y_names <- gsub('_+$', '', clean_y_names)
  
  # # Y.donors columns should match B matrix columns (no unit.tr prefix)
  # colnames(Y.donors) <- clean_y_names
  rownames(Y.donors) <- paste(unit.tr, as.character(aux[,'Time']), sep = ".")
  
  
  Y.donors    <- Y.donors[ , colnames(B)]  # Re-order according to B
  
  ## Create matrix with post-period outcomes for the donors
  sel <- data.bal[rows.co.post, c(outcome.var, "ID", "Time")] # Select rows and columns
  
  aux.post <- stats::reshape(sel,
                             direction = "wide",
                             idvar     = "Time",
                             timevar   = "ID")
  colnames(aux.post)     <- stringr::str_remove(colnames(aux.post), outcome.var)
  colnames(aux.post)     <- stringr::str_remove(colnames(aux.post),"\\.")
  
  # if "Time" is not in numeric format the matrix becomes a matrix of characters,
  # this is why we do everything in one step
  Y.donors.post <- as.matrix(aux.post[, colnames(aux) != "Time", drop=FALSE])
  
  # Y.donors.post columns should match B matrix columns (no unit.tr prefix)
  # colnames(Y.donors.post) <- clean_y_names_post
  rownames(Y.donors.post) <- paste(unit.tr, as.character(aux.post[,'Time']), sep = ".")
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
  
  colnames(aux)     <- stringr::str_remove(colnames(aux), outcome.var)
  colnames(aux)     <- stringr::str_remove(colnames(aux),"\\.")


  P <- as.matrix(aux[, names(aux) != "Time"])
  
  rownames(P) <- paste(unit.tr, as.character(aux[,'Time']), sep = ".")
  
  P <- P[, colnames(B), drop = F]  # Re-order as the matrix B
  
  # If the outcome variable is within the specified features then we need to
  # augment P with the corresponding (eventual) covariates used for adjustment.
  # If instead the outcome variable is not within the specified features
  # P is composed by the outcome variable of the donors only
  
  if (out.in.features == TRUE) {
    # No additional modifications needed - constant terms handled via C matrix
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
  
  A_outcome = A[grepl(outcome.var, rownames(A)), , drop=F]
  A_nonoutcome = A[!(grepl(outcome.var, rownames(A))), ,drop=F]
  # Only order if A_outcome has rows and rownames
  if (nrow(A_outcome) > 0 && !is.null(rownames(A_outcome))) {
    A_outcome = A_outcome[order(rownames(A_outcome)), , drop=F]
  }
  if(nrow(A_nonoutcome)> 0){
    # Only order if A_nonoutcome has rownames
    if (!is.null(rownames(A_nonoutcome))) {
      A_nonoutcome = A_nonoutcome[order(rownames(A_nonoutcome)),, drop=F ]
    }
    A = rbind(A_outcome, A_nonoutcome)
  } else{
    A = A_outcome
  }
  
  # Matrix A kept unchanged - constant terms now handled in separate C matrix
  
  B_outcome = B[grepl(outcome.var, rownames(B)), , drop=F]
  B_nonoutcome = B[!(grepl(outcome.var, rownames(B))), ,drop=F]
  # Only order if B_outcome has rows and rownames
  if (nrow(B_outcome) > 0 && !is.null(rownames(B_outcome))) {
    B_outcome = B_outcome[order(rownames(B_outcome)), , drop=F]
  }
  if(nrow(B_nonoutcome)>0){
    # Only order if B_nonoutcome has rownames
    if (!is.null(rownames(B_nonoutcome))) {
      B_nonoutcome = B_nonoutcome[order(rownames(B_nonoutcome)),, drop=F ]
    }
    B = rbind(B_outcome, B_nonoutcome)
  } else{
    B = B_outcome
  }
  
  
  
  # Matrix B kept unchanged - constant terms now handled in separate C matrix

  # Create C matrix for covariate adjustments (following SCPI approach)
  C <- NULL
  if (constant == TRUE) {
    C <- matrix(1, nrow = nrow(B), ncol = 1)
    colnames(C) <- paste(unit.tr, "0.constant", sep = ".")
    rownames(C) <- rownames(B)
  }

  # Combine matrices following SCPI structure: X = [A, B, C]
  if (!is.null(C)) {
    X <- cbind(A, B, C)
  } else {
    X <- cbind(A, B)
  }
  
  # Create rownames that match the actual dimensions of X
  if (nrow(X) > 0) {
    expected_rownames <- paste(unit.tr, feature.vec, as.character(time.vec), sep = ".")
    if (length(expected_rownames) == nrow(X)) {
      rownames(X) <- expected_rownames
    } else {
      # If dimensions don't match, create simple rownames matching actual rows
      # Dimension mismatch - using simple rownames
      rownames(X) <- paste0("row_", 1:nrow(X))
    }
  }
  
  
  select      <- rowSums(is.na(X)) == 0
  
  X.na        <- X[select, , drop = FALSE]
  if (nrow(X.na) == 0) {
    stop("Current specification has too many missing values and no observations are left!")
  }
  
  j1   <- dim(as.matrix(A))[2]  # Columns of A
  j2   <- j1 + dim(B)[2]        # Columns of B
  j3   <- j2 + if(!is.null(C)) dim(C)[2] else 0  # Columns of C
  
  A.na           <- X.na[, 1:j1, drop = FALSE]
  B.na           <- X.na[, (j1+1):j2, drop = FALSE]
  C.na           <- if(!is.null(C)) X.na[, (j2+1):j3, drop = FALSE] else NULL
  
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
    # Safely set P column names
    if (length(colnames(B.na)) == ncol(P)) {
      colnames(P) <- colnames(B.na)
    } else {
      warning(paste("Final P matrix column name mismatch: P has", ncol(P), "columns but B.na has", length(colnames(B.na)), "column names"))
    }
  } else {
    # Handle column names for P matrix, accounting for potential constant term
    if (constant == TRUE && ncol(P) == length(colnames(B.na)) + 1) {
      # P matrix has extra column for constant term
      colnames(P) <- c(colnames(B.na), paste(unit.tr, "constant", sep = "."))
    } else if (length(colnames(B.na)) == ncol(P)) {
      colnames(P) <- c(colnames(B.na))
    } else {
      warning(paste("Final P matrix column name mismatch (else): P has", ncol(P), "columns but B.na has", length(colnames(B.na)), "column names"))
    }
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
                units.est = unit.tr,
                col.name.unit = id.var,
                col.name.time = time.var)
  
  df.sc <- list(A = A.na,
                B = B.na,
                C = C.na,
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
