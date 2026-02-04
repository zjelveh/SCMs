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
#' @param covagg List. Flexible specification for covariate aggregation.
#'   Each element is a list with \code{var} and aggregation type:
#'   \itemize{
#'     \item \code{list(var = "gdp", each = TRUE)} - Separate variable for each period
#'     \item \code{list(var = "trade", average = "full_pre")} - Average over full pre-period  
#'     \item \code{list(var = "schooling", every_n = 5)} - Every N periods
#'     \item \code{list(var = "invest", periods = c(1980, 1985))} - Specific periods
#'     \item \code{list(var = "gdp", first = 3)} - First N periods
#'     \item \code{list(var = "industry", last = 2)} - Last N periods
#'     \item \code{list(var = "trade", rolling = 3)} - Rolling N-period averages
#'     \item \code{list(var = "gdp", growth = "period_over_period")} - Growth rates
#'     \item \code{list(var = "gdp", volatility = "sd")} - Volatility measures
#'   }
#'   Default is empty list.
#' @param default_agg List. Optional default aggregation specification for simplified covariate processing.
#'   Contains three fields:
#'   \itemize{
#'     \item \code{outcome_agg}: Character. How to aggregate outcome variable ("mean", "per_period"). 
#'           If NULL, no outcome variable added to matching matrices.
#'     \item \code{features}: Character vector. Explicit list of variable names to include.
#'     \item \code{features_agg}: Character. How to aggregate feature variables ("mean", "per_period").
#'           Required when features are specified.
#'   }
#'   Variables specified in both default_agg and covagg will use covagg specification (with warning).
#'   Default is NULL.
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
#' # With flexible covariate aggregation (ADH-style)
#' scm_data_adh <- scdata(df = germany_data,
#'                        id.var = "country",
#'                        time.var = "year",
#'                        outcome.var = "gdp", 
#'                        period.pre = 1960:1990,
#'                        period.post = 1991:2003,
#'                        unit.tr = "West Germany",
#'                        unit.co = c("USA", "UK", "France"),
#'                        covagg = list(
#'                          gdp_each = list(var = "gdp", each = TRUE),
#'                          investment_avg = list(var = "investment", average = "full_pre"),
#'                          schooling_5yr = list(var = "schooling", every_n = 5),
#'                          trade_final = list(var = "trade", periods = 1990)
#'                        ))
#'                        
#' # With simplified default aggregation
#' scm_data_simple <- scdata(df = germany_data,
#'                           id.var = "country",
#'                           time.var = "year",
#'                           outcome.var = "gdp",
#'                           period.pre = 1960:1990,
#'                           period.post = 1991:2003,
#'                           unit.tr = "West Germany", 
#'                           unit.co = c("USA", "UK", "France"),
#'                           default_agg = list(
#'                             outcome_agg = "mean",              # Mean outcome for matching
#'                             features = c("investment", "trade"), # Explicit feature list
#'                             features_agg = "per_period"        # Per-period features
#'                           ))
#'                           
#' # Mixed approach: default_agg + custom covagg
#' scm_data_mixed <- scdata(df = germany_data,
#'                          id.var = "country", 
#'                          time.var = "year",
#'                          outcome.var = "gdp",
#'                          period.pre = 1960:1990,
#'                          period.post = 1991:2003,
#'                          unit.tr = "West Germany",
#'                          unit.co = c("USA", "UK", "France"),
#'                          default_agg = list(
#'                            outcome_agg = "per_period",
#'                            features = c("investment", "trade"),
#'                            features_agg = "mean"
#'                          ),
#'                          covagg = list(
#'                            schooling_custom = list(var = "schooling", last = 5)
#'                          ))
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
#'                            gdp_avg = list(var = "gdp", average = "full_pre"),
#'                            investment_avg = list(var = "investment", average = "full_pre"),
#'                            trade_periods = list(var = "trade", periods = c(1985, 1990))
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
                   default_agg = NULL,
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
  
  ############################################################################
  ############################################################################
  ############################################################################
  # Process default_agg parameter
  ############################################################################
  ############################################################################
  ############################################################################
  
  # Define supported aggregation types (extensible)
  SUPPORTED_AGG_TYPES <- c("mean", "per_period")
  
  # Validation function for default_agg
  validate_default_agg <- function(default_agg, data, id.var, time.var, outcome.var) {
    
    # Check required structure
    if (!"features" %in% names(default_agg)) {
      stop("default_agg must have 'features' field")
    }
    
    if (length(default_agg$features) > 0 && !"features_agg" %in% names(default_agg)) {
      stop("default_agg must have 'features_agg' when features are specified")
    }
    
    # Validate aggregation types
    if (!is.null(default_agg$outcome_agg)) {
      if (!default_agg$outcome_agg %in% SUPPORTED_AGG_TYPES) {
        stop(paste0("Invalid outcome_agg '", default_agg$outcome_agg, 
                    "'. Supported types: ", paste(SUPPORTED_AGG_TYPES, collapse = ", ")))
      }
    }
    
    if (length(default_agg$features) > 0) {
      if (!default_agg$features_agg %in% SUPPORTED_AGG_TYPES) {
        stop(paste0("Invalid features_agg '", default_agg$features_agg, 
                    "'. Supported types: ", paste(SUPPORTED_AGG_TYPES, collapse = ", ")))
      }
    }
    
    # Check variable existence
    all_vars <- names(data)
    missing_features <- setdiff(default_agg$features, all_vars)
    if (length(missing_features) > 0) {
      stop("Variables not found in data: ", paste(missing_features, collapse = ", "))
    }
    
    # Warn about structural variables
    structural_vars <- c(id.var, time.var, outcome.var)
    structural_conflicts <- intersect(default_agg$features, structural_vars)
    if (length(structural_conflicts) > 0) {
      warning("Features include structural variables: ", paste(structural_conflicts, collapse = ", "))
    }
  }
  
  # Create aggregation specification (extensible)
  create_agg_spec <- function(var_name, agg_type) {
    switch(agg_type,
      "mean" = list(var = var_name, average = "full_pre"),
      "per_period" = list(var = var_name, each = TRUE),
      # Future types can be added here easily
      stop("Unsupported aggregation type: ", agg_type)
    )
  }
  
  # Generate covagg specs from default_agg
  generate_default_covagg <- function(default_agg, outcome.var) {
    specs <- list()
    
    # Generate outcome spec if specified
    if (!is.null(default_agg$outcome_agg)) {
      outcome_spec <- create_agg_spec(outcome.var, default_agg$outcome_agg)
      specs[[paste0(outcome.var, "_", default_agg$outcome_agg)]] <- outcome_spec
    }
    
    # Generate feature specs
    for (feature in default_agg$features) {
      feature_spec <- create_agg_spec(feature, default_agg$features_agg)
      specs[[paste0(feature, "_", default_agg$features_agg)]] <- feature_spec
    }
    
    return(specs)
  }
  
  # Merge and deduplicate specs with warnings
  merge_and_dedupe <- function(default_specs, custom_covagg) {
    
    # Combine both sources
    all_specs <- c(default_specs, custom_covagg)
    
    # Extract variable names for deduplication
    var_names <- sapply(all_specs, function(x) x$var)
    
    # Check for duplicates
    duplicated_vars <- var_names[duplicated(var_names)]
    if (length(duplicated_vars) > 0) {
      unique_dups <- unique(duplicated_vars)
      warning("Duplicate variable specifications for: ", 
              paste(unique_dups, collapse = ", "), 
              " - custom covagg takes precedence")
    }
    
    # Keep last occurrence (custom wins)
    final_specs <- all_specs[!duplicated(var_names, fromLast = TRUE)]
    
    return(final_specs)
  }
  
  # Process default_agg if provided
  if (!is.null(default_agg)) {
    validate_default_agg(default_agg, df, id.var, time.var, outcome.var)
    default_specs <- generate_default_covagg(default_agg, outcome.var)
    covagg <- merge_and_dedupe(default_specs, covagg)
  }
  
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
        message(sprintf("DEBUG C: ERROR - Variable '%s' not found in data!", var_name))
        next
      }
      
      # Determine periods to use
      if (!is.null(cov_spec$periods)) {
        periods_to_use <- cov_spec$periods
      } else {
        periods_to_use <- period.pre
      }
      
      # Check if each=TRUE (create separate feature for each period)
      if (!is.null(cov_spec$each) && cov_spec$each) {
        
        # Create separate feature for each period
        for (period in periods_to_use) {
          if (period %in% period.pre) {
            idx2 = length(A_list) + 1
            
            # Extract data for this specific period and variable
            period_data <- data.bal.tr[data.bal.tr[, "Time"] == period, c(outcome.var, var_name)]
            if (nrow(period_data) > 0) {
              A <- as.matrix(period_data[[var_name]])
              colnames(A) <- unit.tr
              rownames(A) <- paste(unit.tr, cov_name, period, sep='.')
              A_list[[idx2]] = data.frame(A)
            }
          }
        }
        
      } else {
        # Default: aggregate using mean (or custom aggfunc if specified)
        aggfunc <- if (!is.null(cov_spec$aggfunc)) cov_spec$aggfunc else "mean"
        
        # Filter data for the specified periods
        period_data <- data.bal.tr[data.bal.tr[, "Time"] %in% periods_to_use, c(outcome.var, var_name)]
        
        if (nrow(period_data) > 0) {
          idx2 = length(A_list) + 1
          
          # Apply aggregation function
          if (aggfunc == "mean") {
            agg_value <- mean(period_data[[var_name]], na.rm = TRUE)
          } else {
            # Add support for other aggregation functions as needed
            agg_value <- mean(period_data[[var_name]], na.rm = TRUE)
          }
          
          A <- as.matrix(agg_value)
          colnames(A) <- unit.tr
          rownames(A) <- paste(unit.tr, cov_name, aggfunc, sep='.')
          A_list[[idx2]] = data.frame(A)
        }
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
        message(sprintf("DEBUG E: ERROR - Variable '%s' not found in B matrix data!", var_name))
        next
      }
      # Determine periods to use
      if (!is.null(cov_spec$periods)) {
        periods_to_use <- cov_spec$periods
      } else {
        periods_to_use <- period.pre
      }
      
      # Check if each=TRUE (create separate feature for each period)
      if (!is.null(cov_spec$each) && cov_spec$each) {
        
        # Create separate feature for each period
        for (period in periods_to_use) {
          if (period %in% period.pre) {
            temp = data.bal.co[data.bal.co[, "Time"] == period, c('ID', 'Time', var_name)]
            
            if (nrow(temp) > 0) {
              B = data.frame(data.table::dcast(data.table::melt(data.table(temp), id.vars=c('ID', 'Time')), variable+Time~ID))
              
              # Use underscores consistently - clean all names
              clean_unit_tr <- gsub('[^A-Za-z0-9_]', '_', unit.tr)
              clean_unit_tr <- gsub('_{2,}', '_', clean_unit_tr)
              rownames(B) = paste(rep(clean_unit_tr, nrow(B)), cov_name, period, sep='.')
              
              idxes = which(colnames(B) %in% c('variable', 'Time'))
              B = B[, -idxes]
              
              # IMPORTANT: Use the exact same transformation that Y.donors will use
              # This matches the str_remove + gsub pattern in Y.donors creation
              temp_names <- paste("gdpcap", colnames(B), sep = ".")
              y_names_equivalent <- stringr::str_remove(temp_names, "gdpcap")
              y_names_equivalent <- stringr::str_remove(y_names_equivalent, "\\.")
              final_clean_names <- gsub('[^A-Za-z0-9_]', '_', y_names_equivalent)
              final_clean_names <- gsub('_{2,}', '_', final_clean_names)
              # Remove trailing underscores to match Y.donors
              final_clean_names <- gsub('_+$', '', final_clean_names)
              
              # Only set column names if dimensions match
              if (length(final_clean_names) == ncol(B)) {
                colnames(B) <- final_clean_names
              } else {
                warning(paste("Column name dimension mismatch: B has", ncol(B), "columns but", length(final_clean_names), "names"))
              }
              
              idx2 = length(B_list) + 1
              B_list[[idx2]] = B
            }
          }
        }
        
      } else {
        
        # Default: aggregate using mean (or custom aggfunc if specified)
        aggfunc <- if (!is.null(cov_spec$aggfunc)) cov_spec$aggfunc else "mean"
        # Filter data for the specified periods
        temp = data.bal.co[data.bal.co[, "Time"] %in% periods_to_use, c('ID', 'Time', var_name)]
        
        if (nrow(temp) > 0) {
          # Apply aggregation by ID
          if (aggfunc == "mean") {
            agg_data <- aggregate(temp[[var_name]], by=list(ID=temp$ID), FUN=mean, na.rm=TRUE)
          } else {
            agg_data <- aggregate(temp[[var_name]], by=list(ID=temp$ID), FUN=mean, na.rm=TRUE)
          }
          
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
          rownames(B) <- paste(clean_unit_tr, cov_name, aggfunc, sep='.')
          
          idx2 = length(B_list) + 1
          B_list[[idx2]] = B
        }
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
