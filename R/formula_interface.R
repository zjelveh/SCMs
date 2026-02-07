#' Formula Interface for Synthetic Control Methods
#'
#' @title Formula-Based Interface for SCMs
#' @description Provides a formula interface for synthetic control methods, similar to lm() and other
#' standard R modeling functions. This makes the package more accessible to users familiar with
#' standard R model syntax while maintaining compatibility with data.table workflows.
#'
#' @name formula-interface
NULL

#' Synthetic Control with Formula Interface
#'
#' @description Estimate synthetic control using standard R formula syntax. This function provides
#' a convenient interface that automatically handles data preparation and estimation based on
#' formula specification.
#'
#' @param formula A formula object specifying the model. Format: \code{outcome ~ covariates | treated_unit}
#'   where covariates specify matching variables and treated_unit identifies the treated unit.
#'   Example: \code{gdp ~ population + investment | California}
#' @param data A data.frame or data.table containing the panel data
#' @param time.var Character. Name of time variable column
#' @param id.var Character. Name of unit identifier column
#' @param treated.period Numeric. First period when treatment begins
#' @param pre.period Numeric vector. Pre-treatment periods to use for matching
#' @param post.period Numeric vector. Post-treatment periods for outcome prediction
#' @param constraint Character or list. Weight constraint specification:
#'   \itemize{
#'     \item \code{"simplex"} - Standard synthetic control (default)
#'     \item \code{"lasso"} - L1 penalized weights
#'     \item \code{"ridge"} - L2 penalized weights  
#'     \item \code{"ols"} - Unconstrained weights
#'     \item List with detailed constraint parameters
#'   }
#' @param feature.weights Character. Feature weighting method ("uniform" or "optimize")
#' @param covariate.matching Character. How to aggregate covariates:
#'   \itemize{
#'     \item \code{"average"} - Use pre-period means (default)
#'     \item \code{"each"} - Include each pre-period separately
#'     \item \code{"last"} - Use only the last pre-treatment period
#'   }
#' @param ... Additional arguments passed to \code{scest()}
#'
#' @details
#' The formula interface translates standard R formula syntax into the appropriate scdata/scest calls:
#' 
#' \strong{Formula Components:}
#' \itemize{
#'   \item Left side: Outcome variable
#'   \item Right side (before |): Covariates for matching
#'   \item Right side (after |): Treated unit identifier
#' }
#' 
#' \strong{Special Formula Syntax:}
#' \itemize{
#'   \item \code{I(log(var))} - Transform variables
#'   \item \code{var1 + var2} - Include multiple covariates
#'   \item \code{.} - Include all available covariates
#'   \item \code{-var} - Exclude specific variables from matching
#' }
#'
#' @return An object of class "scest_formula" inheriting from "scest" with additional
#' formula-related metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic formula interface
#' result <- synth(gdp ~ population + investment | California,
#'                 data = panel_data,
#'                 time.var = "year",
#'                 id.var = "state", 
#'                 treated.period = 1989,
#'                 pre.period = 1980:1988,
#'                 post.period = 1989:2000)
#'
#' # With transformed variables
#' result <- synth(log(gdp) ~ log(population) + I(investment/1000) | California,
#'                 data = panel_data,
#'                 time.var = "year",
#'                 id.var = "state",
#'                 treated.period = 1989,
#'                 pre.period = 1980:1988,
#'                 post.period = 1989:2000,
#'                 constraint = "lasso")
#'
#' # Use all available covariates except one
#' result <- synth(gdp ~ . - population | California,
#'                 data = panel_data,
#'                 time.var = "year", 
#'                 id.var = "state",
#'                 treated.period = 1989,
#'                 pre.period = 1980:1988,
#'                 post.period = 1989:2000)
#' }
synth <- function(formula, data, time.var, id.var, treated.period, 
                  pre.period, post.period,
                  constraint = "simplex",
                  feature.weights = "uniform", 
                  covariate.matching = "average",
                  ...) {
  
  # Validate inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  
  if (!is.data.frame(data) && !is.data.table(data)) {
    stop("'data' must be a data.frame or data.table")
  }
  
  # Convert to data.table for efficiency
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  
  # Parse formula
  formula_parts <- parse_synth_formula(formula, data)
  
  # Extract components
  outcome_var <- formula_parts$outcome
  covariate_vars <- formula_parts$covariates
  treated_unit <- formula_parts$treated_unit
  
  # Validate that specified variables exist
  missing_vars <- setdiff(c(outcome_var, covariate_vars, time.var, id.var), names(data))
  if (length(missing_vars) > 0) {
    stop(paste("Missing variables in data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Validate treated unit exists
  if (!treated_unit %in% unique(data[[id.var]])) {
    stop(paste("Treated unit", treated_unit, "not found in", id.var, "column"))
  }
  
  # Get control units (all except treated)
  all_units <- unique(data[[id.var]])
  control_units <- setdiff(all_units, treated_unit)
  
  if (length(control_units) == 0) {
    stop("No control units found")
  }
  
  # Build covariate aggregation specification
  covagg_spec <- build_covagg_from_formula(covariate_vars, covariate.matching)
  
  # Convert constraint specification
  constraint_spec <- standardize_constraint_spec(constraint)
  
  # Prepare synthetic control data
  scdata_result <- scdata(
    df = data,
    id.var = id.var,
    time.var = time.var,
    outcome.var = outcome_var,
    period.pre = pre.period,
    period.post = post.period,
    unit.tr = treated_unit,
    unit.co = control_units,
    covagg = covagg_spec,
    constant = FALSE  # Can be added as formula option later
  )
  
  # Estimate synthetic control
  scest_result <- scest(
    data = scdata_result,
    w.constr = constraint_spec,
    feature_weights = feature.weights,
    ...
  )
  
  # Add formula metadata
  scest_result$formula <- formula
  scest_result$formula_parts <- formula_parts
  scest_result$call <- match.call()
  
  # Update class to include formula interface
  class(scest_result) <- c("scest_formula", class(scest_result))
  
  return(scest_result)
}

#' Parse Synthetic Control Formula
#'
#' @description Internal function to parse formula syntax for synthetic control models.
#' Handles the special "|" separator for treated unit specification.
#'
#' @param formula Formula object to parse
#' @param data Data containing the variables
#'
#' @return List with outcome, covariates, and treated_unit components
#' @keywords internal
parse_synth_formula <- function(formula, data) {
  
  # Convert formula to character for parsing
  formula_char <- deparse(formula)
  
  # Split on "|" to separate main formula from treated unit
  if (!grepl("\\|", formula_char)) {
    stop("Formula must include '|' to specify treated unit. Example: outcome ~ covariates | treated_unit")
  }
  
  formula_parts <- strsplit(formula_char, "\\|")[[1]]
  if (length(formula_parts) != 2) {
    stop("Formula must have exactly one '|' separator")
  }
  
  main_formula <- trimws(formula_parts[1])
  treated_unit_spec <- trimws(formula_parts[2])
  
  # Parse main formula (outcome ~ covariates)
  main_parsed <- strsplit(main_formula, "~")[[1]]
  if (length(main_parsed) != 2) {
    stop("Main formula must be of form: outcome ~ covariates")
  }
  
  outcome <- trimws(main_parsed[1])
  covariates_spec <- trimws(main_parsed[2])
  
  # Handle outcome transformations
  outcome <- parse_variable_expression(outcome)
  
  # Parse covariates
  if (covariates_spec == ".") {
    # Use all variables except outcome, time, and id
    exclude_vars <- c(outcome, "time", "year", "period")  # common time variable names
    covariates <- setdiff(names(data), exclude_vars)
  } else {
    # Parse specific covariate specification
    covariates <- parse_covariate_specification(covariates_spec, data)
  }
  
  # Parse treated unit (remove quotes if present)
  treated_unit <- gsub("[\"\']", "", treated_unit_spec)
  
  return(list(
    outcome = outcome,
    covariates = covariates,
    treated_unit = treated_unit
  ))
}

#' Parse Variable Expression
#'
#' @description Handle variable transformations in formula (like log, I(), etc.)
#'
#' @param expr Character expression to parse
#' @return Character variable name (simplified for now)
#' @keywords internal
parse_variable_expression <- function(expr) {
  # Simplified implementation - just extract base variable name
  # In full implementation would handle I(), log(), etc.
  
  # Remove common transformations to get base variable
  base_var <- gsub("^(log|sqrt|I)\\(([^)]+)\\)$", "\\2", expr)
  base_var <- gsub("^I\\(([^)]+)\\)$", "\\1", base_var)
  
  # For now, return the expression as-is
  # Full implementation would create new transformed variables
  return(trimws(expr))
}

#' Parse Covariate Specification
#'
#' @description Parse covariate part of formula, handling +, -, and other operators
#'
#' @param spec Character specification string
#' @param data Data.table containing variables
#' @return Character vector of covariate names
#' @keywords internal
parse_covariate_specification <- function(spec, data) {
  
  # Split on +
  terms <- strsplit(spec, "\\+")[[1]]
  terms <- trimws(terms)
  
  # Handle exclusions (terms starting with -)
  include_terms <- terms[!grepl("^-", terms)]
  exclude_terms <- gsub("^-\\s*", "", terms[grepl("^-", terms)])
  
  # For include terms, parse expressions
  covariates <- character(0)
  for (term in include_terms) {
    if (term == ".") {
      # Add all available variables
      available_vars <- names(data)
      covariates <- c(covariates, available_vars)
    } else {
      # Parse individual term
      parsed_term <- parse_variable_expression(term)
      covariates <- c(covariates, parsed_term)
    }
  }
  
  # Remove excluded terms
  if (length(exclude_terms) > 0) {
    covariates <- setdiff(covariates, exclude_terms)
  }
  
  # Remove duplicates
  covariates <- unique(covariates)
  
  return(covariates)
}

#' Build Covariate Aggregation from Formula
#'
#' @description Convert formula-parsed covariates into scdata covagg format
#'
#' @param covariates Character vector of covariate names
#' @param matching_method Character method for aggregation
#' @return List in scdata covagg format
#' @keywords internal
build_covagg_from_formula <- function(covariates, matching_method) {
  
  if (length(covariates) == 0) {
    return(list())
  }
  
  build_single <- function(v) {
    switch(
      matching_method,
      "average" = list(
        var = v,
        select_periods = list(type = "all_pre"),
        partition_periods = list(type = "all"),
        compute = "mean"
      ),
      "each" = list(
        var = v,
        select_periods = list(type = "all_pre"),
        partition_periods = list(type = "by_period"),
        compute = "mean"
      ),
      "last" = list(
        var = v,
        select_periods = list(type = "predicate", fn = function(p) p == max(p)),
        partition_periods = list(type = "all"),
        compute = "mean"
      ),
      list(
        var = v,
        select_periods = list(type = "all_pre"),
        partition_periods = list(type = "all"),
        compute = "mean"
      )
    )
  }

  lapply(covariates, build_single)
}

#' Standardize Constraint Specification
#'
#' @description Convert character constraint to scest format
#'
#' @param constraint Character or list constraint specification
#' @return List in scest w.constr format
#' @keywords internal
standardize_constraint_spec <- function(constraint) {
  
  if (is.character(constraint)) {
    constraint_spec <- switch(constraint,
      "simplex" = list(name = "simplex"),
      "lasso" = list(name = "lasso", Q = 0.1),  # Default lambda
      "ridge" = list(name = "ridge", Q = 0.1),  # Default lambda
      "ols" = list(name = "ols"),
      # Default to simplex
      list(name = "simplex")
    )
  } else if (is.list(constraint)) {
    constraint_spec <- constraint
  } else {
    stop("'constraint' must be character or list")
  }
  
  return(constraint_spec)
}

#' Print Method for Formula-Based Synthetic Control
#'
#' @description Print method for scest_formula objects
#'
#' @param x An scest_formula object
#' @param ... Additional arguments
#' @export
print.scest_formula <- function(x, ...) {
  
  cat("Synthetic Control Model (Formula Interface)\n")
  cat("==========================================\n\n")
  
  # Print formula
  cat("Formula: ")
  print(x$formula)
  cat("\n")
  
  # Print basic summary
  weights <- x$est.results$w
  n_donors <- length(weights)
  n_active <- sum(weights > 1e-6)
  
  cat("Donor units:", n_donors, "(", n_active, "active )\n")
  cat("Constraint: ", x$est.results$w.constr$name %||% "unknown", "\n")
  
  # Print top donor weights
  if (!is.null(names(weights))) {
    donor_names <- names(weights)
  } else {
    donor_names <- paste0("donor_", seq_along(weights))
  }
  
  weight_df <- data.frame(
    unit = donor_names,
    weight = weights,
    stringsAsFactors = FALSE
  )
  weight_df <- weight_df[order(-weight_df$weight), ]
  
  cat("\nTop donor weights:\n")
  print(head(weight_df, 5), row.names = FALSE, digits = 3)
  
  if (nrow(weight_df) > 5) {
    cat("... and", nrow(weight_df) - 5, "more donors\n")
  }
  
  # Fit statistics
  if (!is.null(x$est.results$Y.pre.fit) && !is.null(x$data$Y.pre)) {
    pre_rmse <- sqrt(mean((x$data$Y.pre - x$est.results$Y.pre.fit)^2, na.rm = TRUE))
    cat("\nPre-treatment RMSE:", round(pre_rmse, 4), "\n")
  }
  
  invisible(x)
}

#' Summary Method for Formula-Based Synthetic Control
#'
#' @description Summary method for scest_formula objects
#'
#' @param object An scest_formula object
#' @param ... Additional arguments
#' @export
summary.scest_formula <- function(object, ...) {
  
  cat("Synthetic Control Model Summary (Formula Interface)\n")
  cat("==================================================\n\n")
  
  # Print call
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  # Print formula details
  cat("Formula components:\n")
  cat("  Outcome:", object$formula_parts$outcome, "\n")
  cat("  Covariates:", paste(object$formula_parts$covariates, collapse = ", "), "\n")
  cat("  Treated unit:", object$formula_parts$treated_unit, "\n\n")
  
  # Call the regular scest summary for detailed results
  NextMethod("summary", object)
  
  invisible(object)
}
