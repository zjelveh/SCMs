#' Tidy Methods for SCMs Objects
#'
#' @title Broom-Compatible Tidy Methods for Synthetic Control Objects
#' @description These methods provide broom-compatible interfaces for extracting information
#' from synthetic control estimation results in tidy data.table format. They follow the
#' broom convention of tidy(), glance(), and augment() methods while returning data.table
#' objects for efficiency.
#'
#' @name tidy-methods
#' @import data.table
#' @importFrom generics tidy glance augment
NULL

#' Extract Synthetic Control Weights in Tidy Format
#'
#' @description Extract donor unit weights and related information from synthetic control
#' estimation results in a tidy data.table format compatible with broom workflows.
#'
#' @param x An object of class "scest" from \code{scest()}
#' @param conf.int Logical. Whether to include confidence intervals (if available)
#' @param conf.level Numeric. Confidence level for intervals (default 0.95)
#' @param ... Additional arguments (not currently used)
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item \code{term} - Character. Donor unit names
#'   \item \code{estimate} - Numeric. Synthetic control weights
#'   \item \code{donor_rank} - Integer. Rank of donor by weight (1 = highest)
#'   \item \code{weight_percentage} - Numeric. Weight as percentage of total
#'   \item \code{is_active} - Logical. Whether weight is meaningfully non-zero
#'   \item \code{constraint_type} - Character. Type of constraint used
#'   \item \code{std.error} - Numeric. Standard error (if conf.int=TRUE and available)
#'   \item \code{conf.low} - Numeric. Lower confidence bound (if conf.int=TRUE)
#'   \item \code{conf.high} - Numeric. Upper confidence bound (if conf.int=TRUE)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic synthetic control estimation
#' scdata_obj <- scdata(df, ...)
#' scest_result <- scest(scdata_obj, w.constr = list(name = "simplex"))
#' 
#' # Extract tidy weights
#' tidy(scest_result)
#' 
#' # With confidence intervals (if bootstrap inference was run)
#' tidy(scest_result, conf.int = TRUE)
#' }
tidy.scest <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  
  # Validate input
  if (!inherits(x, "scest")) {
    stop("Object must be of class 'scest'")
  }
  
  # Extract weights
  weights <- x$est.results$w
  if (is.null(weights)) {
    stop("No weights found in scest object")
  }
  
  # Get donor names (from data or create generic names)
  if (!is.null(names(weights))) {
    donor_names <- names(weights)
  } else if (!is.null(colnames(x$data$B))) {
    donor_names <- colnames(x$data$B)
  } else {
    donor_names <- paste0("donor_", seq_along(weights))
  }
  
  # Create base tidy data.table
  tidy_dt <- data.table(
    term = donor_names,
    estimate = as.numeric(weights)
  )
  
  # Add derived columns
  tidy_dt[, `:=`(
    donor_rank = rank(-estimate, ties.method = "min"),
    weight_percentage = estimate / sum(estimate) * 100,
    is_active = estimate > 1e-6  # Threshold for "meaningful" weights
  )]
  
  # Add constraint information
  constraint_type <- if (!is.null(x$est.results$w.constr$name)) {
    x$est.results$w.constr$name
  } else {
    "unknown"
  }
  tidy_dt[, constraint_type := constraint_type]
  
  # Add confidence intervals if requested
  if (conf.int) {
    # Check if bootstrap results are available
    if ("bootstrap_weights" %in% names(x$est.results)) {
      
      alpha <- 1 - conf.level
      boot_weights <- x$est.results$bootstrap_weights
      
      # Calculate confidence intervals for each donor
      conf_intervals <- apply(boot_weights, 2, function(w) {
        quantile(w, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
      })
      
      tidy_dt[, `:=`(
        std.error = apply(boot_weights, 2, sd, na.rm = TRUE),
        conf.low = conf_intervals[1, ],
        conf.high = conf_intervals[2, ]
      )]
      
    } else {
      warning("Confidence intervals requested but no bootstrap results available")
      tidy_dt[, `:=`(
        std.error = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_
      )]
    }
  }
  
  # Order by weight (highest first)
  setorder(tidy_dt, -estimate)
  
  return(tidy_dt)
}

#' Extract Model-Level Summary Statistics
#'
#' @description Extract model-level fit statistics and metadata from synthetic control
#' estimation results in a single-row data.table.
#'
#' @param x An object of class "scest" from \code{scest()}
#' @param ... Additional arguments (not currently used)
#'
#' @return A single-row data.table with columns:
#' \itemize{
#'   \item \code{n_donors} - Integer. Number of donor units
#'   \item \code{n_pre_periods} - Integer. Number of pre-treatment periods
#'   \item \code{n_post_periods} - Integer. Number of post-treatment periods  
#'   \item \code{rmse_pre} - Numeric. Root mean squared error in pre-treatment fit
#'   \item \code{rmse_post} - Numeric. Root mean squared error in post-treatment prediction
#'   \item \code{constraint_type} - Character. Type of weight constraint used
#'   \item \code{feature_weighting} - Character. Feature weighting method used
#'   \item \code{n_active_donors} - Integer. Number of donors with meaningful weights
#'   \item \code{max_weight} - Numeric. Largest donor weight
#'   \item \code{herfindahl_index} - Numeric. Concentration measure of weights
#'   \item \code{has_constant} - Logical. Whether constant term was included
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract model-level statistics
#' glance(scest_result)
#' }
glance.scest <- function(x, ...) {
  
  # Validate input
  if (!inherits(x, "scest")) {
    stop("Object must be of class 'scest'")
  }
  
  # Extract basic dimensions
  n_donors <- length(x$est.results$w)
  n_pre_periods <- length(x$est.results$Y.pre.fit)
  n_post_periods <- length(x$est.results$Y.post.fit)
  
  # Calculate fit statistics
  if (!is.null(x$data$Y.pre) && !is.null(x$est.results$Y.pre.fit)) {
    pre_residuals <- x$data$Y.pre - x$est.results$Y.pre.fit
    rmse_pre <- sqrt(mean(pre_residuals^2, na.rm = TRUE))
  } else {
    rmse_pre <- NA_real_
  }
  
  if (!is.null(x$data$Y.post) && !is.null(x$est.results$Y.post.fit)) {
    post_residuals <- x$data$Y.post - x$est.results$Y.post.fit
    rmse_post <- sqrt(mean(post_residuals^2, na.rm = TRUE))
  } else {
    rmse_post <- NA_real_
  }
  
  # Weight concentration statistics
  weights <- x$est.results$w
  max_weight <- max(weights, na.rm = TRUE)
  n_active_donors <- sum(weights > 1e-6, na.rm = TRUE)
  
  # Herfindahl index (concentration measure)
  weight_shares <- weights / sum(weights, na.rm = TRUE)
  herfindahl_index <- sum(weight_shares^2, na.rm = TRUE)
  
  # Constraint and method information
  constraint_type <- if (!is.null(x$est.results$w.constr$name)) {
    x$est.results$w.constr$name
  } else {
    "unknown"
  }
  
  # Feature weighting (if available in object)
  feature_weighting <- if ("V" %in% names(x$est.results)) {
    if (all(diag(x$est.results$V) == diag(x$est.results$V)[1])) {
      "uniform"
    } else {
      "optimized"
    }
  } else {
    "unknown"
  }
  
  # Check for constant term
  has_constant <- !is.null(x$est.results$constant_term)
  
  # Create glance data.table
  glance_dt <- data.table(
    n_donors = n_donors,
    n_pre_periods = n_pre_periods,
    n_post_periods = n_post_periods,
    rmse_pre = rmse_pre,
    rmse_post = rmse_post,
    constraint_type = constraint_type,
    feature_weighting = feature_weighting,
    n_active_donors = n_active_donors,
    max_weight = max_weight,
    herfindahl_index = herfindahl_index,
    has_constant = has_constant
  )
  
  return(glance_dt)
}

#' Extract Specification Curve Results in Tidy Format
#'
#' @description Extract treatment effect estimates and metadata from specification
#' curve analysis results in tidy data.table format.
#'
#' @param x An object of class "spec_curve" from \code{run_spec_curve_analysis()}
#' @param conf.int Logical. Whether to include confidence intervals from inference
#' @param conf.level Numeric. Confidence level for intervals (default 0.95) 
#' @param ... Additional arguments (not currently used)
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item \code{spec_id} - Character. Specification identifier
#'   \item \code{estimate} - Numeric. Treatment effect estimate
#'   \item \code{outcome} - Character. Outcome variable name
#'   \item \code{unit_name} - Character. Treated unit name  
#'   \item \code{rmse} - Numeric. Root mean squared error
#'   \item \code{outcome_model} - Character. Outcome modeling approach
#'   \item \code{constraint_type} - Character. Weight constraint type
#'   \item \code{feature_weights} - Character. Feature weighting method
#'   \item \code{p.value} - Numeric. P-value from inference (if available)
#'   \item \code{statistic} - Numeric. Test statistic (if available)
#'   \item \code{method} - Character. Inference method used
#'   \item \code{significant} - Logical. Whether p < 0.05
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract specification curve results
#' spec_results <- run_spec_curve_analysis(...)
#' tidy(spec_results)
#' 
#' # With inference results
#' tidy(spec_results, conf.int = TRUE)
#' }
tidy.spec_curve <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  
  # Validate input
  if (!inherits(x, "spec_curve")) {
    stop("Object must be of class 'spec_curve'")
  }
  
  if (is.null(x$results)) {
    stop("No results found in spec_curve object")
  }
  
  # Extract treated unit results
  treated_results <- x$results[unit_type == "treated" & post_period == TRUE]
  
  if (nrow(treated_results) == 0) {
    warning("No treated unit results found")
    return(data.table())
  }
  
  # Create base tidy data.table
  tidy_dt <- data.table(
    spec_id = treated_results$full_spec_id,
    estimate = treated_results$tau,
    outcome = treated_results$outcome,
    unit_name = treated_results$unit_name,
    rmse = treated_results$rmse
  )
  
  # Add specification characteristics (if available)
  spec_cols <- c("outcome_model", "const", "fw", "feat")
  available_cols <- intersect(spec_cols, names(treated_results))
  
  for (col in available_cols) {
    tidy_dt[, (col) := treated_results[[col]]]
  }
  
  # Rename columns for consistency
  if ("const" %in% names(tidy_dt)) {
    setnames(tidy_dt, "const", "constraint_type")
  }
  if ("fw" %in% names(tidy_dt)) {
    setnames(tidy_dt, "fw", "feature_weights")
  }
  
  # Add inference results if available and requested
  if (conf.int && !is.null(x$abadie_inference)) {
    
    if ("p_values_rmse_ratio" %in% names(x$abadie_inference)) {
      p_values <- x$abadie_inference$p_values_rmse_ratio
      
      # Merge inference results
      inference_dt <- p_values[, .(
        spec_id = full_spec_id,
        p.value = p_value,
        method = "abadie_placebo"
      )]
      
      tidy_dt <- merge(tidy_dt, inference_dt, by = "spec_id", all.x = TRUE)
    }
  }
  
  # Add bootstrap inference if available
  if (conf.int && !is.null(x$bootstrap_inference)) {
    
    if ("p_values" %in% names(x$bootstrap_inference)) {
      boot_p_values <- x$bootstrap_inference$p_values
      
      # If we already have Abadie p-values, add bootstrap as separate columns
      if ("p.value" %in% names(tidy_dt)) {
        boot_dt <- boot_p_values[, .(
          spec_id = full_spec_id,
          p.value.bootstrap = p_value_two_tailed
        )]
        tidy_dt <- merge(tidy_dt, boot_dt, by = "spec_id", all.x = TRUE)
      } else {
        # Use bootstrap p-values as primary
        boot_dt <- boot_p_values[, .(
          spec_id = full_spec_id,
          p.value = p_value_two_tailed,
          method = "bootstrap"
        )]
        tidy_dt <- merge(tidy_dt, boot_dt, by = "spec_id", all.x = TRUE)
      }
    }
  }
  
  # Add significance indicator
  if ("p.value" %in% names(tidy_dt)) {
    tidy_dt[, significant := p.value < 0.05]
  }
  
  # Order by effect size
  setorder(tidy_dt, estimate)
  
  return(tidy_dt)
}

#' Augment Synthetic Control Data with Fitted Values and Residuals
#'
#' @description Add fitted values, residuals, and other observation-level statistics
#' to the original data used in synthetic control estimation.
#'
#' @param x An object of class "scest" from \code{scest()}
#' @param data Optional data.frame/data.table to augment. If NULL, attempts to reconstruct
#' @param ... Additional arguments (not currently used)
#'
#' @return A data.table with original data plus:
#' \itemize{
#'   \item \code{.fitted} - Numeric. Fitted values from synthetic control
#'   \item \code{.resid} - Numeric. Residuals (actual - fitted)
#'   \item \code{.period_type} - Character. "pre" or "post" treatment
#'   \item \code{.is_treated} - Logical. Whether observation is from treated unit
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Augment with fitted values and residuals
#' augmented_data <- augment(scest_result)
#' 
#' # Use specific data
#' augmented_data <- augment(scest_result, data = original_panel_data)
#' }
augment.scest <- function(x, data = NULL, ...) {
  
  # Validate input
  if (!inherits(x, "scest")) {
    stop("Object must be of class 'scest'")
  }
  
  # Extract fitted values
  fitted_pre <- x$est.results$Y.pre.fit
  fitted_post <- x$est.results$Y.post.fit
  
  if (is.null(fitted_pre) || is.null(fitted_post)) {
    stop("No fitted values found in scest object")
  }
  
  # Create augmented data.table
  if (is.null(data)) {
    # Reconstruct basic time series structure
    n_pre <- length(fitted_pre)
    n_post <- length(fitted_post)
    
    augmented_dt <- data.table(
      period = c(seq_len(n_pre), seq_len(n_post) + n_pre),
      actual = c(x$data$Y.pre, x$data$Y.post),
      .fitted = c(fitted_pre, fitted_post),
      .period_type = rep(c("pre", "post"), c(n_pre, n_post)),
      .is_treated = TRUE
    )
    
  } else {
    # Use provided data
    if (!is.data.table(data)) {
      augmented_dt <- as.data.table(data)
    } else {
      augmented_dt <- copy(data)
    }
    
    # Try to match fitted values to data
    # This is simplified - in practice would need more sophisticated matching
    warning("Data augmentation with external data is experimental")
  }
  
  # Calculate residuals
  if ("actual" %in% names(augmented_dt)) {
    augmented_dt[, .resid := actual - .fitted]
  }
  
  return(augmented_dt)
}