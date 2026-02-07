#' S3 Methods for Specification Curve Objects
#'
#' @description
#' This file contains S3 methods for objects of class "spec_curve" created by
#' the \code{\link{spec_curve}} function. These methods provide a professional
#' interface for examining and summarizing specification curve analysis results.
#'
#' @name spec_curve-methods
NULL

#' Print Method for Specification Curve Results
#'
#' @description Print method for specification curve objects that provides a concise
#' overview of the analysis results including number of specifications, treatment
#' effect distribution, and available inference results.
#'
#' @rdname spec_curve-methods
#' @param x An object of class "spec_curve" from \code{\link{spec_curve}}
#' @param ... Additional arguments (currently unused)
#'
#' @return No return value, called for side effects (printing)
#'
#' @details
#' The print method displays:
#' \itemize{
#'   \item Number of specifications estimated
#'   \item Number of outcomes analyzed  
#'   \item Treatment effect summary statistics
#'   \item Available inference methods and significance rates
#'   \item Specification dimension summary
#' }
#'
#' @examples
#' \dontrun{
#' # After running specification curve analysis
#' spec_results <- spec_curve(...)
#' print(spec_results)  # or just: spec_results
#' }
#'
#' @export
print.spec_curve <- function(x, ...) {
  cat("\n")
  cat("Specification Curve Analysis Results\n")
  cat("====================================\n\n")
  
  # Basic information
  results_data <- x$results
  n_specs <- length(unique(results_data$full_spec_id))
  n_outcomes <- length(unique(results_data$outcome))
  n_units <- length(unique(results_data$unit_name))
  
  cat("Specifications:       ", n_specs, "\n")
  cat("Outcomes analyzed:    ", n_outcomes, "\n")
  cat("Total observations:   ", nrow(results_data), "\n")
  cat("Units:                ", n_units, "\n")
  
  # Treatment effect summary for treated units only
  treated_effects <- results_data[results_data$unit_type == "treated" & 
                                 results_data$post_period == TRUE, ]
  
  if (nrow(treated_effects) > 0) {
    cat("\nTreatment Effect Summary:\n")
    cat("  Mean:               ", round(mean(treated_effects$tau, na.rm = TRUE), 4), "\n")
    cat("  Median:             ", round(median(treated_effects$tau, na.rm = TRUE), 4), "\n") 
    cat("  Std. Dev:           ", round(sd(treated_effects$tau, na.rm = TRUE), 4), "\n")
    cat("  Range:              [", round(min(treated_effects$tau, na.rm = TRUE), 4), 
        ", ", round(max(treated_effects$tau, na.rm = TRUE), 4), "]\n")
  }
  
  # Inference information
  has_abadie <- !is.null(x$abadie_inference)
  has_bootstrap <- !is.null(x$bootstrap_inference)
  
  if (has_abadie || has_bootstrap) {
    cat("\nInference Methods:\n")
    
    if (has_abadie) {
      cat("  Abadie placebo:     Available\n")
      
      # Report significance rates if p-values available
      if ("p_values_rmse_ratio" %in% names(x$abadie_inference)) {
        treated_pvals <- x$abadie_inference$p_values_rmse_ratio[
          x$abadie_inference$p_values_rmse_ratio$unit_type == "treated", ]
        if (nrow(treated_pvals) > 0) {
          sig_rate <- mean(treated_pvals$p_value < 0.05, na.rm = TRUE) * 100
          cat("    Significant (p<0.05): ", round(sig_rate, 1), "%\n")
        }
      }
    }
    
    if (has_bootstrap) {
      cat("  Bootstrap:          Available\n")
      
      # Report significance rates if p-values available  
      if ("p_values" %in% names(x$bootstrap_inference)) {
        sig_rate <- mean(x$bootstrap_inference$p_values$p_value_two_tailed < 0.05, na.rm = TRUE) * 100
        cat("    Significant (p<0.05): ", round(sig_rate, 1), "%\n")
      }
    }
  } else {
    cat("\nInference Methods:    None\n")
  }
  
  # Specification dimensions
  cat("\nSpecification Dimensions:\n")
  if ("outcome_model" %in% names(results_data)) {
    n_outcome_models <- length(unique(results_data$outcome_model))
    cat("  Outcome models:     ", n_outcome_models, "\n")
  }
  if ("const" %in% names(results_data)) {
    n_constraints <- length(unique(results_data$const))
    cat("  Constraint types:   ", n_constraints, "\n")
  }
  if ("fw" %in% names(results_data)) {
    n_fw <- length(unique(results_data$fw))
    cat("  Feature weights:    ", n_fw, "\n")
  }
  
  cat("\n")
  cat("Use summary() for detailed breakdown\n")
  cat("Use plot() to visualize specification curve\n")
  cat("\n")
}

#' Summary Method for Specification Curve Results
#'
#' @description Detailed summary method for specification curve objects that provides
#' comprehensive information about the analysis results, including breakdown by
#' specification dimensions and detailed inference results.
#'
#' @param object An object of class "spec_curve" from \code{\link{spec_curve}}
#' @param ... Additional arguments (currently unused)
#'
#' @return No return value, called for side effects (printing detailed summary)
#'
#' @details
#' The summary method provides detailed information including:
#' \itemize{
#'   \item Treatment effect statistics by specification dimension
#'   \item Inference results breakdown
#'   \item Model performance statistics (RMSE, fit)
#'   \item Specification robustness indicators
#' }
#'
#' @examples
#' \dontrun{
#' spec_results <- spec_curve(...)
#' summary(spec_results)
#' }
#'
#' @export
summary.spec_curve <- function(object, ...) {
  print(object)
  
  results_data <- object$results
  treated_effects <- results_data[results_data$unit_type == "treated" & 
                                 results_data$post_period == TRUE, ]
  
  if (nrow(treated_effects) == 0) {
    cat("No treated unit effects found for detailed summary.\n")
    return(invisible(object))
  }
  
  cat("\nDetailed Specification Breakdown:\n")
  cat("=================================\n")
  
  # Effect size by outcome model
  if ("outcome_model" %in% names(treated_effects)) {
    cat("\nTreatment Effects by Outcome Model:\n")
    model_summary <- aggregate(tau ~ outcome_model, data = treated_effects, 
                              FUN = function(x) c(mean = mean(x), 
                                                 median = median(x),
                                                 sd = sd(x),
                                                 n = length(x)))
    model_df <- do.call(data.frame, model_summary)
    names(model_df) <- c("Model", "Mean", "Median", "Std.Dev", "N.Specs")
    model_df[,2:4] <- round(model_df[,2:4], 4)
    print(model_df, row.names = FALSE)
  }
  
  # Effect size by constraint type
  if ("const" %in% names(treated_effects)) {
    cat("\nTreatment Effects by Constraint Type:\n")
    constraint_summary <- aggregate(tau ~ const, data = treated_effects,
                                   FUN = function(x) c(mean = mean(x),
                                                      median = median(x), 
                                                      sd = sd(x),
                                                      n = length(x)))
    constraint_df <- do.call(data.frame, constraint_summary)
    names(constraint_df) <- c("Constraint", "Mean", "Median", "Std.Dev", "N.Specs")
    constraint_df[,2:4] <- round(constraint_df[,2:4], 4)
    print(constraint_df, row.names = FALSE)
  }
  
  # Model fit information
  if ("rmse" %in% names(treated_effects)) {
    cat("\nModel Fit Statistics:\n")
    cat("  Mean RMSE:          ", round(mean(treated_effects$rmse, na.rm = TRUE), 4), "\n")
    cat("  Median RMSE:        ", round(median(treated_effects$rmse, na.rm = TRUE), 4), "\n")
    cat("  RMSE Range:         [", round(min(treated_effects$rmse, na.rm = TRUE), 4),
        ", ", round(max(treated_effects$rmse, na.rm = TRUE), 4), "]\n")
  }
  
  # Robustness indicators
  positive_effects <- sum(treated_effects$tau > 0, na.rm = TRUE)
  total_effects <- sum(!is.na(treated_effects$tau))
  
  cat("\nRobustness Indicators:\n")
  cat("  Positive effects:   ", positive_effects, " (", 
      round(positive_effects / total_effects * 100, 1), "%)\n")
  cat("  Negative effects:   ", total_effects - positive_effects, " (",
      round((total_effects - positive_effects) / total_effects * 100, 1), "%)\n")
  
  # Effect size consistency
  if (total_effects > 1) {
    cv <- sd(treated_effects$tau, na.rm = TRUE) / abs(mean(treated_effects$tau, na.rm = TRUE))
    cat("  Coeff. of variation:", round(cv, 3), "\n")
  }
  
  cat("\n")
  invisible(object)
}

#' Plot Method for Specification Curve Results
#'
#' @description Plot method that creates a specification curve visualization by calling
#' the main plotting function with appropriate defaults.
#'
#' @param x An object of class "spec_curve" from \code{\link{spec_curve}}
#' @param ... Additional arguments passed to \code{\link{plot_spec_curve}}
#'
#' @return A ggplot object containing the specification curve plot
#'
#' @details
#' This is a convenience method that calls \code{\link{plot_spec_curve}} with the
#' specification curve data. All arguments supported by \code{plot_spec_curve}
#' can be passed through the \code{...} parameter.
#'
#' @examples
#' \dontrun{
#' spec_results <- spec_curve(...)
#' plot(spec_results)
#' 
#' # With additional options
#' plot(spec_results, show_shap = TRUE, show_pvalues = TRUE)
#' }
#'
#' @seealso \code{\link{plot_spec_curve}} for detailed plotting options
#'
#' @export
plot.spec_curve <- function(x, ...) {
  # Extract treated unit name - look for the treated unit in the data
  treated_rows <- x$results[x$results$unit_type == "treated"]
  treated_units <- unique(treated_rows$unit_name)
  
  if (length(treated_units) == 0) {
    stop("No treated unit found in specification curve results")
  }
  
  if (length(treated_units) > 1) {
    warning("Multiple treated units found. Using the first one: ", treated_units[1])
  }
  
  name_treated_unit <- treated_units[1]
  
  # Call the main plotting function
  plot_spec_curve(long_data = x, 
                  name_treated_unit = name_treated_unit,
                  ...)
}
