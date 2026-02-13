#' Run Comprehensive Specification Curve Analysis
#'
#' @title Run Specification Curve Analysis with Multiple Configurations
#' @description Performs specification curve analysis across multiple modeling choices
#' and configurations, providing a systematic exploration of researcher degrees of freedom.
#'
#' @param dataset Data frame containing the panel data for analysis.
#' @param params List containing ALL analysis parameters including:
#'   \itemize{
#'     \item \code{outcomes} - Character vector of outcome variables (REQUIRED)
#'     \item \code{col_name_unit_name} - Column name for unit identifiers
#'     \item \code{name_treated_unit} - Name of treated unit
#'     \item \code{covagg} - List of covariate specifications
#'     \item \code{treated_period} - Treatment start period
#'     \item \code{min_period} - Minimum period for analysis
#'     \item \code{end_period} - Maximum period for analysis
#'     \item \code{col_name_period} - Column name for time periods
#'     \item \code{feature_weights} - Feature weighting methods
#'     \item \code{donor_sample} - Donor sample selection methods
#'     \item \code{outcome_models} - Outcome modeling approaches
#'     \item \code{constraints} - Weight constraint specifications
#'     \item \code{inference_type} - "placebo", "bootstrap", or "all" (default: "placebo")
#'     \item \code{inference_config} - List with bootstrap_n_replications, verbose, etc.
#'     \item \code{expected_direction} - "negative", "positive", or "two_sided" (default: "negative")
#'   }
#' @param cores Integer. Number of cores for parallel processing (applies to both bootstrap and placebo inference unless overridden in params$inference_config).
#' @param output_dir Character. Directory to save results. If NULL, results are not saved.
#'
#' @return List structure containing:
#'   \itemize{
#'     \item \code{results} - Data.table in long format with all specification results
#'     \item \code{abadie_inference} - Abadie placebo inference results (if inference_type includes "placebo")  
#'     \item \code{bootstrap_inference} - Bootstrap inference results (if inference_type includes "bootstrap")
#'   }
#'
#' @export
#'
#' @examples
#' \donttest{
#' toy <- data.frame(
#'   unit = rep(c("treated", "donor1", "donor2"), each = 6),
#'   year = rep(1:6, 3),
#'   y = c(10, 11, 12, 14, 15, 16,
#'         9, 10, 11, 12, 13, 14,
#'         11, 12, 13, 13, 14, 15)
#' )
#' params <- list(
#'   outcomes = "y",
#'   col_name_unit_name = "unit",
#'   name_treated_unit = "treated",
#'   covagg = list(
#'     baseline = list(
#'       label = "baseline",
#'       operations = list(list(var = "outcome_var", partition_periods = list(type = "by_period")))
#'     )
#'   ),
#'   treated_period = 5,
#'   min_period = 1,
#'   end_period = 6,
#'   col_name_period = "year",
#'   feature_weights = "uniform",
#'   donor_sample = "all",
#'   outcome_models = "none",
#'   constraints = list(list(name = "simplex")),
#'   inference_type = "placebo",
#'   inference_config = list(verbose = FALSE, placebo_cores = 1)
#' )
#' results <- run_spec_curve_analysis(dataset = toy, params = params, cores = 1)
#' names(results)
#' }
run_spec_curve_analysis <- function(dataset, params, cores = 1, output_dir = NULL) {
  # Validate required parameters
  if (is.null(params$outcomes)) {
    stop("params$outcomes must be specified")
  }
  
  # Set up all parameters for spec_curve
  all_params <- params
  all_params$dataset <- dataset
  
  # Handle inference configuration - use params settings or defaults
  if (is.null(all_params$inference_type)) {
    all_params$inference_type <- "placebo"
  }
  
  # Handle expected direction - use params setting or default
  if (is.null(all_params$expected_direction)) {
    all_params$expected_direction <- "negative"
  }
  
  # Set up inference configuration with defaults, applying cores parameter
  default_inference_config <- list(
    bootstrap_n_replications = 1000,
    bootstrap_cores = cores,
    placebo_cores = cores,
    verbose = FALSE
  )
  
  if (is.null(all_params$inference_config)) {
    all_params$inference_config <- default_inference_config
  } else {
    # Merge user config with defaults, applying cores to both bootstrap and placebo if not specified
    user_config <- all_params$inference_config
    if (is.null(user_config$bootstrap_cores)) user_config$bootstrap_cores <- cores
    if (is.null(user_config$placebo_cores)) user_config$placebo_cores <- cores
    all_params$inference_config <- utils::modifyList(default_inference_config, user_config)
  }
  
  actual_outcomes <- all_params$outcomes
  
  # Call spec_curve with all parameters
  results <- do.call(spec_curve, all_params)
  
  # Save results if output_dir is provided
  if (!is.null(output_dir)) {
    save_path <- file.path(output_dir, paste0(params$name_treated_unit, "_", 
                                              paste(actual_outcomes, collapse="_"), "_sc.rdata"))
    save(results, file = save_path)
  }
  
  return(results)
}
