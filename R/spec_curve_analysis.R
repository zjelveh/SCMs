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
#' \dontrun{
# Create analysis configuration (now includes ALL parameters)
#' params <- list(
#'   outcomes = "gdp",
#'   col_name_unit_name = "country",
#'   name_treated_unit = "West Germany",
#'   covagg = list(
#'     gdp_avg = list(var = "gdp", average = "full_pre")
#'   ),
#'   treated_period = 1990,
#'   min_period = 1975,
#'   end_period = 2003,
#'   col_name_period = "year",
#'   feature_weights = c("uniform"),
#'   donor_sample = c("all", "most_similar"),
#'   outcome_models = c("none", "lasso"),
#'   constraints = list(list(name = "simplex"))
#' )
#'
#' # Run analysis with placebo inference (default)
#' results_placebo <- run_spec_curve_analysis(dataset, params)
#' 
#' # Run analysis expecting positive treatment effects
#' params_positive <- params
#' params_positive$expected_direction <- "positive"
#' results_positive <- run_spec_curve_analysis(dataset, params_positive)
#' 
#' # Run bootstrap inference with 4 cores
#' params_bootstrap <- params
#' params_bootstrap$inference_type <- "bootstrap"
#' params_bootstrap$inference_config <- list(bootstrap_n_replications = 1000)
#' results_bootstrap <- run_spec_curve_analysis(dataset, params_bootstrap, cores = 4)
#' 
#' # Run both inference methods for comparison  
#' params_both <- params
#' params_both$inference_type <- "all"
#' params_both$inference_config <- list(
#'   bootstrap_n_replications = 500,  # Fewer for speed when running both
#'   verbose = TRUE
#' )
#' results_both <- run_spec_curve_analysis(dataset, params_both, cores = 2)
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
  
  # Remove output_format parameter (no longer used)
  all_params$output_format <- NULL
  
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

#' Extract Results from Specification Curve Analysis
#'
#' @title Extract and Compile Specification Curve Results
#' @description Extracts treatment effects and RMSE data from specification curve
#' results. Only supports the new structured format from spec_curve.
#'
#' @param results List. Results from \code{spec_curve} in structured format.
#'
#' @return Data table containing compiled results with columns:
#'   \itemize{
#'     \item \code{unit_name} - Unit identifier
#'     \item \code{period} - Time period
#'     \item \code{tau} - Treatment effect estimate
#'     \item \code{post_period} - Logical indicating post-treatment period
#'     \item \code{pre_rmse} - Pre-treatment RMSE
#'     \item \code{outcome} - Outcome variable name
#'     \item \code{const} - Constraint specification
#'     \item \code{fw} - Feature weighting method
#'     \item \code{data_sample} - Donor sample method
#'     \item \code{feat} - Feature specification
#'     \item \code{spec_number} - Specification number
#'     \item \code{full_spec_id} - Full specification identifier
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract results from specification curve analysis
#' compiled_results <- extract_spec_curve_results(spec_curve_results)
#' }
extract_spec_curve_results <- function(results) {
  # Check if input is already in long format
  if (data.table::is.data.table(results)) {
    # Input is already long format, return as-is
    return(results)
  } else if (is.list(results) && "results" %in% names(results)) {
    # Input is new structured format from spec_curve (results, abadie_inference, bootstrap_inference)
    return(results$results)
  } else {
    stop("Unrecognized input format. Expected structured format from spec_curve with 'results' component or long format data.table.")
  }
}

#' Save Specification Curve Results
#'
#' @title Save Compiled Results to File
#' @description Saves compiled specification curve results to CSV format.
#'
#' @param results_df Data table containing compiled results from \code{extract_spec_curve_results}.
#' @param base_path Character. Base path for output file (without extension).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Save results to CSV
#' save_spec_curve_results(compiled_results, "path/to/output/results")
#' }
save_spec_curve_results <- function(results_df, base_path) {
  # Clean up results if needed
  if ("unit_type" %in% names(results_df)) {
    results_df[, unit_type := NULL]
  }
  
  # Save to CSV
  fwrite(results_df, paste0(base_path, ".csv"))
}


