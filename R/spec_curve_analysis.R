#' Run Comprehensive Specification Curve Analysis
#'
#' @title Run Specification Curve Analysis with Multiple Configurations
#' @description Performs specification curve analysis across multiple modeling choices
#' and configurations, providing a systematic exploration of researcher degrees of freedom.
#'
#' @param dataset Data frame containing the panel data for analysis.
#' @param outcomes Character vector of outcome variable names. If NULL, uses outcomes from params.
#' @param params List containing analysis parameters including:
#'   \itemize{
#'     \item \code{outcomes} - Character vector of outcome variables
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
#'   }
#' @param output_dir Character. Directory to save results. If NULL, results are not saved.
#'
#' @return List containing specification curve results for all parameter combinations.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create analysis configuration
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
#' # Run analysis
#' results <- run_spec_curve_analysis(dataset, params = params)
#' }
run_spec_curve_analysis <- function(dataset, outcomes = NULL, params, output_dir = NULL) {
  # Check if outcomes is provided as an argument or in params
  if (is.null(outcomes)) {
    # Use outcomes from params if not explicitly provided
    all_params <- params
    actual_outcomes <- params$outcomes
  } else {
    # If outcomes is explicitly provided, override the one in params
    params_copy <- params
    params_copy$outcomes <- NULL  # Remove outcomes from params to avoid duplication
    
    all_params <- c(
      list(
        dataset = dataset,
        outcomes = outcomes
      ),
      params_copy
    )
    actual_outcomes <- outcomes
  }
  
  # Add dataset if not already in the parameters
  if (!"dataset" %in% names(all_params)) {
    all_params$dataset <- dataset
  }
  
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
#' results and compiles them into a single data frame for further analysis.
#'
#' @param results List containing specification curve results from \code{run_spec_curve_analysis}.
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
#'     \item \code{num_pre_period_years} - Number of pre-treatment periods
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
  sc_results_list <- list()
  
  # Loop through all model specifications and compile results
  for(outcome in names(results)) {
    for(const in names(results[[outcome]])) {
      for(fw in names(results[[outcome]][[const]])) {
        for(feat in names(results[[outcome]][[const]][[fw]])) {
          for(ds in names(results[[outcome]][[const]][[fw]][[feat]])) {
            for(ny in names(results[[outcome]][[const]][[fw]][[feat]][[ds]])) {
              # Get the current result
              current_result <- results[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]
              
              # Skip if NULL
              if (is.null(current_result)) next
              
              # Process inference results
              tryCatch({
                # Check if infer component exists and has the right structure
                if (!is.null(current_result$infer) && 
                    is.list(current_result$infer) && 
                    "inference.results" %in% names(current_result$infer) &&
                    "rmse" %in% names(current_result$infer)) {
                  
                  # Get treatment effects data
                  tau_data <- current_result$infer$inference.results
                  
                  # Get RMSE data
                  rmse_data <- current_result$infer$rmse
                  
                  # Check if both datasets exist
                  if (!is.null(tau_data) && !is.null(rmse_data)) {
                    # Make sure they're data.tables
                    if (!is.data.table(tau_data)) tau_data <- as.data.table(tau_data)
                    if (!is.data.table(rmse_data)) rmse_data <- as.data.table(rmse_data)
                    
                    # Merge the datasets on the common keys
                    merged_data <- merge(
                      tau_data,
                      rmse_data,
                      by = c("unit_name", "unit_type", "outcome_model"),
                      all.x = TRUE
                    )
                    
                    # Add metadata
                    merged_data$outcome <- outcome
                    merged_data$const <- const
                    merged_data$fw <- fw
                    merged_data$data_sample <- ds
                    merged_data$feat <- feat
                    merged_data$num_pre_period_years <- ny
                    
                    # Add to results list
                    sc_results_list[[length(sc_results_list) + 1]] <- merged_data
                  }
                }
              }, error = function(e) {
                warning("Error processing result: ", e$message, 
                        " in spec: ", outcome, "/", const, "/", fw, "/", feat, "/", ds, "/", ny)
              })
            }
          }
        }
      }
    }
  }
  
  # Combine all results into a single dataframe
  if (length(sc_results_list) > 0) {
    sc_results_df <- rbindlist(sc_results_list, fill = TRUE)
    return(sc_results_df)
  } else {
    warning("No valid results found to extract")
    return(NULL)
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