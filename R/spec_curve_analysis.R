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
#' @param inference_type Character. Type of inference to perform: "placebo" (default), "bootstrap", or "all".
#'   - "placebo": Traditional placebo-based inference using control units
#'   - "bootstrap": Bootstrap null hypothesis testing with unit-level resampling  
#'   - "all": Both placebo and bootstrap inference (for comparison)
#' @param inference_config List. Configuration for inference methods. Can include:
#'   \itemize{
#'     \item \code{bootstrap_n_replications} - Number of bootstrap replications (default: 1000)
#'     \item \code{bootstrap_cores} - Cores for bootstrap processing (default: 1)
#'     \item \code{placebo_cores} - Cores for placebo processing (default: 1)
#'     \item \code{verbose} - Verbose inference output (default: FALSE)
#'   }
#'
#' @return Data.table in long format containing specification curve results with inference measures.
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
#' # Run analysis with placebo inference (default)
#' results_placebo <- run_spec_curve_analysis(dataset, params = params)
#' 
#' # Run bootstrap inference only 
#' results_bootstrap <- run_spec_curve_analysis(
#'   dataset, 
#'   params = params,
#'   inference_type = "bootstrap",
#'   inference_config = list(
#'     bootstrap_n_replications = 1000,
#'     bootstrap_cores = 4
#'   )
#' )
#' 
#' # Run both inference methods for comparison
#' results_both <- run_spec_curve_analysis(
#'   dataset,
#'   params = params, 
#'   inference_type = "all",
#'   inference_config = list(
#'     bootstrap_n_replications = 500,  # Fewer for speed when running both
#'     bootstrap_cores = 2,
#'     placebo_cores = 1
#'   )
#' )
#' }
run_spec_curve_analysis <- function(dataset, outcomes = NULL, params, output_dir = NULL,
                                   inference_type = "placebo", inference_config = list()) {
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
  
  # Validate inference_type
  inference_type <- match.arg(inference_type, choices = c("placebo", "bootstrap", "all"))
  
  # Set up inference configuration with defaults
  default_config <- list(
    bootstrap_n_replications = 1000,
    bootstrap_cores = 1,
    placebo_cores = 1,
    verbose = FALSE
  )
  inference_config <- modifyList(default_config, inference_config)
  
  # Configure inference parameters based on type
  if (inference_type == "placebo") {
    all_params$include_bootstrap_inference <- FALSE
    all_params$bootstrap_only <- FALSE
  } else if (inference_type == "bootstrap") {
    all_params$include_bootstrap_inference <- TRUE
    all_params$bootstrap_only <- TRUE
  } else if (inference_type == "all") {
    all_params$include_bootstrap_inference <- TRUE
    all_params$bootstrap_only <- FALSE
  }
  
  # Set bootstrap parameters
  all_params$bootstrap_n_replications <- inference_config$bootstrap_n_replications
  all_params$bootstrap_cores <- inference_config$bootstrap_cores
  
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
#' results and ensures output is in long format. Handles both nested results (traditional)
#' and long format data (new) for backward compatibility.
#'
#' @param results List or data.table. Results from \code{spec_curve} in any supported format.
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
  # Check if input is already in long format
  if (data.table::is.data.table(results)) {
    # Input is already long format, return as-is
    message("Input detected as long format data.table, returning as-is")
    return(results)
  } else if (is.list(results) && "long" %in% names(results)) {
    # Input is a list with long format data (from output_format = "both")
    message("Input detected as list containing long format, extracting long data")
    return(results$long)
  } else if (is.list(results) && length(results) > 0 && all(sapply(results, is.list))) {
    # Input is nested format, process as before
    message("Input detected as nested format, extracting to long format")
  } else {
    stop("Unrecognized input format. Expected nested list, long format data.table, or list with 'long' component.")
  }
  
  sc_results_list <- list()
  
  # Loop through all model specifications and compile results (nested format processing)
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
                    # Make sure they're data.tables (using setDT for efficiency)
                    if (!is.data.table(tau_data)) setDT(tau_data)
                    if (!is.data.table(rmse_data)) setDT(rmse_data)
                    
                    # Merge the datasets on the common keys
                    merged_data <- merge(
                      tau_data,
                      rmse_data,
                      by = c("unit_name", "unit_type", "outcome_model"),
                      all.x = TRUE
                    )
                    
                    # Add Abadie significance measures if available
                    if ("abadie_significance" %in% names(current_result$infer)) {
                      abadie_data <- current_result$infer$abadie_significance
                      
                      # Add p-values
                      if (!is.null(abadie_data$p_values) && nrow(abadie_data$p_values) > 0) {
                        p_vals <- abadie_data$p_values[, .(outcome_model, p_value, treated_ratio, rank, total_units)]
                        merged_data <- merge(merged_data, p_vals, by = "outcome_model", all.x = TRUE)
                      }
                      
                      # Add post/pre ratios
                      if (!is.null(abadie_data$post_pre_ratios) && nrow(abadie_data$post_pre_ratios) > 0) {
                        ratios <- abadie_data$post_pre_ratios[, .(unit_name, outcome_model, unit_type, post_pre_ratio, post_rmspe)]
                        merged_data <- merge(merged_data, ratios, by = c("unit_name", "outcome_model", "unit_type"), all.x = TRUE)
                      }
                      
                      # Add filtered results
                      if (!is.null(abadie_data$filtered_results) && nrow(abadie_data$filtered_results) > 0) {
                        filtered <- abadie_data$filtered_results[, .(outcome_model, p_value_filtered, rank_filtered, total_units_filtered, units_excluded, rmspe_threshold)]
                        merged_data <- merge(merged_data, filtered, by = "outcome_model", all.x = TRUE)
                      }
                    }
                    
                    # Add bootstrap significance measures if available
                    if ("bootstrap_significance" %in% names(current_result$infer)) {
                      bootstrap_data <- current_result$infer$bootstrap_significance
                      
                      # Add bootstrap p-values
                      if (!is.null(bootstrap_data$p_values) && nrow(bootstrap_data$p_values) > 0) {
                        bootstrap_pvals <- bootstrap_data$p_values[, .(outcome_model, 
                                                                      bootstrap_p_value_two_tailed = p_value_two_tailed,
                                                                      bootstrap_p_value_positive = p_value_positive,
                                                                      bootstrap_p_value_negative = p_value_negative,
                                                                      bootstrap_actual_effect = actual_effect,
                                                                      bootstrap_n = n_bootstrap,
                                                                      bootstrap_mean = bootstrap_mean,
                                                                      bootstrap_sd = bootstrap_sd)]
                        merged_data <- merge(merged_data, bootstrap_pvals, by = "outcome_model", all.x = TRUE)
                      }
                      
                      # Add bootstrap effects distribution for plotting null
                      if (!is.null(bootstrap_data$bootstrap_effects) && length(bootstrap_data$bootstrap_effects) > 0) {
                        effects_list <- list()
                        for (oc_model in names(bootstrap_data$bootstrap_effects)) {
                          effects <- bootstrap_data$bootstrap_effects[[oc_model]]
                          if (length(effects) > 0) {
                            effects_list[[oc_model]] <- data.table(
                              outcome_model = oc_model,
                              tau = effects,
                              unit_name = paste0("bootstrap_rep_", seq_along(effects)),
                              unit_type = "bootstrap"
                            )
                          }
                        }
                        if (length(effects_list) > 0) {
                          bootstrap_effects_df <- rbindlist(effects_list, fill = TRUE)
                          # Combine with main data
                          merged_data <- rbindlist(list(merged_data, bootstrap_effects_df), fill = TRUE)
                        }
                      }
                    }
                    
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

#' Extract Abadie Significance Measures from Specification Curve Results
#'
#' @title Extract Abadie-style P-values and RMSPE Ratios 
#' @description Extracts Abadie-style significance measures (post/pre RMSPE ratios,
#' rank-based p-values, and filtering results) from specification curve analysis results.
#' This provides significance testing following Abadie et al. (2010) German reunification methodology.
#' Uses the existing extraction infrastructure for efficiency.
#'
#' @param results List. Results from \code{spec_curve} in nested format.
#' @param include_filtering Logical. Whether to include RMSPE filtering results. Default is TRUE.
#'
#' @return List containing three data tables:
#'   \itemize{
#'     \item \code{p_values} - Rank-based p-values for each specification
#'     \item \code{post_pre_ratios} - Post/pre RMSPE ratios for all units and specifications
#'     \item \code{filtered_results} - Results after RMSPE filtering (optional)
#'   }
#'   Each table includes specification metadata (outcome, constraint, feature weighting, etc.)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract Abadie significance measures
#' abadie_results <- extract_abadie_significance(spec_curve_results)
#' 
#' # View p-values
#' print(abadie_results$p_values)
#' 
#' # View significant specifications (p < 0.05)
#' significant_specs <- abadie_results$p_values[p_value < 0.05]
#' }
extract_abadie_significance <- function(results, include_filtering = TRUE) {
  
  if (!is.list(results) || length(results) == 0) {
    stop("Results must be a non-empty list from spec_curve analysis")
  }
  
  # Use helper function to avoid code duplication
  extract_abadie_component <- function(results, component_name) {
    component_list <- list()
    
    # Reuse the existing nested extraction pattern from extract_spec_curve_results
    for(outcome in names(results)) {
      for(const in names(results[[outcome]])) {
        for(fw in names(results[[outcome]][[const]])) {
          for(feat in names(results[[outcome]][[const]][[fw]])) {
            for(ds in names(results[[outcome]][[const]][[fw]][[feat]])) {
              for(ny in names(results[[outcome]][[const]][[fw]][[feat]][[ds]])) {
                current_result <- results[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]
                
                if (is.null(current_result)) next
                
                # Extract Abadie significance component
                tryCatch({
                  if (!is.null(current_result$infer) && 
                      "abadie_significance" %in% names(current_result$infer) && 
                      !is.null(current_result$infer$abadie_significance[[component_name]]) &&
                      nrow(current_result$infer$abadie_significance[[component_name]]) > 0) {
                    
                    component_data <- copy(current_result$infer$abadie_significance[[component_name]])
                    
                    # Add specification metadata
                    component_data[, `:=`(
                      outcome = outcome,
                      const = const,
                      fw = fw,
                      data_sample = ds,
                      feat = feat,
                      num_pre_period_years = ny
                    )]
                    
                    component_list[[length(component_list) + 1]] <- component_data
                  }
                }, error = function(e) {
                  warning("Error extracting ", component_name, " for spec: ", 
                          outcome, "/", const, "/", fw, "/", feat, "/", ds, "/", ny)
                })
              }
            }
          }
        }
      }
    }
    
    if (length(component_list) > 0) {
      return(rbindlist(component_list, fill = TRUE))
    } else {
      return(data.table())
    }
  }
  
  # Extract each component efficiently
  result <- list(
    p_values = extract_abadie_component(results, "p_values"),
    post_pre_ratios = extract_abadie_component(results, "post_pre_ratios")
  )
  
  if (include_filtering) {
    result$filtered_results <- extract_abadie_component(results, "filtered_results")
  }
  
  # Provide warnings for empty results
  if (nrow(result$p_values) == 0) warning("No p-values found in results")
  if (nrow(result$post_pre_ratios) == 0) warning("No post/pre ratios found in results")
  if (include_filtering && nrow(result$filtered_results) == 0) warning("No filtered results found")
  
  return(result)
}

#' Extract Bootstrap Significance Measures from Specification Curve Results
#'
#' @title Extract Bootstrap-based P-values and Effects Distribution
#' @description Extracts bootstrap-based significance measures from specification curve analysis results.
#' This provides significance testing using the "no effect" bootstrap approach that enforces the null
#' hypothesis and resamples at the unit level to generate a distribution of treatment effects under
#' the null hypothesis.
#'
#' @param results List. Results from \code{spec_curve} in nested format.
#'
#' @return List containing two data tables:
#'   \itemize{
#'     \item \code{p_values} - Bootstrap p-values for each specification with effect sizes and distribution statistics
#'     \item \code{effects_distribution} - Full bootstrap effects distribution for visualization (if available)
#'   }
#'   Each table includes specification metadata (outcome, constraint, feature weighting, etc.)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract bootstrap significance measures
#' bootstrap_results <- extract_bootstrap_significance(spec_curve_results)
#' 
#' # View p-values
#' print(bootstrap_results$p_values)
#' 
#' # View significant specifications (two-tailed p < 0.05)
#' significant_specs <- bootstrap_results$p_values[bootstrap_p_value_two_tailed < 0.05]
#' }
extract_bootstrap_significance <- function(results) {
  
  if (!is.list(results) || length(results) == 0) {
    stop("Results must be a non-empty list from spec_curve analysis")
  }
  
  # Use helper function to avoid code duplication
  extract_bootstrap_component <- function(results, component_name) {
    component_list <- list()
    
    # Reuse the existing nested extraction pattern
    for(outcome in names(results)) {
      for(const in names(results[[outcome]])) {
        for(fw in names(results[[outcome]][[const]])) {
          for(feat in names(results[[outcome]][[const]][[fw]])) {
            for(ds in names(results[[outcome]][[const]][[fw]][[feat]])) {
              for(ny in names(results[[outcome]][[const]][[fw]][[feat]][[ds]])) {
                current_result <- results[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]
                
                if (is.null(current_result)) next
                
                # Extract bootstrap significance component
                tryCatch({
                  if (!is.null(current_result$infer) && 
                      "bootstrap_significance" %in% names(current_result$infer) && 
                      !is.null(current_result$infer$bootstrap_significance[[component_name]])) {
                    
                    component_data <- current_result$infer$bootstrap_significance[[component_name]]
                    
                    if (component_name == "p_values" && nrow(component_data) > 0) {
                      component_data <- copy(component_data)
                    } else if (component_name == "bootstrap_effects" && length(component_data) > 0) {
                      # Convert named list of effect vectors to long format data.table
                      effects_list <- list()
                      for (oc in names(component_data)) {
                        if (length(component_data[[oc]]) > 0) {
                          effects_list[[oc]] <- data.table(
                            outcome_model = oc,
                            bootstrap_effect = component_data[[oc]],
                            bootstrap_id = seq_along(component_data[[oc]])
                          )
                        }
                      }
                      if (length(effects_list) > 0) {
                        component_data <- rbindlist(effects_list)
                      } else {
                        component_data <- data.table()
                      }
                    } else {
                      next  # Skip if no valid data
                    }
                    
                    # Add specification metadata
                    component_data[, `:=`(
                      outcome = outcome,
                      const = const,
                      fw = fw,
                      data_sample = ds,
                      feat = feat,
                      num_pre_period_years = ny
                    )]
                    
                    component_list[[length(component_list) + 1]] <- component_data
                  }
                }, error = function(e) {
                  warning("Error extracting ", component_name, " for spec: ", 
                          outcome, "/", const, "/", fw, "/", feat, "/", ds, "/", ny)
                })
              }
            }
          }
        }
      }
    }
    
    if (length(component_list) > 0) {
      return(rbindlist(component_list, fill = TRUE))
    } else {
      return(data.table())
    }
  }
  
  # Extract each component efficiently
  result <- list(
    p_values = extract_bootstrap_component(results, "p_values"),
    effects_distribution = extract_bootstrap_component(results, "bootstrap_effects")
  )
  
  # Provide warnings for empty results
  if (nrow(result$p_values) == 0) warning("No bootstrap p-values found in results")
  if (nrow(result$effects_distribution) == 0) warning("No bootstrap effects distribution found in results")
  
  return(result)
}