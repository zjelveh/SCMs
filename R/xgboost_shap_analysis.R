# Mask the custom predict function to avoid conflicts with XGBoost predict
predict <- stats::predict

#' XGBoost and SHAP Analysis Configuration
#'
#' @title Create Configuration for XGBoost SHAP Analysis
#' @description Creates configuration objects for running XGBoost models with SHAP
#' analysis on specification curve results.
#'
#' @param dataset_name Character. Name of the dataset for identification.
#' @param file_path Character. Path to CSV file containing specification curve results.
#' @param treated_unit_name Character. Name of the treated unit.
#' @param outcome_filter Character. Outcome variable to filter on (optional).
#' @param spec_train_test_split Logical. Whether to split specifications for train/test.
#' @param spec_split_ratio Numeric. Proportion of specifications for training (0-1).
#' @param spec_features Character vector. Specification features to use as predictors.
#'
#' @return List containing configuration parameters for XGBoost SHAP analysis.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create configuration for homicide analysis
#' config <- create_xgboost_config(
#'   dataset_name = "homicide_study",
#'   file_path = "results/spec_curve_results.csv",
#'   treated_unit_name = "PAPEP0000",
#'   outcome_filter = "num_homicide"
#' )
#' }
create_xgboost_config <- function(dataset_name, 
                                  file_path, 
                                  treated_unit_name,
                                  outcome_filter = NULL,
                                  spec_train_test_split = TRUE,
                                  spec_split_ratio = 0.7,
                                  spec_features = c("outcome_model", "const", "fw", "feat")) {
  
  config <- list(
    dataset_name = dataset_name,
    file_path = file_path,
    treated_unit_name = treated_unit_name,
    outcome_filter = outcome_filter,
    spec_train_test_split = spec_train_test_split,
    spec_split_ratio = spec_split_ratio,
    spec_features = spec_features
  )
  
  return(config)
}

#' Load and Prepare Data for XGBoost Analysis
#'
#' @title Load Specification Curve Results for Machine Learning
#' @description Loads specification curve results and prepares them for XGBoost
#' analysis by marking treated units and cleaning data.
#'
#' @param config List. Configuration object from \code{create_xgboost_config}.
#'
#' @return Data table with loaded and cleaned specification curve results.
#'
#' @export
load_xgboost_data <- function(config) {
  # Load the data
  sc_data <- fread(config$file_path)
  
  # Mark the treated unit
  sc_data[, treated_unit := ifelse(unit_name == config$treated_unit_name, TRUE, FALSE)]
  
  # Clean up data
  if("unit_type" %in% names(sc_data)) sc_data[, unit_type := NULL]
  if("num_pre_period_years" %in% names(sc_data)) sc_data[, num_pre_period_years := NULL]
  
  # Ensure we have a year column
  if(!"year" %in% names(sc_data)) {
    sc_data[, year := period]
  }
  
  return(sc_data)
}

#' Prepare Data for XGBoost Analysis
#'
#' @title Prepare and Aggregate Specification Data
#' @description Filters and aggregates specification curve data for XGBoost analysis,
#' computing average treatment effects by unit and specification.
#'
#' @param sc_data Data table. Loaded specification curve results.
#' @param config List. Configuration object from \code{create_xgboost_config}.
#'
#' @return Data table with aggregated treatment effects by unit and specification.
#'
#' @export
prepare_xgboost_data <- function(sc_data, config) {
  # Filter by outcome if specified
  if (!is.null(config$outcome_filter)) {
    # Ensure outcome column exists before filtering
    if (!config$outcome_filter %in% names(sc_data)) {
      warning(paste("Outcome filter column '", config$outcome_filter, "' not found in data. Skipping outcome filter."), call. = FALSE)
      data_filtered <- sc_data[post_period == TRUE] # Apply only post_period filter
    } else {
      data_filtered <- sc_data[outcome == config$outcome_filter & post_period == TRUE]
    }
  } else {
    data_filtered <- sc_data[post_period == TRUE]
  }
  
  # Check if tau column exists
  if (!"tau" %in% names(data_filtered)) {
    stop("Column 'tau' not found in the filtered data. Cannot calculate average tau.", call. = FALSE)
  }
  
  # Check if spec_features exist
  missing_spec_features <- setdiff(config$spec_features, names(data_filtered))
  if (length(missing_spec_features) > 0) {
    stop(paste("Missing specification features:", paste(missing_spec_features, collapse=", ")), call. = FALSE)
  }
  
  # Group by unit and specification, calculate average tau
  grouping_vars <- c("unit_name", config$spec_features)
  
  # Ensure grouping variables exist
  missing_group_vars <- setdiff(grouping_vars, names(data_filtered))
  if (length(missing_group_vars) > 0) {
    stop(paste("Missing grouping variables for aggregation:", paste(missing_group_vars, collapse=", ")), call. = FALSE)
  }
  
  # Perform aggregation
  avg_tau_data <- data_filtered[, .(avg_tau = mean(tau, na.rm = TRUE)), by = grouping_vars]
  
  # Check for cases where avg_tau might be NaN (if all taus were NA within a group)
  avg_tau_data <- avg_tau_data[!is.nan(avg_tau)]
  
  if (nrow(avg_tau_data) == 0) {
    warning("No data remaining after filtering and averaging tau.", call. = FALSE)
  }
  
  return(avg_tau_data)
}

#' Analyze Single Unit with XGBoost and SHAP
#'
#' @title Run XGBoost and SHAP Analysis for Individual Unit
#' @description Performs XGBoost modeling with leave-one-out cross-validation
#' and calculates SHAP values to explain specification choices.
#'
#' @param unit_data Data table. Data for a single unit from \code{prepare_xgboost_data}.
#' @param unit Character. Unit identifier.
#' @param config List. Configuration object from \code{create_xgboost_config}.
#' @param treated_unit_name Character. Name of treated unit for comparison.
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{results} - Feature importance metrics
#'     \item \code{shapley} - SHAP values in long format
#'     \item \code{predictions} - Model predictions with actual values
#'   }
#'
#' @export
analyze_unit_xgboost <- function(unit_data, unit, config, treated_unit_name) {
  cat(paste("Processing unit:", unit), "\n")
  
  outcome_col_name <- "avg_tau"
  spec_features <- config$spec_features
  
  # Need at least ~3 specs for LOO to be meaningful (train on 2, test on 1)
  min_specs_required <- 3
  if (nrow(unit_data) < min_specs_required) {
    cat("Skipping", unit, "- insufficient number of specifications for LOO CV:", nrow(unit_data), "\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }
  
  # Prepare full dataset matrices (handle NAs first)
  check_cols <- c(spec_features, outcome_col_name)
  unit_data_complete <- na.omit(unit_data, cols = check_cols)
  n_specs <- nrow(unit_data_complete)
  
  if (n_specs < min_specs_required) {
    cat("Skipping", unit, "- insufficient specs after NA removal:", n_specs, "\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }
  
  # Check for variation in outcome
  if(sd(unit_data_complete[[outcome_col_name]], na.rm = TRUE) == 0) {
    cat("Skipping", unit, "- no variation in avg_tau after NA removal\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }
  
  # Create full design matrix X and outcome vector y ONCE
  formula_str <- paste("~ 0 +", paste(spec_features, collapse = "+"))
  X_full <- tryCatch({
    as.matrix(model.matrix(as.formula(formula_str), data = unit_data_complete))
  }, error = function(e) { cat("Error full matrix unit", unit, ":", e$message, "\n"); NULL })
  
  if (is.null(X_full)) {
    cat("Skipping", unit, "- Failed to create full feature matrix.\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }
  y_full <- unit_data_complete[[outcome_col_name]]
  
  # LOO CV loop for Test Predictions
  all_test_preds <- numeric(n_specs)
  all_test_actuals <- y_full # Actual values are just the full vector
  
  cat("Starting LOO CV for", unit, "with", n_specs, "specifications...\n")
  for (i in 1:n_specs) {
    # Define train/test indices for this fold
    train_indices <- setdiff(1:n_specs, i)
    test_indices <- i
    
    # Extract data for this fold
    X_train_fold <- X_full[train_indices, , drop = FALSE]
    y_train_fold <- y_full[train_indices]
    X_test_fold <- X_full[test_indices, , drop = FALSE] # Matrix with 1 row
    
    # Skip fold if training data is too small (unlikely here but safe)
    if (length(y_train_fold) < 2 || sd(y_train_fold, na.rm = TRUE) == 0) {
      cat("Skipping LOO fold", i, "for unit", unit, "- insufficient train data or no variance\n")
      all_test_preds[i] <- NA # Mark prediction as NA for this fold
      next
    }
    
    # Create DMatrix for training fold
    dtrain_fold <- xgb.DMatrix(data = X_train_fold, label = y_train_fold)
    
    # XGBoost parameters
    params <- list( 
      objective = "reg:squarederror", 
      eval_metric = "rmse", 
      eta = 0.05,
      max_depth = 2, 
      min_child_weight = 2, 
      subsample = 0.7, 
      colsample_bytree = 0.7 
    )
    
    # Train model for this fold
    xgb_model_fold <- xgb.train( 
      params = params, 
      data = dtrain_fold, 
      nrounds = 10, 
      verbose = 0
    )
    
    # Predict on the held-out specification
    all_test_preds[i] <- predict(xgb_model_fold, X_test_fold)
    
    # Optional: Add progress indicator
    if (i %% 50 == 0 || i == n_specs) { cat("...", i, "/", n_specs) }
  }
  cat("\nFinished LOO CV for", unit, "\n")
  
  # Calculate LOO Performance Metrics
  valid_preds_indices <- !is.na(all_test_preds)
  if (sum(valid_preds_indices) < 2) { # Need at least 2 valid predictions for correlation
    cat("Warning: Not enough valid LOO predictions for unit", unit, "to calculate test metrics.\n")
    test_correlation <- NA
    test_rmse <- NA
  } else {
    test_correlation <- cor(all_test_preds[valid_preds_indices], all_test_actuals[valid_preds_indices])
    test_rmse <- sqrt(mean((all_test_preds[valid_preds_indices] - all_test_actuals[valid_preds_indices])^2))
  }
  
  # Train Final Model on ALL data for SHAP and Train Metrics
  cat("Training final model for SHAP values for unit", unit, "...\n")
  dtrain_full <- xgb.DMatrix(data = X_full, label = y_full)
  xgb_model_final <- xgb.train( 
    params = params, 
    data = dtrain_full, 
    nrounds = 10, 
    verbose = 0
  )
  
  # Get predictions from final model on full training set
  train_preds_final <- predict(xgb_model_final, X_full)
  train_correlation_final <- cor(train_preds_final, y_full)
  train_rmse_final <- sqrt(mean((train_preds_final - y_full)^2))
  
  # Calculate Shapley values using the FINAL model
  shap_values <- tryCatch({
    explain(
      xgb_model_final, # Use final model
      newdata = X_full,    # Use full data matrix
      pred_wrapper = function(model, newdata) predict(model, newdata),
      nsim = 10, # Adjust as needed
      exact = TRUE # Consider FALSE if too slow
    )
  }, error = function(e) { cat("Error SHAP final model unit", unit, ":", e$message, "\n"); NULL })
  
  # Store prediction data (using LOO predictions for 'Test')
  predictions_dt <- data.table(
    unit = unit,
    actual = y_full, # The actual avg_tau values
    predicted_loo = all_test_preds, # LOO predictions (can have NAs)
    predicted_train_final = train_preds_final, # Predictions from model trained on all data
    is_treated = ifelse(unit == treated_unit_name, treated_unit_name, "Other")
  )
  # Add specification features for identification
  spec_features_dt <- unit_data_complete[, ..spec_features] # Select spec feature columns
  predictions_dt <- cbind(predictions_dt, spec_features_dt)
  
  # Process SHAP and Final Results
  shap_long <- data.table()
  results_dt <- data.table()
  
  if (!is.null(shap_values) && inherits(shap_values, "matrix")) {
    shap_dt <- as.data.table(shap_values)
    shap_dt[, unit := unit]
    shap_dt[, is_treated := ifelse(unit == treated_unit_name, treated_unit_name, "Other")]
    shap_long <- melt(shap_dt, id.vars = c("unit", "is_treated"), variable.name = "feature", value.name = "shapley_value")
    
    results_dt <- data.table(
      unit = unit, 
      feature = colnames(X_full), 
      mean_shap = colMeans(shap_values), 
      sd_shap = apply(shap_values, 2, sd), 
      mean_abs_shap = colMeans(abs(shap_values)),
      # Report train metrics from final model, test metrics from LOO
      train_correlation = train_correlation_final, 
      test_correlation = test_correlation,
      train_rmse = train_rmse_final, 
      test_rmse = test_rmse
    )
  } else {
    cat("SHAP values not calculated or invalid for unit", unit, ". Skipping SHAP results.\n")
    results_dt <- data.table( # Minimal results table
      unit = unit, 
      feature = character(0), 
      mean_shap = numeric(0), 
      sd_shap = numeric(0), 
      mean_abs_shap = numeric(0),
      train_correlation = train_correlation_final, 
      test_correlation = test_correlation,
      train_rmse = train_rmse_final, 
      test_rmse = test_rmse
    )
  }
  
  cat("Completed analysis for", unit,
      "- final train correlation:", round(train_correlation_final, 4),
      "- LOO test correlation:", round(test_correlation, 4), "\n")
  
  return(list(
    results = results_dt,
    shapley = shap_long,
    predictions = predictions_dt
  ))
}

#' Run Complete XGBoost SHAP Analysis
#'
#' @title Execute Complete XGBoost and SHAP Analysis Pipeline
#' @description Runs the complete pipeline for XGBoost modeling and SHAP analysis
#' on specification curve results for all units in the dataset.
#'
#' @param config List. Configuration object from \code{create_xgboost_config}.
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{results} - Feature importance results for all units
#'     \item \code{shapley} - SHAP values for all units
#'     \item \code{predictions} - Predictions for all units
#'     \item \code{config} - Original configuration object
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create configuration and run analysis
#' config <- create_xgboost_config(
#'   dataset_name = "homicide_study",
#'   file_path = "results/spec_curve_results.csv", 
#'   treated_unit_name = "PAPEP0000"
#' )
#' 
#' results <- run_xgboost_shap_analysis(config)
#' }
run_xgboost_shap_analysis <- function(config) {
  cat(paste("\nProcessing dataset:", config$dataset_name), "\n")
  
  # Load and prepare data
  sc_data <- load_xgboost_data(config)
  avg_tau_by_spec <- prepare_xgboost_data(sc_data, config)
  
  # Create data tables to store results
  results_dt <- data.table()
  all_shapley_values <- data.table()
  predictions_dt <- data.table()

  # Get unique units
  unique_units <- sort(unique(avg_tau_by_spec$unit_name))
  
  # Process each unit
  for(unit in unique_units) {
    # Filter data for current unit
    unit_data <- avg_tau_by_spec[unit_name == unit]
    # Analyze unit
    analysis_results <- analyze_unit_xgboost(unit_data, unit, config, config$treated_unit_name)

    # Append results if not NULL
    if(!is.null(analysis_results$results)) {
      results_dt <- rbind(results_dt, analysis_results$results)
      all_shapley_values <- rbind(all_shapley_values, analysis_results$shapley)
      predictions_dt <- rbind(predictions_dt, analysis_results$predictions)
    }
  }
  
  # Calculate feature importance rank within each unit
  if(nrow(results_dt) > 0) {
    results_dt[, rnk_abs_mean := rank(-mean_abs_shap), by = unit]
  }
  
  # Return all results
  return(list(
    results = results_dt,
    shapley = all_shapley_values,
    predictions = predictions_dt,
    config = config
  ))
}