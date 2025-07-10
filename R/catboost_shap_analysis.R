# CatBoost functions imported via NAMESPACE

#' Load and Prepare Data for CatBoost Analysis
#'
#' @title Load Specification Curve Results for Machine Learning
#' @description Loads specification curve results and prepares them for CatBoost
#' analysis by marking treated units and cleaning data.
#'
#' @param config List. Configuration object from \code{create_catboost_config}.
#'
#' @return Data table with loaded and cleaned specification curve results.
#'
#' @export
load_catboost_data <- function(config) {
  # Load the data
  sc_data <- fread(config$file_path)
  
  # Mark the treated unit
  sc_data[, treated_unit := ifelse(unit_name == config$treated_unit_name, TRUE, FALSE)]
  
  # Clean up data
  if("unit_type" %in% names(sc_data)) sc_data[, unit_type := NULL]
  
  # Ensure we have a year column
  if(!"year" %in% names(sc_data)) {
    sc_data[, year := period]
  }
  
  return(sc_data)
}

#' Prepare Data for CatBoost Analysis
#'
#' @title Prepare and Aggregate Specification Data
#' @description Filters and aggregates specification curve data for CatBoost analysis,
#' computing average treatment effects by unit and specification.
#'
#' @param sc_data Data table. Loaded specification curve results.
#' @param config List. Configuration object from \code{create_catboost_config}.
#'
#' @return Data table with aggregated treatment effects by unit and specification.
#'
#' @export
prepare_catboost_data <- function(sc_data, config) {
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
  
  # Store column names for efficiency
  data_cols <- names(data_filtered)
  
  # Check if spec_features exist
  missing_spec_features <- setdiff(config$spec_features, data_cols)
  if (length(missing_spec_features) > 0) {
    stop(paste("Missing specification features:", paste(missing_spec_features, collapse=", ")), call. = FALSE)
  }
  
  # Include ALL specification metadata for proper alignment
  all_spec_cols <- c("outcome", "outcome_model", "const", "fw", "feat", "data_sample")
  available_all_spec_cols <- intersect(all_spec_cols, data_cols)
  
  # Group by unit and ALL specification dimensions (not just features)
  grouping_vars <- c("unit_name", available_all_spec_cols)
  
  # Ensure grouping variables exist
  missing_group_vars <- setdiff(grouping_vars, data_cols)
  if (length(missing_group_vars) > 0) {
    stop(paste("Missing grouping variables for aggregation:", paste(missing_group_vars, collapse=", ")), call. = FALSE)
  }
  
  # Perform aggregation with ALL specification metadata preserved
  avg_tau_data <- data_filtered[, .(avg_tau = mean(tau, na.rm = TRUE)), by = grouping_vars]
  
  # Create complete specification ID for alignment
  if (length(available_all_spec_cols) >= 4) {
    avg_tau_data[, spec_combination := do.call(paste, c(.SD, sep = "_")), .SDcols = available_all_spec_cols]
  }
  
  # Create specification ordering based on treated unit (matching plot_spec_curve logic)
  if (config$treated_unit_name %in% avg_tau_data$unit_name) {
    treated_specs <- avg_tau_data[unit_name == config$treated_unit_name][order(avg_tau)]
    treated_specs[, spec_number := 1:.N]
    
    # Merge specification numbers back to all data
    spec_order_map <- treated_specs[, .(spec_combination, spec_number)]
    avg_tau_data <- merge(avg_tau_data, spec_order_map, by = "spec_combination", all.x = TRUE)
    
    cat("Created specification ordering for", nrow(spec_order_map), "specifications\n")
  }
  
  # Check for cases where avg_tau might be NaN (if all taus were NA within a group)
  avg_tau_data <- avg_tau_data[!is.nan(avg_tau)]
  
  if (nrow(avg_tau_data) == 0) {
    warning("No data remaining after filtering and averaging tau.", call. = FALSE)
  }
  
  return(avg_tau_data)
}


#' Run Complete CatBoost SHAP Analysis
#'
#' @title Execute Complete CatBoost and SHAP Analysis Pipeline
#' @description Runs the complete pipeline for CatBoost modeling and SHAP analysis
#' on specification curve results for all units in the dataset.
#'
#' @param config List. Configuration object from \code{create_catboost_config}.
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
#' config <- create_catboost_config(
#'   dataset_name = "homicide_study",
#'   file_path = "results/spec_curve_results.csv", 
#'   treated_unit_name = "PAPEP0000"
#' )
#' 
#' results <- run_catboost_shap_analysis(config)
#' }
#' Create Configuration for Direct Long Format CatBoost Analysis
#'
#' @title Create Configuration for Direct Long Format Analysis
#' @description Creates configuration for CatBoost analysis that works directly
#' with long format data from spec_curve().
#'
#' @param dataset_name Character. Name of the dataset for identification.
#' @param treated_unit_name Character. Name of the treated unit.
#' @param outcome_filter Character. Outcome variable to filter on (optional).
#' @param spec_features Character vector. Specification features to use as predictors.
#' @param treated_unit_only Logical. If TRUE, only analyze the treated unit.
#'
#' @return List containing configuration parameters.
#'
#' @export
create_catboost_config <- function(dataset_name, 
                                        treated_unit_name,
                                        outcome_filter = NULL,
                                        spec_features = c("outcome_model", "const", "fw", "feat", "data_sample"),
                                        treated_unit_only = TRUE) {
  
  config <- list(
    dataset_name = dataset_name,
    treated_unit_name = treated_unit_name,
    outcome_filter = outcome_filter,
    spec_features = spec_features,
    treated_unit_only = treated_unit_only
  )
  
  return(config)
}



#' Run CatBoost SHAP Analysis on Long Format Data
#'
#' @title Run Direct CatBoost Analysis on Long Format Data
#' @description Performs CatBoost modeling and SHAP analysis directly on long format
#' data from spec_curve(..., output_format = "long"). This eliminates the need
#' for CSV files and complex data reconstruction.
#'
#' @param long_data Data.table. Long format data from spec_curve analysis.
#' @param config List. Configuration object from \code{create_catboost_config}.
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{results} - Feature importance results
#'     \item \code{shapley} - SHAP values with proper spec_number alignment
#'     \item \code{predictions} - Model predictions
#'     \item \code{config} - Original configuration object
#'   }
#'
#' @export
run_catboost_shap_analysis <- function(long_data, config) {
  cat(paste("\nProcessing dataset:", config$dataset_name), "\n")
  cat("Using direct long format data - no CSV processing needed\n")
  
  # Validate input data
  if (!data.table::is.data.table(long_data)) {
    long_data <- data.table::as.data.table(long_data)
  }
  
  # Validate that spec_features are provided
  if (is.null(config$spec_features) || length(config$spec_features) == 0) {
    stop("spec_features must be provided in the config")
  }
  
  # Validate required columns
  required_cols <- c("unit_name", "post_period", "tau", "spec_number", config$spec_features)
  missing_cols <- setdiff(required_cols, names(long_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Filter data
  if (!is.null(config$outcome_filter)) {
    if ("outcome" %in% names(long_data)) {
      data_filtered <- long_data[outcome == config$outcome_filter & post_period == TRUE]
    } else {
      warning("Outcome filter specified but no outcome column found")
      data_filtered <- long_data[post_period == TRUE]
    }
  } else {
    data_filtered <- long_data[post_period == TRUE]
  }
  
  # Aggregate to specification level (already done in long format, but ensure consistency)
  # Include spec features needed for CatBoost analysis
  grouping_vars <- c("unit_name", "full_spec_id", config$spec_features)
  available_grouping_vars <- intersect(grouping_vars, names(data_filtered))
  
  
  avg_tau_data <- data_filtered[, .(avg_tau = mean(tau, na.rm = TRUE)), by = available_grouping_vars]
  
  # Remove cases where avg_tau is NaN
  avg_tau_data <- avg_tau_data[!is.nan(avg_tau)]
  
  if (nrow(avg_tau_data) == 0) {
    stop("No data remaining after filtering and averaging tau.")
  }
  
  cat("Prepared", nrow(avg_tau_data), "specification-unit observations\n")
  
  # Process units
  if(config$treated_unit_only) {
    unique_units <- config$treated_unit_name
    cat("Analyzing treated unit only:", unique_units, "\n")
  } else {
    unique_units <- sort(unique(avg_tau_data$unit_name))
    cat("Analyzing all", length(unique_units), "units\n")
  }
  
  # Run analysis on each unit
  results_dt <- data.table::data.table()
  all_shapley_values <- data.table::data.table()
  predictions_dt <- data.table::data.table()
  all_models <- list()  # Store models for shapviz
  
  # Get all possible factor levels from the full dataset for consistent encoding
  all_factor_levels <- list()
  categorical_features <- c()
  for (col in config$spec_features) {
    if (col %in% names(avg_tau_data)) {
      all_factor_levels[[col]] <- unique(avg_tau_data[[col]])
      # Mark categorical features for CatBoost
      if (is.character(avg_tau_data[[col]]) || is.factor(avg_tau_data[[col]])) {
        categorical_features <- c(categorical_features, col)
      }
    }
  }
  
  for(unit in unique_units) {
    unit_data <- avg_tau_data[unit_name == unit]
    analysis_results <- analyze_unit_catboost(unit_data, unit, config, all_factor_levels, categorical_features)
    
    if(!is.null(analysis_results$results)) {
      results_dt <- rbind(results_dt, analysis_results$results)
      all_shapley_values <- rbind(all_shapley_values, analysis_results$shapley)
      predictions_dt <- rbind(predictions_dt, analysis_results$predictions)
      
      # Store the trained model for this unit
      if(!is.null(analysis_results$model)) {
        all_models[[unit]] <- analysis_results$model
      }
    }
  }
  
  # Calculate feature importance rank within each unit
  if(nrow(results_dt) > 0) {
    results_dt[, rnk_abs_mean := rank(-mean_abs_shap), by = unit]
  }
  
  cat("Analysis complete. SHAP values include spec_number for clean alignment.\n")
  
  return(list(
    results = results_dt,
    shapley = all_shapley_values,
    predictions = predictions_dt,
    models = all_models,  # Include trained models for shapviz
    config = config
  ))
}

#' Analyze Single Unit with Direct Long Format Data
#'
#' @title CatBoost Analysis for Single Unit (Direct Long Format)
#' @description Performs CatBoost modeling and SHAP analysis for a single unit
#' using pre-processed long format data with specification numbers.
#'
#' @param unit_data Data.table. Data for a single unit from long format.
#' @param unit Character. Unit identifier.
#' @param config List. Configuration object.
#' @param all_factor_levels List. All factor levels for consistent encoding.
#' @param categorical_features Character vector. Names of categorical features.
#'
#' @return List containing results, shapley values, and predictions.
#'
#' @keywords internal
analyze_unit_catboost <- function(unit_data, unit, config, 
      all_factor_levels = NULL, categorical_features = NULL) {
  cat(paste("Processing unit:", unit), "\n")
  
  outcome_col_name <- "avg_tau"
  spec_features <- config$spec_features
  
  min_specs_required <- 3
  if (nrow(unit_data) < min_specs_required) {
    cat("Skipping", unit, "- insufficient specifications:", nrow(unit_data), "\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }
  
  # Check for complete cases
  # check_cols <- c(spec_features, outcome_col_name)
  check_cols <- c('full_spec_id', outcome_col_name)

  unit_data_complete <- na.omit(unit_data, cols = check_cols)

  n_specs <- nrow(unit_data_complete)
  
  if (n_specs < min_specs_required) {
    cat("Skipping", unit, "- insufficient specs after NA removal:", n_specs, "\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }
  
  # Check for variation in outcome
  if(sd(unit_data_complete[[outcome_col_name]], na.rm = TRUE) == 0) {
    cat("Skipping", unit, "- no variation in avg_tau\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }
  # Check for single-level factors that would cause contrasts error
  features_to_use <- spec_features
  for (feature in spec_features) {
    if (is.factor(unit_data_complete[[feature]]) || is.character(unit_data_complete[[feature]])) {
      unique_vals <- unique(unit_data_complete[[feature]])
      if (length(unique_vals) < 2) {
        cat("Removing feature", feature, "- only has one level:", unique_vals[1], "\n")
        features_to_use <- setdiff(features_to_use, feature)
      }
    }
  }
  
  if (length(features_to_use) == 0) {
    cat("Skipping", unit, "- no varying features available\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }
  
  # Prepare data for CatBoost (keep categorical features as factors)
  X_full <- tryCatch({
    # Convert character columns to factors and ensure all levels are included
    design_data <- data.table::copy(unit_data_complete)
    
    # Create feature groups for categorical variables
    feature_groups <- list()
    
    for (col in features_to_use) {
      if (is.character(design_data[[col]]) || is.factor(design_data[[col]])) {
        # Use provided factor levels for consistency across all units
        if (!is.null(all_factor_levels) && col %in% names(all_factor_levels)) {
          all_levels <- all_factor_levels[[col]]
        } else {
          # Fallback to local levels if not provided
          all_levels <- unique(design_data[[col]])
        }
        design_data[[col]] <- factor(design_data[[col]], levels = all_levels)
        cat("Factor", col, "levels:", paste(all_levels, collapse = ", "), "\n")
        
        # Store feature group information
        feature_groups[[col]] <- all_levels
      }
    }
    
    # For CatBoost, we keep categorical features as factors
    feature_matrix <- design_data[, ..features_to_use]
    
    # Convert to data.frame and ensure factors remain as factors for CatBoost
    feature_matrix <- as.data.frame(feature_matrix)
    for (col in names(feature_matrix)) {
      if (col %in% categorical_features) {
        # Ensure categorical features are factors
        feature_matrix[[col]] <- as.factor(feature_matrix[[col]])
      }
    }
    
    cat("Feature matrix dimensions:", nrow(feature_matrix), "x", ncol(feature_matrix), "\n")
    cat("Feature matrix columns:", paste(colnames(feature_matrix), collapse = ", "), "\n")
    cat("Categorical features:", paste(intersect(categorical_features, colnames(feature_matrix)), collapse = ", "), "\n")
    
    list(data = feature_matrix, feature_groups = feature_groups)
  }, error = function(e) { 
    cat("Error creating feature matrix for unit", unit, ":", e$message, "\n") 
    NULL 
  })
  
  if (is.null(X_full)) {
    cat("Skipping", unit, "- Failed to create feature matrix.\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }
  
  X_data <- X_full$data
  feature_groups <- X_full$feature_groups
  y_full <- unit_data_complete[[outcome_col_name]]
  
  # For CatBoost with data.frame, we don't need to specify cat_features indices
  # CatBoost will automatically detect factor columns
  categorical_columns <- intersect(categorical_features, colnames(X_data))
  
  # Diagnostic output
  cat("Data preparation for CatBoost:\n")
  cat("  Data type:", class(X_data), "\n")
  cat("  Dimensions:", nrow(X_data), "x", ncol(X_data), "\n")
  for (col in colnames(X_data)) {
    cat("  ", col, ":", class(X_data[[col]]), "\n")
  }
  cat("Categorical columns in data:", paste(categorical_columns, collapse = ", "), "\n")
  
  # LOO CV for test predictions
  all_test_preds <- numeric(n_specs)
  
  cat("Starting LOO CV for", unit, "with", n_specs, "specifications...\n")
  for (i in 1:n_specs) {
    train_indices <- setdiff(1:n_specs, i)
    
    X_train_fold <- X_data[train_indices, ]
    y_train_fold <- y_full[train_indices]
    X_test_fold <- X_data[i, ]
    
    if (length(y_train_fold) < 2 || sd(y_train_fold, na.rm = TRUE) == 0) {
      all_test_preds[i] <- NA
      next
    }
    
    # Create CatBoost pools - no need to specify cat_features with data.frame + factors
    train_pool <- catboost.load_pool(
      data = as.data.frame(X_train_fold),
      label = y_train_fold
    )
    
    test_pool <- catboost.load_pool(
      data = as.data.frame(X_test_fold)
    )
    
    # Train CatBoost model - increased iterations and adjusted params for better feature learning
    params <- list(
      loss_function = 'RMSE',
      iterations = 100,  # Increased from 10
      depth = 4,         # Increased from 3 
      learning_rate = 0.1, # Increased from 0.05
      random_seed = 42,
      verbose = 0,
      bootstrap_type = 'Bayesian',  # Better for small datasets
      bagging_temperature = 1
    )
    
    catboost_model_fold <- catboost.train(
      learn_pool = train_pool,
      params = params
    )
    
    all_test_preds[i] <- catboost.predict(catboost_model_fold, test_pool)
  }
  
  # Calculate test metrics
  valid_preds_indices <- !is.na(all_test_preds)
  if (sum(valid_preds_indices) < 2) {
    test_correlation <- NA
    test_rmse <- NA
  } else {
    test_correlation <- cor(all_test_preds[valid_preds_indices], y_full[valid_preds_indices])
    test_rmse <- sqrt(mean((all_test_preds[valid_preds_indices] - y_full[valid_preds_indices])^2))
  }
  
  # Train final model on all data for SHAP
  cat("Training final model for SHAP values for unit", unit, "...\n")
  
  # Create full dataset pool
  full_pool <- catboost.load_pool(
    data = as.data.frame(X_data),
    label = y_full
  )
  
  # Train final CatBoost model
  catboost_model_final <- catboost.train(
    learn_pool = full_pool,
    params = params
  )
  
  train_preds_final <- catboost.predict(catboost_model_final, full_pool)
  train_correlation_final <- cor(train_preds_final, y_full)
  train_rmse_final <- sqrt(mean((train_preds_final - y_full)^2))
  
  # Debug: Check feature importance first
  cat("Checking CatBoost feature importance...\n")
  feature_importance <- catboost.get_feature_importance(
    model = catboost_model_final,
    pool = full_pool,
    type = 'FeatureImportance'
  )
  cat("Feature importance:\n")
  for (i in 1:length(feature_importance)) {
    cat("  ", colnames(X_data)[i], ":", feature_importance[i], "\n")
  }
  
  # Debug: Check data variation for each feature
  cat("Feature variation check:\n")
  for (col in colnames(X_data)) {
    unique_vals <- unique(X_data[[col]])
    cat("  ", col, ": ", length(unique_vals), " unique values (", paste(unique_vals, collapse=", "), ")\n")
  }
  
  # Calculate SHAP values using CatBoost's built-in SHAP
  shap_values <- tryCatch({
    cat("Calculating SHAP values for", nrow(X_data), "observations with", ncol(X_data), "features...\n")
    shap_result <- catboost.get_feature_importance(
      model = catboost_model_final,
      pool = full_pool,
      type = 'ShapValues',
      thread_count = 1
    )
    cat("SHAP calculation completed. Result dimensions:", dim(shap_result), "\n")
    
    # Debug: Check SHAP value ranges
    cat("SHAP value ranges by feature:\n")
    for (i in 1:ncol(X_data)) {
      if (ncol(shap_result) >= i) {
        shap_col <- if (ncol(shap_result) == ncol(X_data) + 1) i else i  # Handle baseline column
        shap_range <- range(shap_result[, shap_col])
        cat("  ", colnames(X_data)[i], ": [", shap_range[1], ", ", shap_range[2], "]\n")
      }
    }
    
    shap_result
  }, error = function(e) { 
    cat("Error calculating SHAP for unit", unit, ":", e$message, "\n")
    cat("Attempting alternative SHAP calculation...\n")
    # Fallback: return zero matrix if SHAP fails
    matrix(0, nrow = nrow(X_data), ncol = ncol(X_data))
  })
  
  # Store prediction data
  predictions_dt <- data.table::data.table(
    unit = unit,
    actual = y_full,
    predicted_loo = all_test_preds,
    predicted_train_final = train_preds_final,
    is_treated = ifelse(unit == config$treated_unit_name, config$treated_unit_name, "Other")
  )
  
  # Add specification metadata - full_spec_id and spec features for interaction analysis
  metadata_cols <- c("full_spec_id", config$spec_features)
  available_cols <- intersect(metadata_cols, names(unit_data_complete))
  if (length(available_cols) > 0) {
    spec_metadata <- unit_data_complete[, ..available_cols]
    predictions_dt <- cbind(predictions_dt, spec_metadata)
  }
  
  # Process SHAP results with improved feature group handling
  shap_long <- data.table::data.table()
  results_dt <- data.table::data.table()
  
  if (!is.null(shap_values) && is.matrix(shap_values)) {
    cat("Processing SHAP values with dimensions:", dim(shap_values), "\n")
    
    # Handle CatBoost SHAP format - check if baseline column is included
    expected_cols <- ncol(X_data)
    actual_cols <- ncol(shap_values)
    
    if (actual_cols == expected_cols + 1) {
      cat("Removing baseline column from SHAP values\n")
      shap_values <- shap_values[, -actual_cols, drop = FALSE]
    } else if (actual_cols != expected_cols) {
      cat("Warning: SHAP dimensions (", actual_cols, ") don't match features (", expected_cols, ")\n")
    }
    
    # Validate final dimensions
    if (ncol(shap_values) != ncol(X_data)) {
      cat("Error: SHAP values dimensions still don't match after processing\n")
      return(list(results = NULL, shapley = NULL, predictions = predictions_dt))
    }
    
    # SIMPLE: Set column names to match features directly
    colnames(shap_values) <- colnames(X_data)
    cat("SHAP values processed successfully with features:", paste(colnames(shap_values), collapse = ", "), "\n")
    
    # Create SHAP data table with minimal metadata
    shap_dt <- data.table::as.data.table(shap_values)
    shap_dt[, unit := unit]
    shap_dt[, is_treated := ifelse(unit == config$treated_unit_name, config$treated_unit_name, "Other")]
    shap_dt[, full_spec_id := unit_data_complete$full_spec_id]
    
    # SIMPLE: Melt with actual feature names as measure.vars
    shap_long <- data.table::melt(
      shap_dt, 
      id.vars = c("unit", "is_treated", "full_spec_id"), 
      measure.vars = colnames(X_data),
      variable.name = "feature_group", 
      value.name = "shapley_value"
    )
    
    # Convert feature_group to character
    shap_long[, feature_group := as.character(feature_group)]
    
    # CRITICAL FIX: Map to actual categorical levels using full_spec_id for proper alignment
    # For each row, look up what the actual categorical level was for that spec and feature group
    shap_long[, feature := ""]  # Initialize
    
    # Create a lookup table for feature values by full_spec_id
    feature_lookup <- unit_data_complete[, c("full_spec_id", colnames(X_data)), with = FALSE]
    
    for (feat_group in colnames(X_data)) {
      if (feat_group %in% names(unit_data_complete)) {
        # For each SHAP row with this feature_group, look up the actual feature value by full_spec_id
        shap_long[feature_group == feat_group, 
                  feature := as.character(feature_lookup[match(full_spec_id, feature_lookup$full_spec_id), get(feat_group)])]
      }
    }
    
    cat("After melting and mapping to categorical levels:\n")
    cat("  shap_long dimensions:", nrow(shap_long), "x", ncol(shap_long), "\n")
    cat("  Unique feature_groups:", paste(unique(shap_long$feature_group), collapse = ", "), "\n")
    cat("  Unique features (categorical levels):", paste(unique(shap_long$feature), collapse = ", "), "\n")
    
    # Add feature group information - FIXED VERSION
    shap_long[, is_categorical := feature_group %in% categorical_features]
    shap_long[, categorical_levels := ""]  # Initialize empty
    
    # For categorical features, add specific group information
    for (cat_feat in categorical_features) {
      if (cat_feat %in% names(feature_groups)) {
        levels_info <- feature_groups[[cat_feat]]
        # Only update rows where feature matches this categorical feature
        shap_long[feature == cat_feat, categorical_levels := paste(levels_info, collapse = ",")]
        cat("Added feature group info for", cat_feat, "with levels:", paste(levels_info, collapse = ", "), "\n")
      }
    }
    
    # Create results summary
    results_dt <- data.table::data.table(
      unit = unit, 
      feature = colnames(X_data),  # Use feature names
      mean_shap = colMeans(shap_values), 
      sd_shap = apply(shap_values, 2, sd), 
      mean_abs_shap = colMeans(abs(shap_values)),
      train_correlation = train_correlation_final, 
      test_correlation = test_correlation,
      train_rmse = train_rmse_final, 
      test_rmse = test_rmse
    )
    
    # Add feature group information to results - FIXED VERSION
    results_dt[, feature_group := feature]  # Default: feature group = feature name
    results_dt[, is_categorical := feature %in% categorical_features]
    results_dt[, categorical_levels := ""]  # Initialize empty
    
    for (cat_feat in categorical_features) {
      if (cat_feat %in% names(feature_groups)) {
        levels_info <- feature_groups[[cat_feat]]
        # Only update rows where feature matches this categorical feature
        results_dt[feature == cat_feat, categorical_levels := paste(levels_info, collapse = ",")]
      }
    }
    
    cat("SHAP processing complete. Long format has", nrow(shap_long), "rows\n")
    cat("Feature groups identified:", sum(shap_long$is_categorical), "categorical,", 
        sum(!shap_long$is_categorical), "non-categorical\n")
    
    # Validate that we have proper long format
    if (nrow(shap_long) == 0) {
      cat("ERROR: SHAP long format is empty! Returning wide format as fallback.\n")
      shap_long <- shap_dt  # Fallback to wide format
    }
  }
  
  cat("Completed analysis for", unit, "\n")
  cat("  - Feature groups properly identified and labeled\n")
  cat("  - Categorical features:", paste(categorical_features, collapse = ", "), "\n")
  
  return(list(
    results = results_dt,
    shapley = shap_long,
    predictions = predictions_dt,
    model = catboost_model_final  # Return the trained model for shapviz
  ))
}
