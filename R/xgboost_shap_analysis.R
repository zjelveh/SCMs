# XGBoost functions imported via NAMESPACE (xgb.DMatrix, xgb.train, etc.)

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
#' Create Configuration for Direct Long Format XGBoost Analysis
#'
#' @title Create Configuration for Direct Long Format Analysis
#' @description Creates configuration for XGBoost analysis that works directly
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
create_xgboost_config <- function(dataset_name, 
                                        treated_unit_name,
                                        outcome_filter = NULL,
                                        spec_features = c("outcome_model", "const", "fw", "feat", "data_sample"),
                                        treated_unit_only = TRUE,
                                        xgboost_params = NULL) {
  
  # Set default XGBoost parameters if not provided
  # Tuned via LOO CV grid search: d=10, eta=0.05, n=500, sub=0.8
  # achieves LOO R²=0.57 (full, n=280) and 0.37 (refined, n=72)
  if (is.null(xgboost_params)) {
    xgboost_params <- list(
      objective = "reg:squarederror",
      max_depth = 10,
      eta = 0.05,
      nrounds = 500,
      subsample = 0.8,
      colsample_bytree = 0.8,
      nthread = 1,
      seed = 42,
      verbose = 0
    )
  }
  
  config <- list(
    dataset_name = dataset_name,
    treated_unit_name = treated_unit_name,
    outcome_filter = outcome_filter,
    spec_features = spec_features,
    treated_unit_only = treated_unit_only,
    xgboost_params = xgboost_params
  )
  
  return(config)
}



#' Run XGBoost SHAP Analysis on Long Format Data
#'
#' @title Run Direct XGBoost Analysis on Long Format Data
#' @description Performs XGBoost modeling and SHAP analysis directly on long format
#' data from spec_curve(..., output_format = "long"). This eliminates the need
#' for CSV files and complex data reconstruction.
#'
#' @param long_data Data.table. Long format data from spec_curve analysis.
#' @param config List. Configuration object from \code{create_xgboost_config}.
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
run_xgboost_shap_analysis <- function(long_data, config, compute_loo = TRUE) {
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
  # Include spec features needed for XGBoost analysis
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
      # Mark categorical features for XGBoost
      if (is.character(avg_tau_data[[col]]) || is.factor(avg_tau_data[[col]])) {
        categorical_features <- c(categorical_features, col)
      }
    }
  }
  
  for(unit in unique_units) {
    unit_data <- avg_tau_data[unit_name == unit]
    analysis_results <- analyze_unit_xgboost(unit_data, unit, config, all_factor_levels, categorical_features, compute_loo = compute_loo)
    
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
#' @title XGBoost Analysis for Single Unit (Direct Long Format)
#' @description Performs XGBoost modeling and SHAP analysis for a single unit
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
analyze_unit_xgboost <- function(unit_data, unit, config,
      all_factor_levels = NULL, categorical_features = NULL,
      compute_loo = TRUE) {
  cat(paste("Processing unit:", unit), "\n")

  outcome_col_name <- "avg_tau"
  spec_features <- config$spec_features

  min_specs_required <- 3
  if (nrow(unit_data) < min_specs_required) {
    cat("Skipping", unit, "- insufficient specifications:", nrow(unit_data), "\n")
    return(list(results = NULL, shapley = NULL, predictions = NULL))
  }

  # Check for complete cases
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

  # Prepare data for XGBoost: one-hot encode categorical features
  design_data <- data.table::copy(unit_data_complete)
  feature_groups <- list()

  for (col in features_to_use) {
    if (is.character(design_data[[col]]) || is.factor(design_data[[col]])) {
      if (!is.null(all_factor_levels) && col %in% names(all_factor_levels)) {
        all_levels <- all_factor_levels[[col]]
      } else {
        all_levels <- unique(design_data[[col]])
      }
      design_data[[col]] <- factor(design_data[[col]], levels = all_levels)
      feature_groups[[col]] <- all_levels
      cat("Factor", col, "levels:", paste(all_levels, collapse = ", "), "\n")
    }
  }

  # One-hot encode via model.matrix (no intercept)
  X_df <- as.data.frame(design_data[, ..features_to_use])
  X_mat <- model.matrix(~ . - 1, data = X_df)

  # Build mapping: one-hot column → (feature_group, feature_level)
  onehot_mapping <- data.table::data.table(
    onehot_col = colnames(X_mat),
    feature_group = character(ncol(X_mat)),
    feature_level = character(ncol(X_mat))
  )
  for (feat in features_to_use) {
    mask <- grepl(paste0("^", feat), colnames(X_mat))
    onehot_mapping[mask, feature_group := feat]
    onehot_mapping[mask, feature_level := gsub(paste0("^", feat), "", colnames(X_mat)[mask])]
  }

  cat("One-hot matrix dimensions:", nrow(X_mat), "x", ncol(X_mat), "\n")
  cat("Feature groups:", paste(unique(onehot_mapping$feature_group), collapse = ", "), "\n")

  y_full <- unit_data_complete[[outcome_col_name]]

  # Extract XGBoost training params (nrounds is separate from params list)
  # Note: seed must be set via set.seed() in R, not in params (xgboost R package ignores seed param)
  params <- config$xgboost_params
  nrounds <- params$nrounds
  params_train <- params[setdiff(names(params), c("nrounds", "verbose", "seed"))]
  if (!is.null(params$seed)) set.seed(params$seed)

  if (compute_loo) {
    # LOO test predictions
    all_test_preds <- numeric(n_specs)

    cat("Starting LOO test predictions for", unit, "with", n_specs, "specifications...\n")
    for (i in 1:n_specs) {
      y_train_fold <- y_full[-i]

      if (length(y_train_fold) < 2 || sd(y_train_fold, na.rm = TRUE) == 0) {
        all_test_preds[i] <- NA
        next
      }

      dtrain <- xgboost::xgb.DMatrix(data = X_mat[-i, , drop = FALSE], label = y_train_fold)
      dtest <- xgboost::xgb.DMatrix(data = X_mat[i, , drop = FALSE])

      model_fold <- xgboost::xgb.train(
        params = params_train,
        data = dtrain,
        nrounds = nrounds,
        verbose = 0
      )

      all_test_preds[i] <- xgboost:::predict.xgb.Booster(model_fold, dtest)
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
  } else {
    all_test_preds <- rep(NA_real_, n_specs)
    test_correlation <- NA
    test_rmse <- NA
    cat("Skipping LOO test predictions for", unit, "(compute_loo=FALSE)\n")
  }

  # Train final model on all data for SHAP
  cat("Training final model for SHAP values for unit", unit, "...\n")
  dfull <- xgboost::xgb.DMatrix(data = X_mat, label = y_full)
  xgb_model_final <- xgboost::xgb.train(
    params = params_train,
    data = dfull,
    nrounds = nrounds,
    verbose = 0
  )

  train_preds_final <- xgboost:::predict.xgb.Booster(xgb_model_final, dfull)
  train_correlation_final <- cor(train_preds_final, y_full)
  train_rmse_final <- sqrt(mean((train_preds_final - y_full)^2))

  cat("Train R²:", 1 - sum((y_full - train_preds_final)^2) / sum((y_full - mean(y_full))^2), "\n")
  if (compute_loo) {
    valid_preds_indices <- !is.na(all_test_preds)
    cat("LOO R²:", 1 - sum((y_full[valid_preds_indices] - all_test_preds[valid_preds_indices])^2) /
          sum((y_full[valid_preds_indices] - mean(y_full[valid_preds_indices]))^2), "\n")
  }

  # Calculate SHAP values using XGBoost's built-in TreeSHAP (predcontrib=TRUE)
  cat("Calculating SHAP values for", nrow(X_mat), "observations with", ncol(X_mat), "one-hot features...\n")
  shap_matrix <- xgboost:::predict.xgb.Booster(xgb_model_final, dfull, predcontrib = TRUE)

  # Last column is BIAS term — remove it
  shap_matrix <- shap_matrix[, -ncol(shap_matrix), drop = FALSE]
  colnames(shap_matrix) <- colnames(X_mat)

  cat("SHAP matrix dimensions:", dim(shap_matrix), "\n")

  # Store prediction data
  predictions_dt <- data.table::data.table(
    unit = unit,
    actual = y_full,
    predicted_loo = all_test_preds,
    predicted_train_final = train_preds_final,
    is_treated = ifelse(unit == config$treated_unit_name, config$treated_unit_name, "Other")
  )

  # Add specification metadata
  metadata_cols <- c("full_spec_id", config$spec_features)
  available_cols <- intersect(metadata_cols, names(unit_data_complete))
  if (length(available_cols) > 0) {
    spec_metadata <- unit_data_complete[, ..available_cols]
    predictions_dt <- cbind(predictions_dt, spec_metadata)
  }

  # Aggregate one-hot SHAP values back to original categorical features
  # For each (observation, feature_group): sum all one-hot SHAP values,
  # attribute the sum to the active categorical level
  shap_long <- data.table::data.table()

  for (feat in features_to_use) {
    feat_cols <- onehot_mapping[feature_group == feat, onehot_col]

    # Sum SHAP values across all one-hot columns for this feature group
    if (length(feat_cols) == 1) {
      feat_shap_sum <- shap_matrix[, feat_cols]
    } else {
      feat_shap_sum <- rowSums(shap_matrix[, feat_cols, drop = FALSE])
    }

    # Look up the actual categorical level for each observation
    actual_levels <- as.character(unit_data_complete[[feat]])

    feat_dt <- data.table::data.table(
      unit = unit,
      is_treated = ifelse(unit == config$treated_unit_name, config$treated_unit_name, "Other"),
      full_spec_id = unit_data_complete$full_spec_id,
      feature_group = feat,
      feature = actual_levels,
      shapley_value = feat_shap_sum
    )

    shap_long <- rbind(shap_long, feat_dt)
  }

  cat("SHAP long format:", nrow(shap_long), "rows\n")
  cat("  Unique feature_groups:", paste(unique(shap_long$feature_group), collapse = ", "), "\n")
  cat("  Unique features:", paste(unique(shap_long$feature), collapse = ", "), "\n")

  # Add metadata columns for downstream compatibility
  shap_long[, is_categorical := feature_group %in% names(feature_groups)]
  shap_long[, categorical_levels := ""]
  for (cat_feat in names(feature_groups)) {
    shap_long[feature_group == cat_feat,
              categorical_levels := paste(feature_groups[[cat_feat]], collapse = ",")]
  }

  # Create per-feature-group results summary (aggregate one-hot SHAP)
  shap_summary <- shap_long[, .(
    mean_shap = mean(shapley_value),
    sd_shap = sd(shapley_value),
    mean_abs_shap = mean(abs(shapley_value))
  ), by = .(feature_group)]

  results_dt <- data.table::data.table(
    unit = unit,
    feature = shap_summary$feature_group,
    mean_shap = shap_summary$mean_shap,
    sd_shap = shap_summary$sd_shap,
    mean_abs_shap = shap_summary$mean_abs_shap,
    train_correlation = train_correlation_final,
    test_correlation = test_correlation,
    train_rmse = train_rmse_final,
    test_rmse = test_rmse
  )
  results_dt[, feature_group := feature]
  results_dt[, is_categorical := feature %in% names(feature_groups)]
  results_dt[, categorical_levels := ""]
  for (cat_feat in names(feature_groups)) {
    results_dt[feature == cat_feat,
               categorical_levels := paste(feature_groups[[cat_feat]], collapse = ",")]
  }

  cat("Completed analysis for", unit, "\n")

  return(list(
    results = results_dt,
    shapley = shap_long,
    predictions = predictions_dt,
    model = xgb_model_final,
    feature_matrix_onehot = X_mat,
    onehot_mapping = onehot_mapping
  ))
}
