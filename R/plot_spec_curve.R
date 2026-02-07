#' Plot Specification Curve Using Direct Long Format Data
#'
#' @title Plot Specification Curve with Internal SHAP Computation and Filtering
#' @description Creates specification curve plots directly from long format data with integrated
#' SHAP analysis and flexible filtering capabilities. Key features:
#' - Internal SHAP computation using XGBoost (automatic when show_shap=TRUE)
#' - Specification filtering BEFORE SHAP computation and p-value calculation
#' - Perfect alignment between SHAP values and specifications via full_spec_id
#' - FAIL HARD error handling with specific guidance when requirements aren't met
#'
#' @param long_data Data.table or List. Long format data from spec_curve().
#'   Can be either the results data.table directly, or the full structured results
#'   (with results, abadie_inference, bootstrap_inference components).
#' @param name_treated_unit Character. Name of the treated unit.
#' @param outcomes Character vector or NULL. Names of the outcome variables to plot. If NULL, all outcomes in the data will be plotted.
#' @param normalize_outcomes Character. Method for normalizing treatment effects for comparability.
#'   Options: "none" (default), "percent_of_mean", "percent_of_preperiod", "standardized", "rmspe_ratio".
#'   - "percent_of_mean": (tau / mean(pre-period treated)) * 100
#'   - "percent_of_preperiod": (tau / final_preperiod_value) * 100
#'   - "standardized": tau / sd(control_effects)
#'   - "rmspe_ratio": post_RMSPE / pre_RMSPE (Abadie-style)
#' @param rmse_threshold Numeric. Threshold for root mean square error to filter results. Default is Inf.
#' @param shap_values Data.table or NULL. External SHAP values from run_xgboost_shap_analysis().
#'   Should contain columns: unit, full_spec_id, feature_group, feature, shapley_value.
#'   If NULL and show_shap=TRUE, SHAP values will be computed internally. Default is NULL.
#' @param file_path_save Character or NA. File path to save the plot. If NA, plot is not saved. Default is NA.
#' @param width Numeric. Width of the saved plot in inches. Default is 6.
#' @param height Numeric. Height of the saved plot in inches. Default is 10.
#' @param show_pvalues Logical. Whether to include significance coloring for treatment effects.
#'   Uses Abadie p-values by default, or bootstrap p-values if available. Default is TRUE.
#' @param p_threshold Numeric. P-value threshold for significance coloring. Default is 0.05.
#' @param prefer_bootstrap_pvalues Logical. Whether to prefer bootstrap p-values over Abadie p-values
#'   when both are available. Default is FALSE (prefer Abadie).
#' @param test_statistic Character. Which test statistic p-values to use for coloring.
#'   Options: "rmse_ratio" (default), "treatment_effect", "normalized_te".
#'   Only applies to Abadie placebo inference which provides multiple test statistics.
#' @param null_distribution Character. The distribution to plot as the null/comparison.
#'   Options: "placebo" (default) or "bootstrap".
#'   - "placebo": Use placebo effects from control units as the null distribution.
#'   - "bootstrap": Use the bootstrapped null distribution for the treated unit.
#' @param crop_outliers Character or Numeric. Method for cropping outliers in Panel A.
#'   Options: "none" (default), "percentile", "iqr", "mad", or numeric vector c(ymin, ymax).
#'   - "percentile": Crop to 1st-99th percentile of treated unit effects
#'   - "iqr": Crop to Q1 - 1.5*IQR to Q3 + 1.5*IQR
#'   - "mad": Crop to median +/- 3*MAD (robust outlier detection)
#'   - c(ymin, ymax): Manual y-axis limits
#' @param sort_by Character. Method for sorting specifications on x-axis.
#'   Options: "tau" (default), "pvalue", "rmspe_ratio".
#'   - "tau": Sort by treatment effect magnitude (original behavior)
#'   - "pvalue": Sort by statistical significance (most significant first)
#'   - "rmspe_ratio": Sort by post/pre RMSPE ratio
#' @param filter_specs Named list or NULL. Filters to apply to specifications BEFORE SHAP computation and p-value calculation.
#'   Each element should be named by the feature group column and contain allowed values.
#'   Available filters: "constant" (TRUE/FALSE), "const" (simplex/lasso/ridge/pensynth), 
#'   "outcome_model" (none/augsynth/lasso/ridge/ols), "fw" (uniform/optimize), 
#'   "feat" (covariate aggregation labels), "data_sample" (all/most_similar).
#'   Examples: list(constant = "TRUE"), list(const = c("simplex", "lasso")), 
#'   list(constant = "TRUE", outcome_model = c("none", "augsynth")). Default is NULL (no filtering).
#' @param show_shap Logical. Whether to compute and display SHAP values. When TRUE:
#'   - If shap_values is provided: uses external SHAP values
#'   - If shap_values is NULL: computes SHAP values internally using XGBoost
#'   - Requires at least 3 unique specifications for internal computation
#'   - Automatically detects available specification features (outcome_model, const, fw, feat, data_sample, constant)
#'   Default is TRUE.
#' @param shap_config List. Configuration for SHAP computation and significance testing with elements:
#'   - compute_pvalues: Logical. Whether to calculate SHAP significance via placebo inference.
#'     Requires multi-unit SHAP data (treated_unit_only = FALSE in XGBoost config). Default is TRUE.
#'   - pvalue_type: Character. Method for SHAP significance testing: "absolute" or "signed".
#'     "absolute" ranks by |SHAP| values, "signed" uses two-sided rank-based test. Default is "absolute".
#' @param shap_label_type Character. How to display SHAP values in y-axis labels: "absolute" or "signed".
#'   - "absolute": Shows mean |SHAP| values, e.g., "ridge (0.123)"
#'   - "signed": Shows mean SHAP values with sign, e.g., "ridge (+0.089)" or "lasso (-0.045)"
#'   Default is "absolute".
#' @param show_predictions Logical. Whether to display predicted treatment effects from XGBoost models
#'   alongside actual treatment effects in Panel A. Uses leave-one-out cross-validation predictions
#'   for robust evaluation. Default is FALSE.
#' @param predictions Data.table or NULL. External LOO predictions from run_xgboost_shap_analysis().
#'   Should contain columns: unit, full_spec_id, predicted_loo, actual.
#'   If NULL and show_predictions=TRUE, predictions will be computed internally. Default is NULL.
#' @param xgboost_params List. Custom parameters for XGBoost model training. If NULL, uses defaults:
#'   list(objective='reg:squarederror', max_depth=10, eta=0.05, nrounds=500, subsample=0.8,
#'   colsample_bytree=0.8, nthread=1, seed=42, verbose=0). Common parameters to modify:
#'   nrounds (more = better fit but slower), max_depth (complexity), eta (learning rate).
#'   Default is NULL. If tune_xgboost is TRUE, xgboost_params must be NULL.
#' @param tune_xgboost Logical or NULL. If TRUE, tune XGBoost hyperparameters via
#'   k-fold CV on the treated unit before SHAP/LOO. If NULL, defaults to TRUE when
#'   xgboost_params is NULL and FALSE otherwise. Default is NULL.
#' @param xgboost_grid Data.frame, data.table, or list. Optional grid for tuning.
#'   Must include max_depth, eta, nrounds, subsample, colsample_bytree. Default is NULL.
#' @param xgboost_cv_folds Integer. Number of CV folds for tuning. Default is 5.
#'
#' @return List containing:
#' \itemize{
#'   \item final_plot: ggplot object representing the complete specification curve with Panel A (treatment effects) and Panel B (feature groups/SHAP)
#'   \item panel_a: ggplot object for Panel A only (treatment effects)
#'   \item panel_b: ggplot object for Panel B only (feature groups/SHAP)
#'   \item plot_data_p1: data.table with Panel A plotting data including columns: Unit Name, Estimate, RMSE, Specification, unit_type, p_value (if available)
#'   \item plot_data_p2: data.table with Panel B plotting data including columns: Specification, feature_group, feature, shapley_value (if SHAP computed), shap_pvalue (if significance computed)
#'   \item computed_shap: Complete results from internal SHAP computation (NULL if external shap_values provided or show_shap=FALSE). 
#'     Contains: results (feature importance), shapley (SHAP values), predictions (model predictions), models (trained XGBoost models), config (SHAP configuration)
#'   \item spec_curve_pvals: List with specification curve-level p-values calculated on filtered data:
#'     median_tau_pvalues (median treatment effect ranks), stouffer_pvalues (Stouffer's Z-method ranks). NULL if no inference data available
#'   \item filtered_specs: Integer count of specifications remaining after filtering (before filtering if filter_specs=NULL)
#'   \item feature_groups_displayed: Character vector of feature groups shown in Panel B (only groups with variation in filtered data)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # ===== RECOMMENDED WORKFLOW: Internal SHAP Computation =====
#' 
#' # 1. Generate specification curve data
#' spec_results <- run_spec_curve_analysis(
#'   dataset = your_data,
#'   params = your_params,
#'   inference_type = "placebo"
#' )
#'
#' # 2. Basic plot with internal SHAP computation (most common use case)
#' plot_basic <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   show_shap = TRUE,  # Automatically computes SHAP internally
#'   shap_config = list(
#'     compute_pvalues = TRUE,
#'     pvalue_type = "absolute"
#'   )
#' )
#'
#' # ===== SPECIFICATION FILTERING EXAMPLES =====
#'
#' # Filter by constant terms only
#' plot_constant_only <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   filter_specs = list(constant = "TRUE"),
#'   show_shap = TRUE
#' )
#'
#' # Filter by multiple criteria (AND logic)
#' plot_filtered <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   filter_specs = list(
#'     constant = "TRUE",
#'     const = c("simplex", "lasso"),
#'     outcome_model = c("none", "augsynth")
#'   ),
#'   show_shap = TRUE
#' )
#'
#' # Compare specific constraint methods
#' plot_constraints <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   filter_specs = list(const = c("simplex", "lasso", "ridge")),
#'   show_shap = TRUE,
#'   shap_config = list(compute_pvalues = FALSE)  # Disable SHAP significance
#' )
#'
#' # ===== DIFFERENT TEST STATISTICS =====
#'
#' # Using treatment effect p-values
#' plot_treatment_effect <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   test_statistic = "treatment_effect",
#'   show_shap = TRUE
#' )
#'
#' # Using RMSE ratio p-values (default)
#' plot_rmse_ratio <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   test_statistic = "rmse_ratio",
#'   show_shap = TRUE
#' )
#'
#' # ===== SHAP SIGNIFICANCE TESTING =====
#'
#' # SHAP with signed significance testing
#' plot_shap_signed <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   show_shap = TRUE,
#'   shap_config = list(
#'     compute_pvalues = TRUE,
#'     pvalue_type = "signed"  # Two-sided SHAP significance
#'   )
#' )
#'
#' # ===== EXTERNAL SHAP WORKFLOW (Advanced) =====
#'
#' # For advanced users who want external control over SHAP computation
#' shap_config_external <- create_xgboost_config(
#'   dataset_name = "example",
#'   treated_unit_name = "TREATED_ID",
#'   treated_unit_only = FALSE  # Multi-unit for significance testing
#' )
#' external_shap <- run_xgboost_shap_analysis(spec_results, shap_config_external)
#'
#' plot_external_shap <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = external_shap$shapley,  # Provide external SHAP
#'   show_shap = TRUE
#' )
#'
#' # ===== OTHER OPTIONS =====
#'
#' # Disable SHAP entirely
#' plot_no_shap <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   show_shap = FALSE
#' )
#'
#' # Crop outliers and sort by p-values
#' plot_customized <- plot_spec_curve(
#'   long_data = spec_results,
#'   name_treated_unit = "TREATED_ID",
#'   show_shap = TRUE,
#'   crop_outliers = "percentile",
#'   sort_by = "pvalue",
#'   normalize_outcomes = "standardized"
#' )
#'
#' # ===== EXTRACTING INDIVIDUAL PANELS =====
#'
#' # Extract and display only Panel A (treatment effects)
#' 
#' # Extract and save Panel B (specification features/SHAP)
#'
#' }
plot_spec_curve <- function(
    long_data,
    name_treated_unit,
    outcomes = NULL,
    normalize_outcomes = "none",
    rmse_threshold = Inf,
    shap_values = NULL,
    file_path_save = NA,
    width = 6,
    height = 10,
    show_pvalues = FALSE,
    p_threshold = 0.05,
    prefer_bootstrap_pvalues = FALSE,
    test_statistic = "rmse_ratio",
    null_distribution = "placebo",
    crop_outliers = "none",
    sort_by = "tau",
    filter_specs = NULL,
    show_shap = TRUE,
    shap_config = list(
        compute_pvalues = TRUE,
        pvalue_type = "absolute"
    ),
    shap_label_type = "absolute",
    show_predictions = FALSE,
    predictions = NULL,
    xgboost_params = NULL,
    tune_xgboost = NULL,
    xgboost_grid = NULL,
    xgboost_cv_folds = 5
) {

    # Libraries imported via NAMESPACE
    
    # Define color constants to avoid redundancy
    SHAP_COLORS <- list(
        red = "#CA0020",      # Red for negative SHAP
        gray = "#969696",     # Gray for near-zero SHAP
        blue = "#0571B0",     # Blue for positive SHAP
        light_gray = "#F0F0F0" # Light gray for low absolute SHAP
    )
    
    # Helper function to map SHAP values to colors using same scale as legend
    map_shap_to_color <- function(shap_values, shap_label_type, all_shap_values = NULL) {
        if (is.null(all_shap_values)) {
            all_shap_values <- shap_values
        }
        
        if (shap_label_type == "absolute") {
            # Map absolute values using same scale as gradient: light_gray to blue
            abs_values <- abs(shap_values)
            abs_range <- range(abs(all_shap_values), na.rm = TRUE)
            if (abs_range[2] == abs_range[1]) {
                return(rep(SHAP_COLORS$gray, length(shap_values)))
            }
            # Normalize to 0-1 range
            normalized <- (abs_values - abs_range[1]) / (abs_range[2] - abs_range[1])
            # Interpolate between light_gray and blue
            colors <- grDevices::colorRamp(c(SHAP_COLORS$light_gray, SHAP_COLORS$blue))(normalized)
            return(grDevices::rgb(colors[,1], colors[,2], colors[,3], maxColorValue = 255))
        } else {
            # Map signed values using same scale as gradient2: red to gray to blue (midpoint = 0)
            shap_range <- range(all_shap_values, na.rm = TRUE)
            if (shap_range[2] == shap_range[1]) {
                return(rep(SHAP_COLORS$gray, length(shap_values)))
            }
            
            # Create symmetric range around 0 for proper gradient2 behavior
            max_abs <- max(abs(shap_range), na.rm = TRUE)
            
            # Vectorized color mapping for signed values
            result_colors <- character(length(shap_values))
            
            for (i in seq_along(shap_values)) {
                val <- shap_values[i]
                if (is.na(val)) {
                    result_colors[i] <- SHAP_COLORS$gray
                } else if (val == 0) {
                    result_colors[i] <- SHAP_COLORS$gray
                } else if (val < 0) {
                    # Map negative values from red to gray
                    proportion <- abs(val) / max_abs
                    color_matrix <- grDevices::colorRamp(c(SHAP_COLORS$gray, SHAP_COLORS$red))(proportion)
                    result_colors[i] <- grDevices::rgb(color_matrix[1], color_matrix[2], color_matrix[3], maxColorValue = 255)
                } else {
                    # Map positive values from gray to blue
                    proportion <- val / max_abs
                    color_matrix <- grDevices::colorRamp(c(SHAP_COLORS$gray, SHAP_COLORS$blue))(proportion)
                    result_colors[i] <- grDevices::rgb(color_matrix[1], color_matrix[2], color_matrix[3], maxColorValue = 255)
                }
            }
            
            return(result_colors)
        }
    }

    # Helper function for filtering specifications
    apply_spec_filters <- function(data, filter_specs) {
        if (!data.table::is.data.table(data)) {
            data <- data.table::as.data.table(data)
        }
        
        filtered_data <- data.table::copy(data)
        
        # Apply each filter sequentially
        for (filter_col in names(filter_specs)) {
            filter_values <- filter_specs[[filter_col]]
            
            # Check if the filter column exists in the data
            if (!filter_col %in% names(filtered_data)) {
                warning(paste("Filter column", filter_col, "not found in data. Available columns:",
                             paste(names(filtered_data), collapse = ", ")))
                next
            }
            
            # Apply the filter
            initial_rows <- nrow(filtered_data)
            filtered_data <- filtered_data[get(filter_col) %in% filter_values]
            final_rows <- nrow(filtered_data)
            
            # Provide feedback about filtering effect
            message(paste("Filter", filter_col, ":", initial_rows, "->", final_rows, "specifications"))
            
            # Check if filtering eliminated all data
            if (nrow(filtered_data) == 0) {
                stop(paste("Filter", filter_col, "eliminated all specifications. No data remaining."))
            }
        }
        
        message(paste("Total specifications after filtering:", nrow(filtered_data)))
        return(filtered_data)
    }

    # Input validation
    valid_test_statistics <- c("rmse_ratio", "treatment_effect", "normalized_te")
    if (!test_statistic %in% valid_test_statistics) {
        stop("test_statistic must be one of: ", paste(valid_test_statistics, collapse = ", "))
    }
    # Validate filter_specs structure
    if (!is.null(filter_specs) && !is.list(filter_specs)) {
        stop("filter_specs must be a named list or NULL")
    }
    
    # Validate shap_config structure
    if (!is.list(shap_config)) {
        stop("shap_config must be a list")
    }
    
    # Set default shap_config values
    if (is.null(shap_config$compute_pvalues)) {
        shap_config$compute_pvalues <- TRUE
    }
    if (is.null(shap_config$pvalue_type)) {
        shap_config$pvalue_type <- "absolute"
    }
    
    # Validate shap_config$pvalue_type
    valid_shap_pvalue_types <- c("absolute", "signed")
    if (!shap_config$pvalue_type %in% valid_shap_pvalue_types) {
        stop("shap_config$pvalue_type must be one of: ", paste(valid_shap_pvalue_types, collapse = ", "))
    }
    
    # Validate shap_label_type
    valid_shap_label_types <- c("absolute", "signed")
    if (!shap_label_type %in% valid_shap_label_types) {
        stop("shap_label_type must be one of: ", paste(valid_shap_label_types, collapse = ", "))
    }

    # Input validation and data extraction
    if (is.list(long_data) && "results" %in% names(long_data)) {
        # Handle structured results from spec_curve
        main_data <- data.table::as.data.table(long_data$results)

        # Merge inference results for the specified test statistic
        p_values_key <- paste0("p_values_", test_statistic)
        test_stats_key <- paste0("test_statistics_", test_statistic)


        if (!is.null(long_data$abadie_inference[[p_values_key]])) {
            abadie_pvals <- long_data$abadie_inference[[p_values_key]]
            # limit to treated unit
            abadie_pvals = abadie_pvals[unit_type=='treated']
            main_data = merge(main_data, abadie_pvals[, .(full_spec_id, p_value)],
                by = "full_spec_id", all.x = TRUE)
        }

        # Merge test statistics values using unit-level keys to avoid cartesian joins
        if (!is.null(long_data$abadie_inference[[test_stats_key]])) {
            test_stats_data <- long_data$abadie_inference[[test_stats_key]]
            # Include unit_type to ensure unique merging (avoids cartesian joins)
            merge_cols <- intersect(c("full_spec_id", "unit_name", "unit_type"), names(test_stats_data))
            if (length(merge_cols) >= 2) {  # Need at least full_spec_id + unit identifier
                # Get the appropriate column name for the test statistic value
                if ("test_statistic_value" %in% names(test_stats_data)) {
                    value_col <- "test_statistic_value"
                } else if ("post_pre_ratio" %in% names(test_stats_data)) {
                    value_col <- "post_pre_ratio"
                } else {
                    value_col <- NULL
                }

                if (!is.null(value_col)) {
                    main_data = merge(main_data, test_stats_data[, c(merge_cols, value_col), with = FALSE],
                        by = merge_cols, all.x = TRUE)
                    # Rename to standardized column name for backward compatibility
                    if (value_col == "test_statistic_value" && test_statistic == "rmse_ratio") {
                        setnames(main_data, "test_statistic_value", "post_pre_ratio")
                    }
                }
            } else {
                warning("Insufficient merge columns for test statistics - skipping merge")
            }
        }

        if (!is.null(long_data$bootstrap_inference$p_values)) {
            bootstrap_pvals <- long_data$bootstrap_inference$p_values
            main_data = merge(main_data, bootstrap_pvals[, .(full_spec_id, p_value_two_tailed )],
                by = "full_spec_id", all.x = TRUE)
        }


        # Add bootstrap iteration data for null distribution plotting
        if (!is.null(long_data$bootstrap_inference$iteration_data)) {
            bootstrap_iteration_data <- long_data$bootstrap_inference$iteration_data

            # Check if bootstrap iteration data exists and has data
            if (length(bootstrap_iteration_data) > 0) {
                # Combine all outcome models' bootstrap data
                bootstrap_combined <- rbindlist(bootstrap_iteration_data, fill = TRUE)

                if (nrow(bootstrap_combined) > 0) {
                    # Add necessary columns to match main_data structure
                    if (!"rmse" %in% names(bootstrap_combined)) {
                        bootstrap_combined[, rmse := NA_real_]
                    }

                    # Bootstrap data needs specification metadata - get unique spec metadata from main_data
                    spec_metadata <- unique(main_data[, .(full_spec_id, outcome, outcome_model, const, fw, feat, data_sample)])

                    # For each bootstrap row, we need to replicate it across all specifications for that outcome_model
                    bootstrap_expanded <- list()
                    for (i in 1:nrow(bootstrap_combined)) {
                        boot_row <- bootstrap_combined[i]
                        matching_specs <- spec_metadata[outcome_model == boot_row$outcome_model]

                        if (nrow(matching_specs) > 0) {
                            # Replicate this bootstrap observation for each matching specification
                            expanded_rows <- boot_row[rep(1, nrow(matching_specs))]
                            expanded_rows <- cbind(expanded_rows, matching_specs)
                            bootstrap_expanded[[i]] <- expanded_rows
                        }
                    }

                    if (length(bootstrap_expanded) > 0) {
                        bootstrap_final <- rbindlist(bootstrap_expanded, fill = TRUE)
                        main_data <- rbindlist(list(main_data, bootstrap_final), fill = TRUE)
                    }
                }
            }
        }

        # Keep main_data as the processed results data
    } else {
        # Handle direct data.table input (not structured from spec_curve)
        if (!data.table::is.data.table(long_data)) {
            long_data <- data.table::as.data.table(long_data)
        }
        main_data <- long_data
    }

    # Check required columns
    required_cols <- c("unit_name", "tau", "post_period", "full_spec_id")
    missing_cols <- setdiff(required_cols, names(main_data))
    if (length(missing_cols) > 0) {
        stop(paste("Missing required columns in main_data:", paste(missing_cols, collapse = ", ")))
    }

    # Check if p-value data is available
    has_abadie_pvalues <- "p_value" %in% names(main_data)
    has_bootstrap_pvalues <- "p_value_two_tailed" %in% names(main_data)
    has_pvalues <- has_abadie_pvalues || has_bootstrap_pvalues


    # Determine which p-values to use
    # When using bootstrap null distribution, prefer bootstrap p-values
    if (null_distribution == "bootstrap" && has_bootstrap_pvalues) {
        use_bootstrap_pvalues <- TRUE
    } else if (prefer_bootstrap_pvalues && has_bootstrap_pvalues) {
        use_bootstrap_pvalues <- TRUE
    } else if (!has_abadie_pvalues && has_bootstrap_pvalues) {
        use_bootstrap_pvalues <- TRUE
    } else {
        use_bootstrap_pvalues <- FALSE
    }

    p_value_column <- if (use_bootstrap_pvalues) "p_value_two_tailed" else "p_value"

    # Apply filtering pipeline in proper order
    sc_results_df <- main_data
    
    # 1. Filter by outcomes (if specified)
    if (!is.null(outcomes) && "outcome" %in% names(sc_results_df)) {
        sc_results_df <- sc_results_df[outcome %in% outcomes]
    }

    # 2. Filter by RMSE if threshold provided
    if ("rmse" %in% names(sc_results_df) && rmse_threshold < Inf) {
        sc_results_df <- sc_results_df[rmse < rmse_threshold]
    }
    
    # 3. Apply specification filters BEFORE p-value calculations
    if (!is.null(filter_specs)) {
        sc_results_df <- apply_spec_filters(sc_results_df, filter_specs)
    }

    if (nrow(sc_results_df) == 0) {
        stop("No data remaining after filtering")
    }

    # Calculate average effects for plotting (post-treatment period only)
    # Define common grouping variables for aggregation (only include columns that exist)
    agg_grouping_vars_candidates <- c('unit_name', 'unit_type', 'full_spec_id', 'outcome',
                          'outcome_model', 'const', 'fw', 'feat', 'data_sample', 'constant')
    agg_grouping_vars <- intersect(agg_grouping_vars_candidates, names(sc_results_df))
    
    # Perform aggregation with conditional columns
    if ("post_pre_ratio" %in% names(sc_results_df)) {
        average_effect_df <- sc_results_df[post_period == TRUE, 
                                         list(tau = mean(tau, na.rm = TRUE),
                                              rmse = mean(rmse, na.rm = TRUE),
                                              post_pre_ratio = mean(post_pre_ratio, na.rm = TRUE)), 
                                         by = agg_grouping_vars]
    } else {
        average_effect_df <- sc_results_df[post_period == TRUE, 
                                         list(tau = mean(tau, na.rm = TRUE),
                                              rmse = mean(rmse, na.rm = TRUE)), 
                                         by = agg_grouping_vars]
    }


    # Normalize outcomes if requested
    if (normalize_outcomes != "none") {

        if (normalize_outcomes == "standardized") {
            # Standardize by comparison unit standard deviation (control or bootstrap)
            comparison_units <- average_effect_df[unit_type %in% c("control", "bootstrap")]
            sd_outcome <- sd(comparison_units$tau, na.rm = TRUE)
            if (sd_outcome > 0) {
                average_effect_df[, tau := tau / sd_outcome]
            }
            y_label <- "Standardized Treatment Effect"

        } else if (normalize_outcomes == "percent_of_mean") {
            # For percentage normalization, we need the outcome scale
            # Since tau = actual - synthetic, use RMSE as a proxy for outcome scale
            treated_rmse <- average_effect_df[unit_name == name_treated_unit, mean(rmse, na.rm = TRUE)]
            if (!is.na(treated_rmse) && treated_rmse > 0) {
                # Normalize by typical outcome scale (approximated by RMSE)
                average_effect_df[unit_name == name_treated_unit, tau := (tau / treated_rmse) * 100]
            }
            y_label <- "Treatment Effect (% of Outcome Scale)"

        } else if (normalize_outcomes == "percent_of_preperiod") {
            # Similar approach using RMSE as scale
            treated_rmse <- average_effect_df[unit_name == name_treated_unit, mean(rmse, na.rm = TRUE)]
            if (!is.na(treated_rmse) && treated_rmse > 0) {
                average_effect_df[unit_name == name_treated_unit, tau := (tau / treated_rmse) * 100]
            }
            y_label <- "Treatment Effect (% of Outcome Scale)"
        } else if (normalize_outcomes == "rmspe_ratio") {
            # Use post/pre RMSPE ratio if available
            if ("post_pre_ratio" %in% names(average_effect_df)) {
                average_effect_df[, tau := post_pre_ratio]
            }
            y_label <- "Post/Pre RMSPE Ratio"
        } else {
            y_label <- "Average Treatment Effect"
        }
    } else {
        y_label <- "Average Treatment Effect"
    }

    # Calculate y-axis limits for outlier cropping
    y_limits <- NULL
    if (crop_outliers != "none") {
        treated_effects <- average_effect_df[unit_name == name_treated_unit, tau]

        if (is.numeric(crop_outliers) && length(crop_outliers) == 2) {
            # Manual limits
            y_limits <- crop_outliers

        } else if (crop_outliers == "percentile") {
            # 5th to 95th percentile
            y_limits <- quantile(treated_effects, c(0.01, 0.99), na.rm = TRUE)

        } else if (crop_outliers == "iqr") {
            # Interquartile range with 1.5*IQR extension
            q1 <- quantile(treated_effects, 0.25, na.rm = TRUE)
            q3 <- quantile(treated_effects, 0.75, na.rm = TRUE)
            iqr <- q3 - q1
            y_limits <- c(q1 - 1.5 * iqr, q3 + 1.5 * iqr)

        } else if (crop_outliers == "mad") {
            # Median +/- 3 * Median Absolute Deviation
            med <- median(treated_effects, na.rm = TRUE)
            mad_val <- mad(treated_effects, na.rm = TRUE)
            y_limits <- c(med - 3 * mad_val, med + 3 * mad_val)
        }
    }

    # --- 1. Data Preparation ---

    # A. Prepare data for Panel A and create numbered specifications
    panel_a_data <- copy(average_effect_df)

    # Create specification numbering based on treated unit effects (sorted by normalized tau)
    treated_specs <- panel_a_data[unit_name == name_treated_unit][order(tau)]
    treated_specs[, Specification := 1:.N]
    spec_mapping <- treated_specs[, .(full_spec_id, Specification)]

    # Apply specification numbering to all data
    panel_a_data <- merge(panel_a_data, spec_mapping, by = "full_spec_id", all.x = TRUE)
    setnames(panel_a_data,
             old = c("unit_name", "tau", "rmse"),
             new = c("Unit Name", "Estimate", "RMSE"))

    # Add p-values if requested
    if (has_pvalues && show_pvalues) {
        p_value_data <- unique(sc_results_df[unit_name == name_treated_unit & !is.na(get(p_value_column)),
                                             .(full_spec_id, p_value = get(p_value_column))])
        p_value_data <- merge(p_value_data, spec_mapping, by = "full_spec_id", all.x = TRUE)

        panel_a_data <- merge(panel_a_data, p_value_data[, .(Specification, p_value)], by = "Specification", all.x = TRUE)
    }
    


    # B. Prepare data for Panel B (SHAP & Specification Details)
    if ("post_pre_ratio" %in% names(average_effect_df)) {
        average_effect_df[, post_pre_ratio := NULL]
    }
    
    # Convert constant column to character to avoid melt type warning
    if ("constant" %in% names(average_effect_df)) {
        average_effect_df[, constant := as.character(constant)]
    }
    
    panel_b_data <- melt(average_effect_df[unit_name == name_treated_unit],
                         id.vars = c('unit_name', 'full_spec_id', 'tau', 'rmse', 'unit_type'),
                         variable.name = 'feature_group', value.name = 'feature')

    # Apply the same specification numbering
    panel_b_data <- merge(panel_b_data, spec_mapping, by = "full_spec_id", all.x = TRUE)

    # Handle SHAP computation and merging
    computed_shap <- NULL
    
    # Check if we need to compute SHAP internally
    if (show_shap && is.null(shap_values)) {
        if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("Package 'xgboost' is required for internal SHAP computation. Install with: install.packages('xgboost') or provide shap_values.", call. = FALSE)
        }
        message("Internal SHAP computation requested. Running XGBoost SHAP analysis...")
        
        # Create default configuration for internal SHAP computation
        # Use the most common specifications for internal analysis
        default_spec_features <- c("outcome_model", "const", "fw", "feat", "data_sample")
        
        # Add 'constant' to spec features if it exists in the data and has variation
        if ("constant" %in% names(sc_results_df)) {
            constant_variation <- data.table::uniqueN(sc_results_df$constant) > 1
            if (constant_variation) {
                default_spec_features <- c(default_spec_features, "constant")
                message("Adding 'constant' to SHAP features - detected variation in constant terms")
            }
        }
        
        # Filter to available spec features only
        available_spec_features <- intersect(default_spec_features, names(sc_results_df))
        
        if (length(available_spec_features) == 0) {
            stop("Internal SHAP computation failed: No valid specification features found in data. ",
                 "Available columns: ", paste(names(sc_results_df), collapse = ", "), ". ",
                 "Either provide shap_values parameter or ensure data contains specification features.")
        }
        
        # Detect outcome filter from filtered data  
        outcome_filter <- NULL
        if ("outcome" %in% names(sc_results_df)) {
            unique_outcomes <- unique(sc_results_df$outcome)
            if (length(unique_outcomes) == 1) {
                outcome_filter <- unique_outcomes[1]
            }
        }
        
        # Create configuration for internal SHAP computation
        shap_config_internal <- create_xgboost_config(
            dataset_name = paste0("internal_shap_", Sys.time()),
            treated_unit_name = name_treated_unit,
            outcome_filter = outcome_filter,
            spec_features = available_spec_features,
            treated_unit_only = TRUE,  # Default to treated unit only for efficiency
            xgboost_params = xgboost_params,
            tune_xgboost = tune_xgboost,
            xgboost_grid = xgboost_grid,
            xgboost_cv_folds = xgboost_cv_folds
        )
        
        # Create long format data structure expected by run_xgboost_shap_analysis
        # Add spec_number for proper alignment 
        long_format_data <- copy(sc_results_df)
        
        # Check minimum specifications required for XGBoost
        treated_unit_specs <- unique(sc_results_df[unit_name == name_treated_unit, full_spec_id])
        n_unique_specs <- length(treated_unit_specs)
        
        if (n_unique_specs < 3) {
            stop("Internal SHAP computation failed: Insufficient specifications for XGBoost analysis. ",
                 "Found ", n_unique_specs, " unique specifications, but at least 3 are required. ",
                 "Either provide external shap_values, set show_shap=FALSE, or run with more specification variations.")
        }
        
        # Add spec_number if not present (based on treated unit ordering like in plot creation)
        if (!"spec_number" %in% names(long_format_data)) {
            treated_specs <- long_format_data[unit_name == name_treated_unit & post_period == TRUE]
            if (nrow(treated_specs) > 0) {
                treated_specs_ordered <- treated_specs[order(tau)]
                treated_specs_ordered[, spec_number := 1:.N]
                spec_number_mapping <- treated_specs_ordered[, .(full_spec_id, spec_number)]
                
                # Apply to all data
                long_format_data <- merge(long_format_data, spec_number_mapping, 
                                        by = "full_spec_id", all.x = TRUE)
            } else {
                stop("Internal SHAP computation failed: No treated unit data found for spec_number creation")
            }
        }
        
        # Run internal SHAP analysis
        shap_results_internal <- tryCatch({
            run_xgboost_shap_analysis(long_format_data, shap_config_internal, compute_loo = show_predictions)
        }, error = function(e) {
            stop("Internal SHAP computation failed: ", e$message, ". ",
                 "Either provide shap_values parameter or set show_shap=FALSE.")
        })
        
        # Extract SHAP values and store computed results
        if (!is.null(shap_results_internal) && !is.null(shap_results_internal$shapley)) {
            shap_values <- shap_results_internal$shapley
            computed_shap <- shap_results_internal
            message("Internal SHAP computation completed successfully. ", 
                    nrow(shap_values), " SHAP observations computed.")
        } else {
            stop("Internal SHAP computation failed: No SHAP values returned. ",
                 "Either provide shap_values parameter or set show_shap=FALSE.")
        }
    }
    
    # Merge SHAP values if available (either provided or computed)
    if (!is.null(shap_values)) {
        # Validate SHAP data structure - FAIL HARD if incorrect
        required_shap_cols <- c('unit', 'full_spec_id', 'feature_group', 'feature', 'shapley_value')
        missing_shap_cols <- setdiff(required_shap_cols, names(shap_values))
        if (length(missing_shap_cols) > 0) {
            stop("Invalid SHAP data structure. Missing required columns: ", 
                 paste(missing_shap_cols, collapse = ", "), ". ",
                 "SHAP data must contain columns: ", paste(required_shap_cols, collapse = ", "))
        }
        
        # Check for multi-unit SHAP data and calculate significance if requested
        unique_units <- unique(shap_values$unit)
        has_control_shaps <- length(unique_units) > 1 && any(unique_units != name_treated_unit)

        if (has_control_shaps && shap_config$compute_pvalues) {
            message("Found SHAP values for ", length(unique_units), " units (",
                    sum(unique_units != name_treated_unit), " controls). Calculating SHAP significance using '", shap_config$pvalue_type, "' method...")
            message("To disable SHAP significance testing, set shap_config$compute_pvalues = FALSE")

            # Calculate SHAP significance on filtered data
            shap_significance <- calculate_shap_significance(shap_values, name_treated_unit, pvalue_method = shap_config$pvalue_type)
            # Merge significance data with SHAP values
            shap_values <- merge(shap_values, shap_significance,
                                 by = c("full_spec_id", "feature_group"), all.x = TRUE)

            # Provide user feedback about significance results
            treated_significance <- shap_significance[!is.na(shap_pvalue)]
            if (nrow(treated_significance) > 0) {
                n_significant <- sum(treated_significance$shap_pvalue < 0.05, na.rm = TRUE)
                n_marginal <- sum(treated_significance$shap_pvalue >= 0.05 & treated_significance$shap_pvalue < 0.10, na.rm = TRUE)
                total_tests <- nrow(treated_significance)
                message("SHAP significance: ", n_significant, " significant (p<0.05), ",
                        n_marginal, " marginal (p<0.10) out of ", total_tests, " tests")
            }
        } else if (has_control_shaps && !shap_config$compute_pvalues) {
            message("Multi-unit SHAP data found but significance testing disabled by shap_config$compute_pvalues = FALSE")
        }

        # Debug: Show data structures before merging
        message("DEBUG: SHAP merging diagnostics...")
        message("Panel B data structure (first 5 rows for treated unit):")
        panel_b_treated <- panel_b_data[unit_name == name_treated_unit][1:min(5, .N)]
        if (nrow(panel_b_treated) > 0) {
            message("  Feature groups: ", paste(unique(panel_b_treated$feature_group), collapse = ", "))
            message("  Sample features: ", paste(head(panel_b_treated$feature, 5), collapse = ", "))
        }
        
        message("SHAP values structure (first 5 rows for treated unit):")
        shap_treated <- shap_values[unit == name_treated_unit][1:min(5, .N)]
        if (nrow(shap_treated) > 0) {
            message("  Feature groups: ", paste(unique(shap_treated$feature_group), collapse = ", "))
            message("  Sample features: ", paste(head(shap_treated$feature, 5), collapse = ", "))
        }

        # FIXED: Aggregate SHAP values to specification level before merging
        # SHAP values are at individual feature level, but panel_b needs them at specification level
        message("Aggregating SHAP values to specification level...")
        
        # Aggregate SHAP values by specification (full_spec_id) and feature_group 
        # This sums up all SHAP contributions within each specification dimension
        shap_aggregated <- shap_values[, .(
            shapley_value = sum(shapley_value, na.rm = TRUE),
            n_features = .N
        ), by = .(unit, full_spec_id, feature_group)]
        
        message("SHAP aggregation: ", nrow(shap_values), " individual feature SHAP values -> ", 
                nrow(shap_aggregated), " specification-level SHAP values")
        
        # For each aggregated SHAP value, we need to determine the corresponding 'feature' value
        # This is the actual categorical value for that specification dimension
        panel_b_treated <- panel_b_data[unit_name == name_treated_unit]
        
        if (nrow(panel_b_treated) > 0) {
            # Create lookup for feature values by full_spec_id and feature_group
            feature_lookup <- unique(panel_b_treated[, .(full_spec_id, feature_group, feature)])
            
            # Merge feature values into aggregated SHAP data
            shap_aggregated <- merge(shap_aggregated, feature_lookup,
                                   by = c("full_spec_id", "feature_group"),
                                   all.x = TRUE)
        }
        
        # Perform the merge with aggregated SHAP values
        panel_b_data_before_merge <- copy(panel_b_data)
        panel_b_data <- merge(panel_b_data, shap_aggregated,
                              by.x = c('unit_name', 'full_spec_id', 'feature_group', 'feature'),
                              by.y = c('unit', 'full_spec_id', 'feature_group', 'feature'),
                              all.x = TRUE)
        
        # Debug: Check merge success rate on mergeable rows only.
        # Panel B may contain feature groups that are intentionally absent from SHAP
        # (e.g., non-varying groups removed before model fitting), so the denominator
        # should be restricted to groups that exist in SHAP for the treated unit.
        treated_before <- panel_b_data_before_merge[unit_name == name_treated_unit]
        treated_after <- panel_b_data[unit_name == name_treated_unit]
        shap_groups_treated <- unique(shap_aggregated[unit == name_treated_unit]$feature_group)
        
        if (length(shap_groups_treated) == 0) {
            warning("No treated-unit SHAP feature groups available after aggregation. ",
                    "Cannot assess SHAP merge coverage.")
        } else {
            treated_mergeable_before <- treated_before[feature_group %in% shap_groups_treated]
            treated_mergeable_after <- treated_after[feature_group %in% shap_groups_treated]
            mergeable_rows_before <- nrow(treated_mergeable_before)
            mergeable_rows_with_shap <- nrow(treated_mergeable_after[!is.na(shapley_value)])
            
            if (mergeable_rows_before == 0) {
                warning("No mergeable Panel B rows found for treated unit across SHAP feature groups. ",
                        "Feature group naming mismatch may exist.")
            } else {
                merge_success_rate <- mergeable_rows_with_shap / mergeable_rows_before * 100
                message("SHAP merge success (mergeable rows): ",
                        mergeable_rows_with_shap, "/", mergeable_rows_before,
                        " (", round(merge_success_rate, 1), "%)")
                
                if (merge_success_rate < 50) {
                    warning("Low SHAP merge success rate on mergeable rows (",
                            round(merge_success_rate, 1), "%). ",
                            "Feature group/value mismatch may still exist between panel_b_data and aggregated SHAP values.")
                } else {
                    message("SHAP merging successful on mergeable rows.")
                }
            }
        }
    }

    setnames(panel_b_data, c('tau', 'unit_name', 'rmse'),
             c('Estimate', 'Unit Name', 'RMSE'))

    # Add predicted treatment effects if requested and available (after SHAP computation)
    # Use external predictions if provided, otherwise fall back to internally computed
    predictions_source <- if (!is.null(predictions)) predictions
                          else if (!is.null(computed_shap)) computed_shap$predictions
                          else NULL
    if (show_predictions && !is.null(predictions_source)) {
        # Extract predictions for treated unit
        predictions_data <- predictions_source[unit == name_treated_unit]
        
        if (nrow(predictions_data) > 0 && "predicted_loo" %in% names(predictions_data)) {
            # Merge predictions with specification mapping
            if ("full_spec_id" %in% names(predictions_data)) {
                pred_with_spec <- merge(predictions_data[, list(full_spec_id, predicted_loo)], 
                                      spec_mapping, by = "full_spec_id", all.x = TRUE)
                
                # Add to panel data
                panel_a_data <- merge(panel_a_data, pred_with_spec[, list(Specification, Predicted = predicted_loo)], 
                                    by = "Specification", all.x = TRUE)
            }
        }
    }

    # C. Re-sort specifications if requested (default is by tau, which is already done)
    if (sort_by != "tau") {
        # Create new ordering based on sort_by parameter
        treated_data <- panel_a_data[unit_type == "treated"]

        if (sort_by == "pvalue" && "p_value" %in% names(treated_data)) {
            treated_data <- treated_data[order(p_value)]
        } else if (sort_by == "rmspe_ratio" && "post_pre_ratio" %in% names(treated_data)) {
            treated_data <- treated_data[order(post_pre_ratio)]
        }

        # Create new specification mapping
        treated_data[, new_specification := 1:.N]
        spec_reorder <- treated_data[, .(old_spec = Specification, new_specification)]

        # Apply new ordering to both panel datasets
        panel_a_data <- merge(panel_a_data, spec_reorder, by.x = "Specification", by.y = "old_spec", all.x = TRUE)
        panel_a_data[, Specification := new_specification][, new_specification := NULL]

        panel_b_data <- merge(panel_b_data, spec_reorder, by.x = "Specification", by.y = "old_spec", all.x = TRUE)
        panel_b_data[, Specification := new_specification][, new_specification := NULL]
    }

    # D. Finalize Panel A data
    # Select the null distribution
    null_unit_type <- match.arg(null_distribution, choices = c("placebo", "bootstrap"))
    null_effects <- panel_a_data[unit_type == ifelse(null_unit_type == "placebo", "control", "bootstrap")]

    if (nrow(null_effects) == 0 && null_unit_type == "bootstrap") {
        warning("No bootstrap distribution found. Falling back to placebo units.")
        null_effects <- panel_a_data[unit_type == "control"]
    }

    treated_effects <- panel_a_data[unit_type == "treated"]
    plot_data_p1 <- rbindlist(list(treated_effects, null_effects), use.names = TRUE, fill = TRUE)

    # Set factor levels for plotting order
    plot_data_p1[, unit_type := factor(unit_type, levels = c("treated", "control", "bootstrap"))]

    # E. Finalize Panel B data
    # Clean up feature names for display
    panel_b_data[feature_group=='const', feature_group:= 'Weight\nMethod']
    panel_b_data[feature_group=='outcome', feature_group:=  'Outcome']
    panel_b_data[feature_group=='outcome_model', feature_group:= 'Outcome\nModel']
    panel_b_data[feature_group== 'fw', feature_group:= 'V Weights']
    panel_b_data[feature_group== 'feat', feature_group:= 'Features']
    panel_b_data[feature_group== 'data_sample', feature_group:='Donor Pool']
    panel_b_data[feature_group== 'const', feature_group:= 'Weight\nMethod']
    panel_b_data[feature_group== 'constant', feature_group:= 'Constant\nTerm']
    
    # Transform constraint names for display (preserve descriptive names in Features group)
    panel_b_data[feature == 'simplex' & feature_group=='Weight\nMethod', feature := "Original"]
    panel_b_data[feature == 'lasso' & feature_group=='Weight\nMethod', feature := "Penalty Lasso"]
    panel_b_data[feature == 'ridge' & feature_group=='Weight\nMethod', feature := "Penalty Ridge"]
    panel_b_data[feature == 'ols' & feature_group=='Weight\nMethod', feature := "OLS Weights"]
    
    # Transform constant term display values
    panel_b_data[feature == 'FALSE' & feature_group=='Constant\nTerm', feature := "No Constant"]
    panel_b_data[feature == 'TRUE' & feature_group=='Constant\nTerm', feature := "With Constant"]
    
    # DO NOT transform feature names in the 'Features' group - preserve descriptive names as-is

    # Dynamic feature group detection - only show groups with variation in filtered data
    feature_counts <- panel_b_data[, .(n_unique = data.table::uniqueN(feature)), by = feature_group]
    groups_to_keep <- feature_counts[n_unique > 1, feature_group]
    
    # Validate that feature groups have variation - FAIL HARD if not
    if (length(groups_to_keep) == 0) {
        stop("No feature groups have variation in filtered data. This indicates insufficient specification diversity. ",
             "Available feature groups: ", paste(unique(panel_b_data$feature_group), collapse = ", "), ". ",
             "Either expand your specification parameters or remove restrictive filters.")
    }
    
    message("Displaying feature groups with variation: ", paste(groups_to_keep, collapse = ", "))
    plot_data_p2 <- panel_b_data[feature_group %in% groups_to_keep]


    # --- 2. Create the Final Plots ---
    data.table::setorder(plot_data_p1, -unit_type)

    # Create Panel A with conditional coloring for p-values
    if (has_pvalues && show_pvalues && "p_value" %in% names(plot_data_p1)) {
        treated_data <- plot_data_p1[unit_type == "treated"]
        control_data <- plot_data_p1[unit_type %in% c("control", "bootstrap")]

        pvalue_style='adf'
        if (pvalue_style == "continuous") {
            # Continuous p-value coloring with significance boundary
            treated_data[, is_significant := p_value < p_threshold]

            p1 <- ggplot() +
                geom_point(data = control_data, aes(x = Specification, y = Estimate),
                           color = "gray60", alpha = 0.3, shape = 19, size = 2) +
                # Add significance boundary line for treated effects
                geom_point(data = treated_data[is_significant == TRUE],
                           aes(x = Specification, y = Estimate),
                           color = "black", alpha = 0.9, shape = 21, stroke = 1.2, size = 2) +
                geom_point(data = treated_data, aes(x = Specification, y = Estimate,
                                                    color = -log10(pmax(p_value, 1e-10))),
                           alpha = 0.8, shape = 19, size = 2) +
                geom_hline(yintercept = 0, alpha = 0.5, linetype = 'dashed') +
                scale_color_gradient2(
                    name = "-log10(p-value)",
                    low = "#d73027",      # Red for high p-values (non-significant)
                    mid = "#fee08b",      # Yellow for moderate p-values
                    high = "#1a9850",     # Green for low p-values (significant)
                    midpoint = -log10(p_threshold),   # Threshold at p=0.05
                    breaks = c(0, -log10(c(0.5, 0.1, 0.05, 0.01, 0.001))),
                    labels = c("1.0", "0.5", "0.1", "0.05", "0.01", "0.001"),
                    guide = guide_colorbar(
                        title.position = "top",
                        barwidth = grid::unit(4, "cm"),
                        barheight = grid::unit(0.5, "cm")
                    )
                )
        } else {
            # Categorical p-value coloring (original behavior)
            treated_data[, significance_category := fcase(
                p_value < 0.01, "p < 0.01",
                p_value < 0.05, "p < 0.05",
                p_value < 0.10, "p < 0.10",
                p_value >= 0.10, "p >= 0.10",
                default = "Unknown"
            )]

            p1 <- ggplot() +
                geom_point(data = control_data, aes(x = Specification, y = Estimate),
                           color = "gray60", alpha = 0.3, shape = 19, size = 2) +
                geom_point(data = treated_data, aes(x = Specification, y = Estimate, color = significance_category),
                           alpha = 0.8, shape = 19, size = 2) +
                geom_hline(yintercept = 0, alpha = 0.5, linetype = 'dashed') +
                scale_color_manual(
                    name = "p-value",
                    values = c("p < 0.01" = "#08519c", "p < 0.05" = "#3182bd",
                               "p < 0.10" = "#fd8d3c", "p >= 0.10" = "#d94701",
                               "Unknown" = "#999999"),
                    breaks = c("p < 0.01", "p < 0.05",
                               "p < 0.10", "p >= 0.10")
                )
        }
    } else {
        p1 <- ggplot(plot_data_p1, aes(x = Specification, y = Estimate, fill = unit_type,
                                       color = unit_type, alpha = unit_type)) +
            geom_point(shape = 21, size = 2) +
            geom_hline(yintercept = 0, alpha = 0.5, linetype = 'dashed') +
            scale_fill_manual(name = "Unit Type", values = c(control = "gray60", treated = "#1f78b4", bootstrap = "gray60")) +
            scale_color_manual(name = "Unit Type", values = c(control = "gray60", treated = "#1f78b4", bootstrap = "gray60")) +
            scale_alpha_manual(name = "Unit Type", values = c(treated = 0.8, control = 0.3, bootstrap = 0.3))
    }

    # Add predicted treatment effects if available (BEFORE actual points so they're underneath)
    if (show_predictions && "Predicted" %in% names(plot_data_p1)) {
        # Add predicted values for treated unit only
        treated_predictions <- plot_data_p1[unit_type == "treated" & !is.na(Predicted)]
        
        if (nrow(treated_predictions) > 0) {
            p1 <- p1 +
                # Add subtle connecting lines first (so they're in background)
                geom_segment(data = treated_predictions, 
                           aes(x = Specification, xend = Specification, 
                               y = Estimate, yend = Predicted),
                           color = "red", alpha = 0.2, linetype = "dotted") +
                # Add predicted points (more transparent, underneath actual points)
                geom_point(data = treated_predictions, aes(x = Specification, y = Predicted),
                          color = "red", shape = 4, size = 2.5, stroke = 1.2, alpha = 0.4)
        }
    }
    
    # Re-plot actual treated unit points on top to ensure they're visible
    if (show_predictions && "Predicted" %in% names(plot_data_p1)) {
        treated_actual <- plot_data_p1[unit_type == "treated"]
        if (nrow(treated_actual) > 0) {
            p1 <- p1 +
                geom_point(data = treated_actual, aes(x = Specification, y = Estimate),
                          color = "#1f78b4", shape = 21, size = 2.2, fill = "#1f78b4", alpha = 0.9, stroke = 0.5)
        }
    }

    # Add treatment/control/prediction legend via dummy data
    legend_data <- data.table::data.table(
        x = NA_real_, y = NA_real_,
        Legend = factor(c("Treated", "Control/Placebo"),
                        levels = c("Treated", "Control/Placebo"))
    )
    if (show_predictions && "Predicted" %in% names(plot_data_p1)) {
        legend_data <- rbind(legend_data, data.table::data.table(
            x = NA_real_, y = NA_real_,
            Legend = factor("XGBoost LOO-CV", levels = c("Treated", "Control/Placebo", "XGBoost LOO-CV"))
        ))
        legend_data[, Legend := factor(Legend, levels = c("Treated", "Control/Placebo", "XGBoost LOO-CV"))]
    }

    p1 <- p1 +
        geom_point(data = legend_data[Legend == "Treated"],
                   aes(x = x, y = y, shape = Legend), color = "#1f78b4", fill = "#1f78b4",
                   size = 2, alpha = 0.8, na.rm = TRUE) +
        geom_point(data = legend_data[Legend == "Control/Placebo"],
                   aes(x = x, y = y, shape = Legend), color = "gray60",
                   size = 2, alpha = 0.3, na.rm = TRUE) +
        {if (show_predictions && "Predicted" %in% names(plot_data_p1))
            geom_point(data = legend_data[Legend == "XGBoost LOO-CV"],
                       aes(x = x, y = y, shape = Legend), color = "red",
                       size = 2.5, stroke = 1.2, alpha = 0.4, na.rm = TRUE)
        } +
        scale_shape_manual(
            name = NULL,
            values = c("Treated" = 19, "Control/Placebo" = 19,
                        "XGBoost LOO-CV" = 4),
            drop = FALSE
        ) +
        guides(
            color = "none",
            fill = "none",
            alpha = "none",
            shape = guide_legend(
                override.aes = list(
                    color = c("#1f78b4", "gray60", "red")[seq_len(nrow(legend_data))],
                    fill = c("#1f78b4", "gray60", NA)[seq_len(nrow(legend_data))],
                    alpha = c(0.8, 0.3, 0.4)[seq_len(nrow(legend_data))],
                    size = c(2, 2, 2.5)[seq_len(nrow(legend_data))],
                    stroke = c(0, 0, 1.2)[seq_len(nrow(legend_data))]
                )
            )
        ) +

        theme_minimal() +
        theme(
            legend.position = "none",
            axis.line.x = element_line(color = "black", linewidth = 0.5),
            axis.line.y = element_blank(),
            axis.text = element_text(colour = "black"),
            axis.title.y = element_text(margin = ggplot2::margin(r = 2))
        ) +
        labs(x = NULL, y = NULL) +
        # Add subtitle explaining predictions if shown
        {if (show_predictions && "Predicted" %in% names(plot_data_p1))
            labs(subtitle = paste0(y_label, " | Red x = XGBoost predictions (LOO-CV)"))
        else
            NULL
        }

    # Calculate specification curve p-values on filtered data
    # This ensures p-values reflect the actual data being plotted
    spec_curve_pvals <- NULL
    if (is.list(long_data)) {
        # Get expected direction from long_data or default to negative
        expected_direction <- if (!is.null(long_data$expected_direction)) long_data$expected_direction else "negative"

        # Get abadie inference if available
        abadie_inference <- if (!is.null(long_data$abadie_inference)) long_data$abadie_inference else NULL

        # Calculate p-values on the filtered data (sc_results_df after all filtering)
        spec_curve_pvals <- calculate_spec_curve_pvalues_filtered(
            filtered_results = sc_results_df,
            abadie_inference = abadie_inference,
            expected_direction = expected_direction,
            name_treated_unit = name_treated_unit
        )
    }

    # Add specification curve p-values annotation if calculated
    if (!is.null(spec_curve_pvals)) {

        # Extract treated unit p-values for both test statistics
        treated_median_pval <- NULL
        treated_stouffer_pval <- NULL

        # Get median tau p-value for treated unit
        if (!is.null(spec_curve_pvals$median_tau_pvalues)) {
            median_pvals <- spec_curve_pvals$median_tau_pvalues
            treated_median_data <- median_pvals[unit_type == "treated"]
            if (nrow(treated_median_data) > 0) {
                treated_median_pval <- treated_median_data$pval_rank_med[1]
            }
        }

        # Get Stouffer's method p-value for treated unit
        if (!is.null(spec_curve_pvals$stouffer_pvalues) && nrow(spec_curve_pvals$stouffer_pvalues) > 0) {
            stouffer_pvals <- spec_curve_pvals$stouffer_pvalues
            treated_stouffer_data <- stouffer_pvals[unit_type == "treated"]
            if (nrow(treated_stouffer_data) > 0) {
                treated_stouffer_pval <- treated_stouffer_data$pval_rank_z[1]
            }
        }

        # Create annotation text if we have at least one p-value
        if (!is.null(treated_median_pval) || !is.null(treated_stouffer_pval)) {
            annotation_lines <- c()

            if (!is.null(treated_median_pval)) {
                treated_median_rank <- treated_median_data$rank_med[1]
                # Calculate total units from the full median_pvals data
                total_median_units <- nrow(median_pvals)
                annotation_lines <- c(annotation_lines, sprintf("Median tau: p = %.3f (rank %d/%d)", treated_median_pval, treated_median_rank, total_median_units))
            }

            if (!is.null(treated_stouffer_pval)) {
                treated_stouffer_rank <- treated_stouffer_data$rank_z[1]
                # Calculate total units from the full stouffer_pvals data
                total_stouffer_units <- nrow(stouffer_pvals)
                annotation_lines <- c(annotation_lines, sprintf("Stouffer: p = %.3f (rank %d/%d)", treated_stouffer_pval, treated_stouffer_rank, total_stouffer_units))
            }

            spec_annotation_text <- paste(annotation_lines, collapse = "\n")

            # Add annotation to upper-right of Panel A with semi-transparent background
            p1 <- p1 +
                annotate("label",
                        x = Inf, y = Inf,
                        label = spec_annotation_text,
                        hjust = 1.05, vjust = 1.2,
                        size = 2.8, color = "#000000",
                        fontface = "plain",
                        fill = "white", alpha = 0.8,
                        linewidth = 0.3)
        }
    }

    # Apply y-axis limits if outlier cropping is requested
    if (!is.null(y_limits)) {
        p1 <- p1 + ylim(y_limits[1], y_limits[2])
    }

    # Create Panel B (Specification Choices and Shapley Values)
    # Panel B: SHAP visualization with optional significance transparency
    has_shap_significance <- "shap_pvalue" %in% names(plot_data_p2) && "shap_significance_level" %in% names(plot_data_p2)

    # Add SHAP values to feature labels if SHAP values are available
    if ("shapley_value" %in% names(plot_data_p2)) {
        # Calculate SHAP summary values for each individual feature for treated unit
        if (shap_label_type == "absolute") {
            feature_shap_summary <- plot_data_p2[`Unit Name` == name_treated_unit & !is.na(shapley_value), 
                                               .(shap_summary = mean(abs(shapley_value), na.rm = TRUE)), 
                                               by = .(feature_group, feature)]
        } else {  # signed
            feature_shap_summary <- plot_data_p2[`Unit Name` == name_treated_unit & !is.na(shapley_value), 
                                               .(shap_summary = mean(shapley_value, na.rm = TRUE)), 
                                               by = .(feature_group, feature)]
        }
        
        # Get all SHAP values for proper color scaling
        all_treated_shap <- plot_data_p2[`Unit Name` == name_treated_unit & !is.na(shapley_value), shapley_value]
        
        # Map SHAP summary values to colors using the same scale as the legend
        feature_shap_summary[, shap_color := map_shap_to_color(shap_summary, shap_label_type, all_treated_shap)]
        
        # Create enhanced feature labels with colored SHAP values in HTML format
        feature_shap_summary[, feature_with_shap := paste0(
            feature, " (<span style='color:", shap_color, "'>", 
            ifelse(shap_summary >= 0, "+", ""), 
            round(shap_summary, 1), "</span>)"
        )]
        
        # Merge back with plot data to add enhanced labels
        plot_data_p2 <- merge(plot_data_p2, 
                             feature_shap_summary[, .(feature_group, feature, feature_with_shap)], 
                             by = c("feature_group", "feature"), 
                             all.x = TRUE)
        
        # Use enhanced labels where available, fallback to original feature names
        plot_data_p2[, feature_display := ifelse(!is.na(feature_with_shap), feature_with_shap, feature)]
        
        # Sort features alphabetically within each feature group for better readability
        plot_data_p2[, feature_display := factor(feature_display, levels = sort(unique(feature_display)))]
        
        # Create color variable that matches the label type
        if (shap_label_type == "absolute") {
            plot_data_p2[, shap_color_value := abs(shapley_value)]
        } else {  # signed
            plot_data_p2[, shap_color_value := shapley_value]
        }
    } else {
        # No SHAP values - use original feature names
        plot_data_p2[, feature_display := feature]
        
        # Sort features alphabetically within each feature group for better readability
        plot_data_p2[, feature_display := factor(feature_display, levels = sort(unique(feature_display)))]
        
        plot_data_p2[, shap_color_value := NA]  # No SHAP coloring
    }

    # --- Build numeric y-position mapping for flat Panel B ---
    # Assign sequential y positions within each feature group, with gaps between groups.
    # Groups are built bottom-to-top (rev of alphabetical) so first group appears at top.
    fg_order <- sort(unique(as.character(plot_data_p2$feature_group)))

    y_pos <- 1.0
    within_group_step <- 0.55  # tighter spacing within categories
    between_group_gap <- 1.8   # wider gap between categories
    group_boundaries <- list()
    y_mapping <- data.table::data.table(
        feature_group = character(0),
        feature_display_chr = character(0),
        y_numeric = numeric(0)
    )

    for (fg in rev(fg_order)) {
        features_in_group <- sort(unique(as.character(
            plot_data_p2[feature_group == fg, feature_display]
        )))
        y_start <- y_pos
        for (i in seq_along(features_in_group)) {
            y_mapping <- rbind(y_mapping, data.table::data.table(
                feature_group = fg,
                feature_display_chr = features_in_group[i],
                y_numeric = y_pos
            ))
            if (i < length(features_in_group)) y_pos <- y_pos + within_group_step
        }
        y_end <- y_pos
        group_boundaries[[fg]] <- list(start = y_start, end = y_end,
                                       mid = (y_start + y_end) / 2)
        y_pos <- y_pos + between_group_gap  # wider gap between groups
    }

    # Merge y_numeric into plot_data_p2
    plot_data_p2[, feature_display_chr := as.character(feature_display)]
    plot_data_p2 <- merge(plot_data_p2,
                          y_mapping[, .(feature_group, feature_display_chr, y_numeric)],
                          by = c("feature_group", "feature_display_chr"),
                          all.x = TRUE)

    # Create ordered factor for discrete y-axis (element_markdown works with discrete scales)
    ordered_levels <- y_mapping[order(y_numeric), feature_display_chr]
    plot_data_p2[, feature_display_ordered := factor(feature_display_chr, levels = ordered_levels)]

    # Compute separator positions midway between adjacent groups (in y_numeric space)
    sorted_mapping <- y_mapping[order(y_numeric)]
    n_levels <- nrow(sorted_mapping)
    separator_positions <- numeric(0)
    for (i in seq_len(n_levels - 1)) {
        if (sorted_mapping$feature_group[i] != sorted_mapping$feature_group[i + 1]) {
            separator_positions <- c(separator_positions,
                                     (sorted_mapping$y_numeric[i] + sorted_mapping$y_numeric[i + 1]) / 2)
        }
    }

    # Group midpoints in y_numeric space (for category label NPC coords)
    group_disc_boundaries <- list()
    for (fg in fg_order) {
        y_vals <- sorted_mapping[feature_group == fg, y_numeric]
        group_disc_boundaries[[fg]] <- list(start = min(y_vals), end = max(y_vals),
                                            mid = (min(y_vals) + max(y_vals)) / 2)
    }

    # Check if SHAP values are available in the data
    if ("shapley_value" %in% names(plot_data_p2)) {
        # Prepare alpha mapping if significance is available
        if (has_shap_significance) {
            plot_data_p2[, shap_alpha_value := pmax(0.2, 1.0 - pmin(shap_pvalue, 1.0))]
        }

        # Create base plot with continuous y-axis (tighter within-group, wider between-group)
        if (has_shap_significance) {
            p2 <- ggplot(plot_data_p2, aes(x = Specification, y = y_numeric)) +
                geom_point(aes(color = shap_color_value, alpha = shap_alpha_value), shape = 15, size = 2.5)
        } else {
            p2 <- ggplot(plot_data_p2, aes(x = Specification, y = y_numeric)) +
                geom_point(aes(color = shap_color_value), shape = 15, size = 2.5)
        }

        # Add color scale based on shap_label_type
        if (shap_label_type == "absolute") {
            p2 <- p2 + scale_color_gradient(
                name = "|SHAP Value| = Change in Treatment Effect",
                low = SHAP_COLORS$light_gray,
                high = SHAP_COLORS$blue,
                guide = guide_colorbar(title.position = "top")
            )
        } else {  # signed
            # Force symmetric limits so the color scale is always centered on 0
            shap_max_abs <- max(abs(plot_data_p2$shap_color_value), na.rm = TRUE)
            p2 <- p2 + scale_color_gradient2(
                name = "SHAP Value = Change in Treatment Effect",
                low = SHAP_COLORS$red,
                mid = SHAP_COLORS$gray,
                high = SHAP_COLORS$blue,
                midpoint = 0,
                limits = c(-shap_max_abs, shap_max_abs),
                guide = guide_colorbar(title.position = "top")
            )
        }

        # Add alpha scale if significance is available
        if (has_shap_significance) {
            p2 <- p2 + scale_alpha_continuous(
                name = "SHAP Significance\n(Solid = Low p-value)",
                range = c(0.2, 1.0),
                breaks = c(0.2, 0.5, 0.95, 0.99, 0.999),
                labels = c("p~1.0", "p~0.5", "p~0.05", "p~0.01", "p~0.0")
            )
        }
    } else {
        # No SHAP values - create basic specification plot
        p2 <- ggplot(plot_data_p2, aes(x = Specification, y = y_numeric)) +
            geom_point(color = "#666666", shape = 15, size = 2.5)
    }

    # Add separator lines between feature groups (discrete positions)
    if (length(separator_positions) > 0) {
        p2 <- p2 + geom_hline(yintercept = separator_positions, color = "gray80", linewidth = 0.3)
    }

    # Add remaining plot elements (flat — no facet_grid)
    # Y-axis text is blank here; richtext labels are added as grobs after conversion.
    p2 <- p2 +
        scale_y_continuous(
            breaks = sorted_mapping$y_numeric,
            labels = sorted_mapping$feature_display_chr,
            expand = expansion(add = 0.5)
        ) +
        theme_minimal() +
        theme(
            axis.line.x = element_line(color = "black", linewidth = 0.5),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = "bottom",
            legend.key.width = grid::unit(1.5, "cm"),
            legend.box = "horizontal"
        ) +
        labs(x = "Specification Number", y = "")


    # --- 3. Combine Plots ---
    # Extract legend from Panel B before removing it for grob assembly
    p2_with_legend <- p2  # preserve for return object
    p2_legend <- cowplot::get_legend(p2)
    p2 <- p2 + theme(legend.position = "none")

    # Extract Panel A legend before suppressing it
    p1_for_legend <- p1 + theme(legend.position = "right")
    p1_legend <- cowplot::get_legend(p1_for_legend)

    # Convert both panels to grobs
    p2_grob <- ggplot2::ggplotGrob(p2)
    p1_grob <- ggplot2::ggplotGrob(p1)

    # --- Build label grobs for Panel B ---
    # Compute NPC coordinates matching Panel B's continuous y-scale
    # (y_numeric range with expansion(add = 0.5))
    y_min <- min(sorted_mapping$y_numeric) - 0.5
    y_max <- max(sorted_mapping$y_numeric) + 0.5
    y_range <- y_max - y_min
    to_npc <- function(y_val) (y_val - y_min) / y_range

    # Category labels (bold feature_group names at group midpoints)
    cat_labels_text <- fg_order
    cat_y_npc <- sapply(fg_order, function(fg) to_npc(group_disc_boundaries[[fg]]$mid))
    cat_grob <- grid::textGrob(
        label = cat_labels_text,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(cat_y_npc, "npc"),
        hjust = 0.5, vjust = 0.5,
        gp = grid::gpar(fontsize = 9, fontface = "bold")
    )

    # Subcategory labels (feature names with colored SHAP values, rendered as richtext)
    # Wrap long labels onto two lines using <br> (richtext supports HTML)
    subcat_y_npc <- to_npc(sorted_mapping$y_numeric)
    subcat_labels <- gsub(" \\+ ", "<br>+ ", sorted_mapping$feature_display_chr)
    subcat_grob <- gridtext::richtext_grob(
        text = subcat_labels,
        x = grid::unit(1, "npc"),
        y = grid::unit(subcat_y_npc, "npc"),
        hjust = 1, vjust = 0.5,
        gp = grid::gpar(fontsize = 8),
        default.units = "npc",
        padding = grid::unit(c(0, 4, 0, 0), "pt"),
        margin = grid::unit(c(0, 2, 0, 0), "pt"),
        box_gp = grid::gpar(col = NA, fill = NA)
    )

    # --- Insert label columns into Panel B grob ---
    p2_panel_idx <- which(p2_grob$layout$name == "panel")
    p2_panel_row_t <- min(p2_grob$layout$t[p2_panel_idx])
    p2_panel_row_b <- max(p2_grob$layout$b[p2_panel_idx])

    # Insert category column (col 1) and subcategory column (col 2) at left
    p2_grob <- gtable::gtable_add_cols(p2_grob, grid::unit(1.5, "cm"), pos = 0)
    p2_grob <- gtable::gtable_add_cols(p2_grob, grid::unit(4.0, "cm"), pos = 1)

    p2_grob <- gtable::gtable_add_grob(p2_grob, cat_grob,
        t = p2_panel_row_t, b = p2_panel_row_b,
        l = 1, r = 1, clip = "off", name = "cat_labels")
    p2_grob <- gtable::gtable_add_grob(p2_grob, subcat_grob,
        t = p2_panel_row_t, b = p2_panel_row_b,
        l = 2, r = 2, clip = "off", name = "subcat_labels")

    # --- Align Panel A columns with Panel B ---
    # Panel B now has one extra column (category labels) that Panel A lacks.
    # Insert dummy columns in Panel A so both panels occupy the same column.
    p2_panel_col <- min(p2_grob$layout$l[p2_grob$layout$name == "panel"])
    p1_panel_col <- min(p1_grob$layout$l[grepl("^panel$", p1_grob$layout$name)])

    n_extra <- p2_panel_col - p1_panel_col
    if (n_extra > 0) {
        for (i in seq_len(n_extra)) {
            p1_grob <- gtable::gtable_add_cols(p1_grob, grid::unit(0, "cm"), pos = 0)
        }
    }

    # Pad right side so both grobs have the same column count
    while (ncol(p1_grob) < ncol(p2_grob)) {
        p1_grob <- gtable::gtable_add_cols(p1_grob, grid::unit(0, "cm"), pos = -1)
    }
    while (ncol(p2_grob) < ncol(p1_grob)) {
        p2_grob <- gtable::gtable_add_cols(p2_grob, grid::unit(0, "cm"), pos = -1)
    }

    # Copy Panel B widths to Panel A — panels now occupy the same column
    p1_grob$widths <- p2_grob$widths

    # Move Panel A's ylab and axis-l into the label columns so they don't
    # consume extra width. Panel B has blank y-axis so those columns are zero-width;
    # we relocate the grobs rather than inflating columns.
    p1_panel_ids <- which(grepl("^panel", p1_grob$layout$name))
    p1_panel_t <- min(p1_grob$layout$t[p1_panel_ids])
    p1_panel_b <- max(p1_grob$layout$b[p1_panel_ids])

    # Move ylab into category column (col 1)
    ylab_ids <- which(grepl("^ylab-l$|^ylab$|^ylab-l-.*|^ylab-.*", p1_grob$layout$name))
    if (length(ylab_ids)) {
        p1_grob$layout$l[ylab_ids] <- 1
        p1_grob$layout$r[ylab_ids] <- 1
    }
    # Move axis-l (tick labels) into subcategory column (col 2), right-aligned
    axisl_ids <- which(grepl("^axis-l", p1_grob$layout$name))
    if (length(axisl_ids)) {
        p1_grob$layout$l[axisl_ids] <- 2
        p1_grob$layout$r[axisl_ids] <- 2
    }

    # Place Panel A legend in the left margin columns (where Panel B has category labels)
    if (!is.null(p1_legend)) {
        p1_grob <- gtable::gtable_add_grob(
            p1_grob, p1_legend,
            t = p1_panel_t, b = p1_panel_b,
            l = 1, r = 2,
            clip = "off", name = "panelA_legend"
        )
    }

    # Combine panels
    final_plot <- gridExtra::gtable_rbind(p1_grob, p2_grob, size = "last")

    # Set Panel A panel row height so it gets ~30% of figure
    p1_panel_row <- unique(p1_grob$layout$t[p1_panel_ids])
    for (r in p1_panel_row) {
        final_plot$heights[r] <- grid::unit(0.7, "null")
    }

    # Add legend row at bottom
    if (!is.null(p2_legend)) {
        final_plot <- gtable::gtable_add_rows(final_plot, grid::unit(2.0, "cm"), pos = -1)
        final_plot <- gtable::gtable_add_grob(final_plot, p2_legend,
            t = nrow(final_plot), l = 1, r = ncol(final_plot),
            clip = "off", name = "legend_bottom")
    }

    # Save plot if file path is provided
    if (!is.na(file_path_save)) {
        pdf(file_path_save, width = width, height = height)
        grid::grid.draw(final_plot)
        dev.off()
    }

    # Prepare comprehensive return object
    return_object <- list(
        final_plot = final_plot,
        panel_a = p1,
        panel_b = p2_with_legend,
        plot_data_p1 = plot_data_p1,
        plot_data_p2 = plot_data_p2,
        computed_shap = computed_shap,
        spec_curve_pvals = spec_curve_pvals,
        filtered_specs = nrow(sc_results_df),
        feature_groups_displayed = groups_to_keep
    )
    
    return(return_object)
}

#' Calculate SHAP Significance via Placebo Inference
#'
#' @title Calculate Statistical Significance of SHAP Values
#' @description Performs placebo-style inference on SHAP values by ranking the treated unit's
#' feature importance against control units for each specification and feature group.
#'
#' @param shap_data Data.table. SHAP values for multiple units with columns: unit, full_spec_id,
#'   feature_group, feature, shapley_value
#' @param treated_unit_name Character. Name of the treated unit
#' @param pvalue_method Character. Method for calculating p-values: "absolute" or "signed".
#'   - "absolute": Ranks based on absolute SHAP values.
#'   - "signed": Ranks based on signed SHAP values, distinguishing positive and negative importance.
#'
#' @return Data.table with columns: full_spec_id, feature_group, shap_pvalue, shap_rank,
#'   total_units, shap_significance_level
#'
#' @details
#' For each specification x feature_group combination:
#' \itemize{
#'   \item Extracts relevant SHAP values for all units
#'   \item Ranks units based on the chosen pvalue_method ("absolute" or "signed")
#'   \item Calculates p-value based on treated unit's rank position
#'   \item Assigns significance levels for visualization mapping
#' }
calculate_shap_significance <- function(shap_data, treated_unit_name, pvalue_method = "absolute") { # NEW PARAMETER

    if (!data.table::is.data.table(shap_data)) {
        shap_data <- data.table::as.data.table(shap_data)
    }

    # Input validation for pvalue_method
    valid_pvalue_methods <- c("absolute", "signed")
    if (!pvalue_method %in% valid_pvalue_methods) {
        stop("pvalue_method must be one of: ", paste(valid_pvalue_methods, collapse = ", "))
    }

    # Calculate significance for each specification x feature_group combination
    significance_results <- shap_data[, {

        treated_shap <- shapley_value[unit == treated_unit_name]

        if (length(treated_shap) == 0) {
            # Treated unit not found for this spec x feature_group
            list(
                shap_pvalue = NA_real_,
                shap_rank = NA_integer_,
                total_units = .N,
                shap_significance_level = NA_real_
            )
        } else {
            if (pvalue_method == "absolute") {
                # Original logic: Rank by absolute SHAP value (most important = rank 1)
                shap_rank_val <- NA_integer_
                abs_shap_values <- abs(shapley_value)
                unit_ranks <- rank(-abs_shap_values, ties.method = "min")
                treated_rank <- unit_ranks[unit == treated_unit_name][1]
                p_value <- treated_rank / .N # One-sided test for unusual importance

            } else if (pvalue_method == "signed") {
                # Signed significance: Compare treated unit's signed SHAP to controls
                # This is more complex as it implies a two-sided test around zero,
                # or separate one-sided tests for positive/negative.
                # A common approach for "signed" significance in placebo inference
                # is to rank within the positive values and within the negative values.
                # For simplicity here, we can consider a 'two-tailed' rank against zero.

                # Rank negative SHAP values (more negative = lower rank)
                # Rank positive SHAP values (more positive = lower rank of absolute)
                # Then combine for a two-sided p-value.

                # Let's use a simpler approach that mimics a t-test idea with ranks:
                # How extreme is the treated unit's value compared to the distribution of controls?

                # Rank all values. If a value is positive and large, its rank should be low.
                # If a value is negative and large, its rank should be low.
                # We need to rank by distance from zero for the magnitude, and then consider sign.

                # For signed:
                # 1. Rank units by their signed SHAP values in ascending order.
                # 2. Get the rank of the treated unit.
                # 3. Calculate p-value based on this rank, assuming a null of zero.
                # This is effectively checking if the treated unit's SHAP value is unusually high OR unusually low.

                # Rank based on the actual signed value (ascending rank: smallest value = 1)
                signed_ranks_asc <- rank(shapley_value, ties.method = "average")
                treated_signed_rank_asc <- signed_ranks_asc[unit == treated_unit_name][1]

                # Rank based on the actual signed value (descending rank: largest value = 1)
                signed_ranks_desc <- rank(-shapley_value, ties.method = "average")
                treated_signed_rank_desc <- signed_ranks_desc[unit == treated_unit_name][1]

                # The p-value is the smaller of (rank_asc / N) or (rank_desc / N)
                # This approximates a two-sided p-value based on ranks.
                p_value <- min(treated_signed_rank_asc / .N, treated_signed_rank_desc / .N) * 2 # Multiply by 2 for two-sided test
                p_value <- pmin(p_value, 1.0) # Cap p-value at 1.0

                # For signed significance, the rank `shap_rank` could be less intuitive
                # as it's not simply "most important". We'll just store the one used for p-value.
                shap_rank_val <- min(treated_signed_rank_asc, treated_signed_rank_desc)

            } else {
                # Should not happen due to validation, but as a safeguard
                p_value <- NA_real_
                shap_rank_val <- NA_integer_
            }

            # Create significance level for visualization mapping
            significance_level <- ifelse(p_value < 0.05, 1.0,
                                         ifelse(p_value < 0.10, 0.6, 0.3))

            list(
                shap_pvalue = p_value,
                shap_rank = shap_rank_val, # Use the most extreme rank for reporting
                total_units = .N,
                shap_significance_level = significance_level
            )
        }

    }, by = .(full_spec_id, feature_group)]

    return(significance_results)
}

#' Calculate Specification Curve P-values on Filtered Data
#'
#' @title Calculate Cross-Specification P-values for Filtered Results
#' @description Calculates specification curve-level p-values by aggregating treatment effects
#' across filtered specifications and ranking units. Implements two test statistics:
#' 1) Median treatment effect across specifications, 2) Average Z-score (Stouffer's method).
#' This function operates on the filtered data that will actually be plotted.
#'
#' @param filtered_results Data.table. Filtered results data (post-filtering by outcomes, RMSE, etc.)
#' @param abadie_inference List. Abadie inference results from original spec_curve output
#' @param expected_direction Character. Expected direction of treatment effect: "negative", "positive", or "two_sided"
#' @param name_treated_unit Character. Name of the treated unit for filtering inference results
#'
#' @return List containing specification curve p-values:
#' \itemize{
#'   \item median_tau_pvalues: P-values based on median treatment effect across filtered specifications
#'   \item stouffer_pvalues: P-values based on average Z-score across filtered specifications (Stouffer's method)
#' }
#'
#' @details
#' This function calculates p-values on the subset of data that will actually be visualized,
#' ensuring that the statistical inference matches the displayed results. Key features:
#' \itemize{
#'   \item Works on filtered data (by RMSE threshold, outcomes, etc.)
#'   \item Calculates median treatment effect for each unit across filtered specifications only
#'   \item Converts individual specification p-values to Z-scores and averages them
#'   \item Ranks all units based on filtered results and assigns p-values
#'   \item Accounts for expected direction when ranking (negative vs positive effects)
#' }
calculate_spec_curve_pvalues_filtered <- function(filtered_results, abadie_inference = NULL, expected_direction = "negative", name_treated_unit = NULL) {

    if (nrow(filtered_results) == 0) {
        return(list(
            median_tau_pvalues = data.table(),
            stouffer_pvalues = data.table()
        ))
    }

    # Test Statistic 1: Median Treatment Effect Across Filtered Specifications
    # Calculate average tau for each unit across filtered specs (post-period only)
    avg_tau_by_spec <- filtered_results[post_period == TRUE, .(
        ave_tau = mean(tau)
    ), by = .(full_spec_id, unit_name, unit_type)]

    # Get median of those averages for each unit across filtered specifications
    median_tau_by_unit <- avg_tau_by_spec[, .(
        median_tau = median(ave_tau),
        n_specs = .N
    ), by = .(unit_name, unit_type)]

    # Rank units based on expected direction and assign p-values
    if (expected_direction == "negative") {
        # Most negative gets rank 1 (lowest p-value)
        median_tau_by_unit[, rank_med := rank(median_tau, ties.method = "min")]
    } else if (expected_direction == "positive") {
        # Most positive gets rank 1 (lowest p-value)
        median_tau_by_unit[, rank_med := rank(-median_tau, ties.method = "min")]
    } else {
        # Two-sided: rank by absolute value (most extreme gets rank 1)
        median_tau_by_unit[, rank_med := rank(-abs(median_tau), ties.method = "min")]
    }

    # Calculate p-values
    median_tau_by_unit[, pval_rank_med := rank_med / .N]

    # Test Statistic 2: Average Z-score (Stouffer's Method)
    # Filter Abadie inference to match filtered specifications
    stouffer_results <- data.table()

    if (!is.null(abadie_inference) && "p_values_rmse_ratio" %in% names(abadie_inference)) {

        abadie_pvals <- as.data.table(abadie_inference$p_values_rmse_ratio)

        # Get the specifications that remain after filtering
        filtered_spec_ids <- unique(filtered_results$full_spec_id)

        # Filter Abadie p-values to match filtered specifications
        abadie_pvals_filtered <- abadie_pvals[full_spec_id %in% filtered_spec_ids]

        if (nrow(abadie_pvals_filtered) > 0) {
            # Merge treatment effects with filtered p-values
            tau_pval_merged <- merge(
                avg_tau_by_spec,
                abadie_pvals_filtered[, .(full_spec_id, unit_name, p_value_one_sided)],
                by = c('full_spec_id', 'unit_name'),
                all.x = TRUE
            )

            # Convert p-values to Z-scores
            tau_pval_merged[, zscore := qnorm(p_value_one_sided, mean = 0, sd = 1, lower.tail = TRUE)]

            # Handle extreme p-values
            tau_pval_merged[p_value_one_sided == 1, zscore := 4.25]
            tau_pval_merged[p_value_one_sided == 0, zscore := -4.25]
            tau_pval_merged[is.na(zscore), zscore := 0]    # Handle any remaining NAs

            # Calculate mean Z-score for each unit across filtered specifications
            mean_z_by_unit <- tau_pval_merged[, .(
                mean_z = mean(zscore, na.rm = TRUE),
                n_specs_with_pvals = sum(!is.na(zscore))
            ), by = .(unit_name, unit_type)]

            # Only proceed if we have units with p-values
            if (nrow(mean_z_by_unit) > 0) {
                # Rank units by mean Z-score based on expected direction
                if (expected_direction == "negative") {
                    # Most negative Z-score gets rank 1 (most significant in negative direction)
                    mean_z_by_unit[, rank_z := rank(mean_z, ties.method = "min")]
                } else if (expected_direction == "positive") {
                    # Most positive Z-score gets rank 1 (most significant in positive direction)
                    mean_z_by_unit[, rank_z := rank(-mean_z, ties.method = "min")]
                } else {
                    # Two-sided: most extreme (absolute) Z-score gets rank 1
                    mean_z_by_unit[, rank_z := rank(-abs(mean_z), ties.method = "min")]
                }

                # Calculate p-values
                mean_z_by_unit[, pval_rank_z := rank_z / .N]

                stouffer_results <- mean_z_by_unit
            }
        }
    }

    # Return both test statistics
    return(list(
        median_tau_pvalues = median_tau_by_unit,
        stouffer_pvalues = stouffer_results
    ))
}

