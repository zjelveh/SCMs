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
#' @param curve_stat Character vector. Curve-level statistics to compute for placebo-in-space
#'   inference. Options:
#'   \itemize{
#'     \item \code{"abs_median_tau"}: \eqn{|\mathrm{median}_s(\tau_s)|}
#'     \item \code{"median_abs_tau"}: \eqn{\mathrm{median}_s(|\tau_s|)}
#'     \item \code{"consistency_ratio"}: \eqn{|\sum_s \tau_s| / \sum_s |\tau_s|}
#'   }
#'   Default is \code{"abs_median_tau"}.
#' @param weighting Character vector. Weighting modes to evaluate for curve-level statistics.
#'   Options: \code{"none"} and \code{"pre_rmspe_percentile"}.
#'   Default is \code{c("none", "pre_rmspe_percentile")}.
#'   Weighting is applied generically to each requested curve statistic.
#' @param two_sided Logical. If \code{TRUE} (default), curve-level p-values use absolute-tail
#'   placebo comparison: \code{(1 + sum(abs(placebo) >= abs(observed))) / (J + 1)}.
#'   If \code{FALSE}, one-sided comparison is used with \code{expected_direction} from
#'   \code{long_data} (\code{"negative"} or \code{"positive"} required).
#' @param grid_policy Character. How to handle incomplete placebo specification grids when
#'   computing curve-level p-values. Options:
#'   \itemize{
#'     \item \code{"strict"} (default): fail if any placebo unit is missing/has extra specs.
#'     \item \code{"drop_incomplete_units"}: keep full treated spec grid, drop placebo units
#'       that do not exactly match it.
#'     \item \code{"intersect_specs"}: keep all units, restrict to specs shared by all units.
#'   }
#' @param min_placebos Integer. Minimum number of placebo units required after applying
#'   \code{grid_policy}. Default is 1.
#' @param min_specs Integer. Minimum number of specifications required after applying
#'   \code{grid_policy}. Default is 2.
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
#' @param shap_label_type Character. How to display SHAP values in y-axis labels: "absolute" or "signed".
#'   - "absolute": Shows mean |SHAP| values, e.g., "ridge (0.123)"
#'   - "signed": Shows mean SHAP values with sign, e.g., "ridge (+0.089)" or "lasso (-0.045)"
#'   Default is "absolute".
#' @param richtext_feature_labels Logical. Whether to render Panel B feature labels
#'   using rich HTML text (colored SHAP values). Set to FALSE for plain text labels,
#'   which can be more robust in some SVG renderers (e.g., GitHub). Default is TRUE.
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
#'   \item plot_data_p2: data.table with Panel B plotting data including columns: Specification, feature_group, feature, shapley_value (if SHAP computed)
#'   \item computed_shap: Complete results from internal SHAP computation (NULL if external shap_values provided or show_shap=FALSE). 
#'     Contains: results (feature importance), shapley (SHAP values), predictions (model predictions), models (trained XGBoost models), config (SHAP configuration)
#'   \item spec_curve_pvals: List with specification curve-level p-values calculated on filtered data:
#'     \code{treated_summary} (treated-unit curve estimates and p-values) and
#'     \code{stats_by_unit} (curve statistics for treated and placebo units), plus
#'     \code{attrition_report} summarizing unit/spec retention under \code{grid_policy}.
#'     NULL if no inference data available.
#'   \item filtered_specs: Integer count of specifications remaining after filtering (before filtering if filter_specs=NULL)
#'   \item feature_groups_displayed: Character vector of feature groups shown in Panel B (only groups with variation in filtered data)
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' specs <- data.table::data.table(
#'   full_spec_id = paste0("s", 1:4),
#'   outcome_model = c("none", "ridge", "none", "ridge"),
#'   const = c("simplex", "lasso", "simplex", "lasso"),
#'   fw = c("uniform", "uniform", "optimize", "optimize"),
#'   feat = c("f1", "f1", "f2", "f2")
#' )
#' long_data <- data.table::CJ(
#'   unit_name = c("treated", "donor1"),
#'   full_spec_id = as.character(specs[["full_spec_id"]]),
#'   post_period = c(FALSE, TRUE),
#'   sorted = FALSE
#' )
#' long_data[, unit_type := ifelse(unit_name == "treated", "treated", "control")]
#' long_data <- merge(long_data, specs, by = "full_spec_id", all.x = TRUE)
#' long_data[, outcome := "y"]
#' long_data[, rmse := 0.2 + as.numeric(factor(full_spec_id)) * 0.01]
#' long_data[, tau := ifelse(
#'   unit_name == "treated",
#'   as.numeric(factor(full_spec_id)) * 0.5 + ifelse(post_period, 0.2, 0),
#'   as.numeric(factor(full_spec_id)) * 0.1 + ifelse(post_period, 0.05, 0)
#' )]
#' out <- plot_spec_curve(
#'   long_data = long_data,
#'   name_treated_unit = "treated",
#'   show_shap = FALSE,
#'   show_pvalues = FALSE,
#'   curve_stat = "abs_median_tau",
#'   weighting = "none"
#' )
#' out$panel_a
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
    curve_stat = "abs_median_tau",
    weighting = c("none", "pre_rmspe_percentile"),
    two_sided = TRUE,
    grid_policy = "strict",
    min_placebos = 1L,
    min_specs = 2L,
    test_statistic = "rmse_ratio",
    null_distribution = "placebo",
    crop_outliers = "none",
    sort_by = "tau",
    filter_specs = NULL,
    show_shap = TRUE,
    shap_label_type = "absolute",
    richtext_feature_labels = TRUE,
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
    
    # Validate shap_label_type
    valid_shap_label_types <- c("absolute", "signed")
    if (!shap_label_type %in% valid_shap_label_types) {
        stop("shap_label_type must be one of: ", paste(valid_shap_label_types, collapse = ", "))
    }

    # Validate curve-level inference settings
    valid_curve_stats <- c("abs_median_tau", "median_abs_tau", "consistency_ratio")
    if (!is.character(curve_stat) || length(curve_stat) < 1) {
        stop("curve_stat must be a non-empty character vector.")
    }
    curve_stat <- unique(curve_stat)
    invalid_curve_stats <- setdiff(curve_stat, valid_curve_stats)
    if (length(invalid_curve_stats) > 0) {
        stop("curve_stat contains invalid values: ", paste(invalid_curve_stats, collapse = ", "),
             ". Valid options are: ", paste(valid_curve_stats, collapse = ", "))
    }

    valid_weighting <- c("none", "pre_rmspe_percentile")
    if (!is.character(weighting) || length(weighting) < 1) {
        stop("weighting must be a non-empty character vector.")
    }
    weighting <- unique(weighting)
    invalid_weighting <- setdiff(weighting, valid_weighting)
    if (length(invalid_weighting) > 0) {
        stop("weighting contains invalid values: ", paste(invalid_weighting, collapse = ", "),
             ". Valid options are: ", paste(valid_weighting, collapse = ", "))
    }
    if (!is.logical(two_sided) || length(two_sided) != 1 || is.na(two_sided)) {
        stop("two_sided must be a single TRUE/FALSE value.")
    }
    valid_grid_policy <- c("strict", "drop_incomplete_units", "intersect_specs")
    if (!is.character(grid_policy) || length(grid_policy) != 1 || !grid_policy %in% valid_grid_policy) {
        stop("grid_policy must be one of: ", paste(valid_grid_policy, collapse = ", "))
    }
    if (!is.numeric(min_placebos) || length(min_placebos) != 1 || is.na(min_placebos) ||
        is.infinite(min_placebos) || min_placebos < 1 || min_placebos %% 1 != 0) {
        stop("min_placebos must be a single integer >= 1.")
    }
    if (!is.numeric(min_specs) || length(min_specs) != 1 || is.na(min_specs) ||
        is.infinite(min_specs) || min_specs < 2 || min_specs %% 1 != 0) {
        stop("min_specs must be a single integer >= 2.")
    }
    min_placebos <- as.integer(min_placebos)
    min_specs <- as.integer(min_specs)

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
                   size = 2, alpha = 0.8, na.rm = TRUE, inherit.aes = FALSE) +
        geom_point(data = legend_data[Legend == "Control/Placebo"],
                   aes(x = x, y = y, shape = Legend), color = "gray60",
                   size = 2, alpha = 0.3, na.rm = TRUE, inherit.aes = FALSE) +
        {if (show_predictions && "Predicted" %in% names(plot_data_p1))
            geom_point(data = legend_data[Legend == "XGBoost LOO-CV"],
                       aes(x = x, y = y, shape = Legend), color = "red",
                       size = 2.5, stroke = 1.2, alpha = 0.4, na.rm = TRUE, inherit.aes = FALSE)
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
    # and avoids unnecessary inference work when p-values are hidden.
    spec_curve_pvals <- NULL
    if (show_pvalues && is.list(long_data)) {
        expected_direction <- if (!is.null(long_data$expected_direction)) long_data$expected_direction else NULL

        if (!two_sided) {
            if (is.null(expected_direction) || !expected_direction %in% c("negative", "positive")) {
                stop("two_sided=FALSE requires long_data$expected_direction to be either 'negative' or 'positive'.")
            }
        } else if (is.null(expected_direction)) {
            expected_direction <- "two_sided"
        }

        # Calculate p-values on the filtered data (sc_results_df after all filtering)
        spec_curve_pvals <- calculate_spec_curve_pvalues_filtered(
            filtered_results = sc_results_df,
            curve_stat = curve_stat,
            weighting = weighting,
            two_sided = two_sided,
            expected_direction = expected_direction,
            grid_policy = grid_policy,
            min_placebos = min_placebos,
            min_specs = min_specs
        )
    }

    # Add specification curve p-values annotation if calculated
    if (!is.null(spec_curve_pvals)) {
        treated_curve <- spec_curve_pvals$treated_summary

        if (!is.null(treated_curve) && nrow(treated_curve) > 0) {
            req_cols <- c("curve_statistic", "weighting", "p_value", "n_extreme", "n_placebos")
            missing_cols <- setdiff(req_cols, names(treated_curve))
            if (length(missing_cols) > 0) {
                stop(
                    "treated_summary is missing required columns for curve-level annotation: ",
                    paste(missing_cols, collapse = ", ")
                )
            }

            treated_curve <- data.table::copy(treated_curve)
            stat_levels <- c("abs_median_tau", "median_abs_tau", "consistency_ratio")
            treated_curve[, stat_rank := match(curve_statistic, stat_levels)]
            if (any(is.na(treated_curve$stat_rank))) {
                bad_stats <- unique(treated_curve[is.na(stat_rank), curve_statistic])
                stop("Unknown curve_statistic value(s) in treated_summary: ", paste(bad_stats, collapse = ", "))
            }
            treated_curve[, weight_rank := ifelse(weighting == "none", 1L, 2L)]
            data.table::setorder(treated_curve, stat_rank, weight_rank)

            format_rank_value <- function(row_dt) {
                rank_num <- as.integer(row_dt$n_extreme) + 1L
                rank_den <- as.integer(row_dt$n_placebos) + 1L
                sprintf(
                    "%d/%d (p = %.3f)",
                    rank_num,
                    rank_den,
                    row_dt$p_value
                )
            }

            annotation_lines <- c()

            median_rows <- treated_curve[curve_statistic == "abs_median_tau"]
            if (nrow(median_rows) > 0) {
                row_none <- median_rows[weighting == "none"]
                row_wt <- median_rows[weighting == "pre_rmspe_percentile"]
                if (nrow(row_none) > 0 && nrow(row_wt) > 0) {
                    annotation_lines <- c(
                        annotation_lines,
                        sprintf(
                            "|median(tau)|: %s (wt: %s)",
                            format_rank_value(row_none[1]),
                            format_rank_value(row_wt[1])
                        )
                    )
                } else {
                    row_one <- median_rows[1]
                    weight_tag <- if (row_one$weighting == "pre_rmspe_percentile") " (wt)" else ""
                    annotation_lines <- c(
                        annotation_lines,
                        sprintf(
                            "|median(tau)|%s: %s",
                            weight_tag,
                            format_rank_value(row_one)
                        )
                    )
                }
            }

            if (length(annotation_lines) > 0) {
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
    }

    # Apply y-axis limits if outlier cropping is requested
    if (!is.null(y_limits)) {
        p1 <- p1 + ylim(y_limits[1], y_limits[2])
    }

    # Create Panel B (Specification Choices and Shapley Values)
    # Panel B: SHAP visualization

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
        
        # Create feature labels with SHAP summaries.
        if (isTRUE(richtext_feature_labels)) {
            # HTML labels preserve SHAP sign color in rendered text.
            feature_shap_summary[, feature_with_shap := paste0(
                feature, " (<span style='color:", shap_color, "'>",
                ifelse(shap_summary >= 0, "+", ""),
                round(shap_summary, 1), "</span>)"
            )]
        } else {
            # Plain-text labels are more robust for SVG display on GitHub.
            feature_shap_summary[, feature_with_shap := paste0(
                feature, " (",
                ifelse(shap_summary >= 0, "+", ""),
                round(shap_summary, 1),
                ")"
            )]
        }
        
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
    if (isTRUE(richtext_feature_labels)) {
        within_group_step <- 0.70  # richer labels can sit closer with HTML wrapping
        between_group_gap <- 1.65  # moderate inter-group spacing
    } else {
        # Plain-text SVG labels render as multi-line text; use larger row spacing.
        within_group_step <- 1.15
        between_group_gap <- 2.05  # moderate inter-group spacing
    }
    group_boundaries <- list()
    y_mapping <- data.table::data.table(
        feature_group = character(0),
        feature_display_chr = character(0),
        y_numeric = numeric(0)
    )

    for (fg in rev(fg_order)) {
        # Reverse within-group order relative to prior behavior so top-to-bottom
        # is the opposite alphabetical direction (e.g., optimize above uniform).
        features_in_group <- sort(unique(as.character(
            plot_data_p2[feature_group == fg, feature_display]
        )), decreasing = TRUE)
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
        # Create base plot with continuous y-axis (tighter within-group, wider between-group)
        p2 <- ggplot(plot_data_p2, aes(x = Specification, y = y_numeric)) +
            geom_point(aes(color = shap_color_value), shape = 15, size = 2.5)

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
        x = grid::unit(0.56, "npc"),
        y = grid::unit(cat_y_npc, "npc"),
        hjust = 0.5, vjust = 0.5,
        gp = grid::gpar(fontsize = 9, fontface = "bold")
    )

    # Subcategory labels (feature names + SHAP summary).
    subcat_y_npc <- to_npc(sorted_mapping$y_numeric)
    if (isTRUE(richtext_feature_labels)) {
        # Wrap labels onto two lines for richtext rendering.
        subcat_labels <- gsub(" \\+ ", "<br>+ ", sorted_mapping$feature_display_chr)
        subcat_grob <- gridtext::richtext_grob(
            text = subcat_labels,
            x = grid::unit(1, "npc"),
            y = grid::unit(subcat_y_npc, "npc"),
            hjust = 1, vjust = 0.5,
            gp = grid::gpar(fontsize = 8, lineheight = 1.1),
            default.units = "npc",
            padding = grid::unit(c(1, 8, 1, 1), "pt"),
            margin = grid::unit(c(1, 3, 1, 1), "pt"),
            box_gp = grid::gpar(col = NA, fill = NA)
        )
    } else {
        # Plain text fallback avoids HTML rendering issues in some SVG viewers.
        subcat_labels <- gsub(" \\(", "\n(", sorted_mapping$feature_display_chr)
        subcat_grob <- grid::textGrob(
            label = subcat_labels,
            x = grid::unit(1, "npc"),
            y = grid::unit(subcat_y_npc, "npc"),
            hjust = 1, vjust = 0.5,
            gp = grid::gpar(fontsize = 8, lineheight = 1.1)
        )
    }

    # --- Insert label columns into Panel B grob ---
    p2_panel_idx <- which(p2_grob$layout$name == "panel")
    p2_panel_row_t <- min(p2_grob$layout$t[p2_panel_idx])
    p2_panel_row_b <- max(p2_grob$layout$b[p2_panel_idx])

    # Insert category column (col 1) and subcategory column (col 2) at left
    p2_grob <- gtable::gtable_add_cols(p2_grob, grid::unit(1.9, "cm"), pos = 0)
    p2_grob <- gtable::gtable_add_cols(p2_grob, grid::unit(4.8, "cm"), pos = 1)

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
    # Add a small global left gutter so category labels are not flush with image edge.
    final_plot <- gtable::gtable_add_cols(final_plot, grid::unit(0.45, "cm"), pos = 0)

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

# Validate numeric vectors for curve-level statistics (fail hard).
validate_curve_numeric_vector <- function(x, arg_name, min_length = 2L) {
    if (!is.numeric(x)) {
        stop(arg_name, " must be numeric.")
    }
    if (length(x) < min_length) {
        stop(arg_name, " must have at least ", min_length, " elements.")
    }
    bad_idx <- which(is.na(x) | is.nan(x) | is.infinite(x))
    if (length(bad_idx) > 0) {
        idx_msg <- paste(utils::head(bad_idx, 10), collapse = ", ")
        stop(arg_name, " contains NA/NaN/Inf at indices: ", idx_msg)
    }
}

# Midrank percentile transform in (0,1).
percentile_rank <- function(x, arg_name = "x") {
    validate_curve_numeric_vector(x, arg_name, min_length = 2L)
    ranks <- base::rank(x, ties.method = "average")
    percentiles <- ranks / (length(x) + 1)
    if (any(!(percentiles > 0 & percentiles < 1))) {
        stop("percentile_rank produced values outside (0,1) for ", arg_name)
    }
    percentiles
}

validate_curve_weights <- function(weights, tau_length, arg_name = "weights") {
    if (is.null(weights)) {
        return(invisible(NULL))
    }
    if (!is.numeric(weights)) {
        stop(arg_name, " must be numeric.")
    }
    if (length(weights) != tau_length) {
        stop(arg_name, " must have length ", tau_length, ".")
    }
    bad_idx <- which(is.na(weights) | is.nan(weights) | is.infinite(weights))
    if (length(bad_idx) > 0) {
        idx_msg <- paste(utils::head(bad_idx, 10), collapse = ", ")
        stop(arg_name, " contains NA/NaN/Inf at indices: ", idx_msg)
    }
    if (any(weights <= 0)) {
        stop(arg_name, " must be strictly positive.")
    }
}

compute_pre_rmspe_percentile_weights <- function(pre_rmspe) {
    validate_curve_numeric_vector(pre_rmspe, "pre_rmspe", min_length = 2L)
    pw <- percentile_rank(pre_rmspe, arg_name = "pre_rmspe")
    w <- 1 - pw
    validate_curve_weights(w, length(pre_rmspe), arg_name = "pre_rmspe_weights")
    w_sum <- sum(w)
    if (!is.finite(w_sum) || w_sum <= 0) {
        stop("sum(pre_rmspe_weights) must be positive and finite.")
    }
    w
}

curve_stat_median <- function(tau, weights = NULL) {
    validate_curve_numeric_vector(tau, "tau", min_length = 2L)
    if (is.null(weights)) {
        return(stats::median(tau))
    }

    validate_curve_weights(weights, length(tau), arg_name = "weights")
    o <- order(tau)
    tau_sorted <- tau[o]
    w_sorted <- weights[o]
    w_sum <- sum(w_sorted)
    if (!is.finite(w_sum) || w_sum <= 0) {
        stop("sum(weights) must be positive and finite for weighted median.")
    }

    cdf <- cumsum(w_sorted) / w_sum
    idx <- which(cdf >= 0.5)[1]
    if (is.na(idx)) {
        stop("Failed to locate weighted median index.")
    }
    if (cdf[idx] > 0.5 || idx == length(tau_sorted)) {
        return(tau_sorted[idx])
    }
    # If the CDF lands exactly on 0.5, average adjacent points (median-compatible tie handling).
    mean(c(tau_sorted[idx], tau_sorted[idx + 1L]))
}

curve_stat_abs_median_tau <- function(tau, weights = NULL) {
    med <- curve_stat_median(tau, weights = weights)
    out <- abs(med)
    if (!is.finite(out)) {
        stop("curve_stat_abs_median_tau produced a non-finite value.")
    }
    out
}

curve_stat_median_abs_tau <- function(tau, weights = NULL) {
    validate_curve_numeric_vector(tau, "tau", min_length = 2L)
    tau_abs <- abs(tau)
    out <- curve_stat_median(tau_abs, weights = weights)
    if (!is.finite(out)) {
        stop("curve_stat_median_abs_tau produced a non-finite value.")
    }
    out
}

curve_stat_consistency_ratio <- function(tau, weights = NULL) {
    validate_curve_numeric_vector(tau, "tau", min_length = 2L)
    validate_curve_weights(weights, length(tau), arg_name = "weights")

    if (is.null(weights)) {
        numer <- abs(sum(tau))
        denom <- sum(abs(tau))
    } else {
        w_sum <- sum(weights)
        if (!is.finite(w_sum) || w_sum <= 0) {
            stop("sum(weights) must be positive and finite for weighted consistency_ratio.")
        }
        numer <- abs(sum(weights * tau))
        denom <- sum(weights * abs(tau))
    }

    if (!is.finite(denom) || denom <= 0) {
        stop("consistency_ratio denominator must be positive and finite.")
    }
    out <- numer / denom
    if (!is.finite(out)) {
        stop("curve_stat_consistency_ratio produced a non-finite value.")
    }
    if (out < -1e-12 || out > 1 + 1e-12) {
        stop("curve_stat_consistency_ratio produced value outside [0,1]: ", out)
    }
    min(max(out, 0), 1)
}

apply_curve_statistic <- function(curve_statistic, tau, weights = NULL) {
    if (!is.character(curve_statistic) || length(curve_statistic) != 1) {
        stop("curve_statistic must be a single character value.")
    }
    if (curve_statistic == "abs_median_tau") {
        return(curve_stat_abs_median_tau(tau, weights = weights))
    }
    if (curve_statistic == "median_abs_tau") {
        return(curve_stat_median_abs_tau(tau, weights = weights))
    }
    if (curve_statistic == "consistency_ratio") {
        return(curve_stat_consistency_ratio(tau, weights = weights))
    }
    stop("Unhandled curve statistic: ", curve_statistic)
}

compute_placebo_curve_pvalue <- function(observed_stat, placebo_stats, two_sided = TRUE,
                                         expected_direction = "two_sided", stat_name = "curve_stat") {
    if (!is.numeric(observed_stat) || length(observed_stat) != 1 || is.na(observed_stat) ||
        is.nan(observed_stat) || is.infinite(observed_stat)) {
        stop("Observed ", stat_name, " must be a finite numeric scalar.")
    }
    validate_curve_numeric_vector(placebo_stats, paste0("placebo_", stat_name), min_length = 1L)
    j <- length(placebo_stats)

    if (two_sided) {
        n_extreme <- sum(abs(placebo_stats) >= abs(observed_stat))
    } else if (expected_direction == "negative") {
        n_extreme <- sum(placebo_stats <= observed_stat)
    } else if (expected_direction == "positive") {
        n_extreme <- sum(placebo_stats >= observed_stat)
    } else {
        stop("For one-sided curve p-values, expected_direction must be 'negative' or 'positive'.")
    }

    list(
        p_value = (1 + n_extreme) / (j + 1),
        n_extreme = n_extreme,
        n_placebos = j
    )
}

#' Calculate Specification Curve P-values on Filtered Data
#'
#' @title Calculate Placebo-in-Space Curve-Level Inference
#' @description Calculates curve-level statistics and p-values on filtered specification-curve
#' results using placebo-in-space comparison.
#'
#' @param filtered_results Data.table. Filtered results data (post-filtering by outcomes, RMSE, etc.).
#' @param curve_stat Character vector. Curve-level statistics to compute:
#'   \code{"abs_median_tau"}, \code{"median_abs_tau"}, \code{"consistency_ratio"}.
#' @param weighting Character vector. Weighting mode(s) applied to each requested
#'   curve-level statistic: \code{"none"} or \code{"pre_rmspe_percentile"}.
#' @param two_sided Logical. Whether to compute two-sided p-values using absolute-value tails.
#' @param expected_direction Character. Expected sign direction for one-sided tests
#'   (\code{"negative"} or \code{"positive"}). Ignored when \code{two_sided = TRUE}.
#' @param grid_policy Character. Grid mismatch policy:
#'   \code{"strict"}, \code{"drop_incomplete_units"}, \code{"intersect_specs"}.
#' @param min_placebos Integer. Minimum number of placebo units required after applying
#'   \code{grid_policy}. Must be >= 1.
#' @param min_specs Integer. Minimum number of specifications required after applying
#'   \code{grid_policy}. Must be >= 2.
#'
#' @return List with:
#' \itemize{
#'   \item \code{treated_summary}: Treated-unit curve statistics and p-values.
#'   \item \code{stats_by_unit}: Curve statistics for treated and placebo units.
#'   \item \code{attrition_report}: Unit/spec retention summary under \code{grid_policy}.
#' }
#'
#' @details
#' For each unit and specification, this function uses the signed per-spec effect summary
#' \code{tau_s}. Two-sided p-values always compare absolute tails against placebo-unit
#' curve statistics:
#' \deqn{p = (1 + \#\{|T_j| \ge |T_{treated}|\})/(J+1).}
#' The supported curve statistics are:
#' \itemize{
#'   \item \code{abs_median_tau}: \eqn{|\mathrm{median}_s(\tau_s)|}
#'   \item \code{median_abs_tau}: \eqn{\mathrm{median}_s(|\tau_s|)}
#'   \item \code{consistency_ratio}: \eqn{|\sum_s \tau_s| / \sum_s |\tau_s|}
#' }
#' For weighted variants, pre-RMSPE-percentile weights are defined within unit as
#' \eqn{w_s = 1 - rank(\mathrm{RMSPE}^{pre}_s)/(S+1)} and applied directly to the
#' statistic computation.
calculate_spec_curve_pvalues_filtered <- function(
    filtered_results,
    curve_stat = c("abs_median_tau", "median_abs_tau", "consistency_ratio"),
    weighting = c("none", "pre_rmspe_percentile"),
    two_sided = TRUE,
    expected_direction = "two_sided",
    grid_policy = "strict",
    min_placebos = 1L,
    min_specs = 2L
) {
    if (!data.table::is.data.table(filtered_results)) {
        filtered_results <- data.table::as.data.table(filtered_results)
    }
    if (nrow(filtered_results) == 0) {
        stop("filtered_results is empty. Cannot compute curve-level inference.")
    }

    valid_curve_stats <- c("abs_median_tau", "median_abs_tau", "consistency_ratio")
    curve_stat <- unique(curve_stat)
    invalid_curve_stats <- setdiff(curve_stat, valid_curve_stats)
    if (length(invalid_curve_stats) > 0) {
        stop("Invalid curve_stat values: ", paste(invalid_curve_stats, collapse = ", "))
    }
    if (!is.character(weighting) || length(weighting) < 1) {
        stop("weighting must be a non-empty character vector.")
    }
    weighting <- unique(weighting)
    valid_weighting <- c("none", "pre_rmspe_percentile")
    invalid_weighting <- setdiff(weighting, valid_weighting)
    if (length(invalid_weighting) > 0) {
        stop("Invalid weighting values: ", paste(invalid_weighting, collapse = ", "))
    }
    if (!is.logical(two_sided) || length(two_sided) != 1 || is.na(two_sided)) {
        stop("two_sided must be a single TRUE/FALSE value.")
    }
    if (!two_sided && !expected_direction %in% c("negative", "positive")) {
        stop("For one-sided curve p-values, expected_direction must be 'negative' or 'positive'.")
    }
    valid_grid_policy <- c("strict", "drop_incomplete_units", "intersect_specs")
    if (!is.character(grid_policy) || length(grid_policy) != 1 || !grid_policy %in% valid_grid_policy) {
        stop("grid_policy must be one of: ", paste(valid_grid_policy, collapse = ", "))
    }
    if (!is.numeric(min_placebos) || length(min_placebos) != 1 || is.na(min_placebos) ||
        is.infinite(min_placebos) || min_placebos < 1 || min_placebos %% 1 != 0) {
        stop("min_placebos must be a single integer >= 1.")
    }
    if (!is.numeric(min_specs) || length(min_specs) != 1 || is.na(min_specs) ||
        is.infinite(min_specs) || min_specs < 2 || min_specs %% 1 != 0) {
        stop("min_specs must be a single integer >= 2.")
    }
    min_placebos <- as.integer(min_placebos)
    min_specs <- as.integer(min_specs)

    required_cols <- c("full_spec_id", "unit_name", "unit_type", "post_period", "tau")
    missing_cols <- setdiff(required_cols, names(filtered_results))
    if (length(missing_cols) > 0) {
        stop("filtered_results is missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    analysis_dt <- filtered_results[unit_type %in% c("treated", "control")]
    if (nrow(analysis_dt) == 0) {
        stop("No treated/control rows available after filtering.")
    }

    tau_by_spec <- analysis_dt[post_period == TRUE, .(
        tau_s = mean(tau)
    ), by = .(full_spec_id, unit_name, unit_type)]
    if (nrow(tau_by_spec) == 0) {
        stop("No post-period rows available for curve-level inference.")
    }
    bad_tau <- tau_by_spec[is.na(tau_s) | is.nan(tau_s) | is.infinite(tau_s)]
    if (nrow(bad_tau) > 0) {
        bad_msg <- paste(utils::head(
            paste0(bad_tau$unit_name, "@", bad_tau$full_spec_id), 10
        ), collapse = ", ")
        stop("Non-finite tau_s detected for: ", bad_msg)
    }

    treated_units <- unique(tau_by_spec[unit_type == "treated", unit_name])
    if (length(treated_units) != 1) {
        stop("Expected exactly one treated unit in filtered results; found ", length(treated_units), ".")
    }
    treated_unit <- treated_units[[1]]

    placebo_units <- unique(tau_by_spec[unit_type == "control", unit_name])
    if (length(placebo_units) < 1) {
        stop("At least one placebo/control unit is required for curve-level inference.")
    }

    treated_spec_ids <- sort(unique(tau_by_spec[unit_name == treated_unit, full_spec_id]))
    treated_s <- length(treated_spec_ids)
    if (treated_s < min_specs) {
        stop("Treated unit has only ", treated_s, " specifications; require at least min_specs = ", min_specs, ".")
    }

    unit_spec_sets <- tau_by_spec[, .(
        spec_count = uniqueN(full_spec_id),
        spec_ids = list(sort(unique(full_spec_id)))
    ), by = .(unit_name, unit_type)]

    units_initial <- uniqueN(unit_spec_sets$unit_name)
    placebos_initial <- uniqueN(unit_spec_sets[unit_type == "control", unit_name])
    specs_initial <- treated_s

    spec_match <- data.table::copy(unit_spec_sets)
    spec_match[, grid_matches_treated :=
        (spec_count == treated_s) &
        vapply(spec_ids, identical, logical(1), treated_spec_ids)
    ]

    dropped_units <- data.table::data.table(
        unit_name = character(),
        unit_type = character(),
        reason = character()
    )
    kept_units <- character()
    final_spec_ids <- treated_spec_ids

    if (grid_policy == "strict") {
        bad_units_dt <- spec_match[unit_type == "control" & !grid_matches_treated]
        if (nrow(bad_units_dt) > 0) {
            bad_units <- paste(utils::head(bad_units_dt$unit_name, 10), collapse = ", ")
            stop("All units must share the exact same specification grid. Problematic units: ", bad_units)
        }
        kept_units <- spec_match$unit_name
    } else if (grid_policy == "drop_incomplete_units") {
        kept_units <- spec_match[unit_type == "treated" | grid_matches_treated, unit_name]
        dropped_units <- spec_match[unit_type == "control" & !grid_matches_treated, .(
            unit_name = unit_name,
            unit_type = unit_type,
            reason = "grid_mismatch_vs_treated"
        )]
    } else if (grid_policy == "intersect_specs") {
        unit_spec_lists <- unit_spec_sets$spec_ids
        final_spec_ids <- Reduce(intersect, unit_spec_lists)
        final_spec_ids <- sort(final_spec_ids)
        kept_units <- spec_match$unit_name
    } else {
        stop("Unhandled grid_policy: ", grid_policy)
    }

    s <- length(final_spec_ids)
    if (s < min_specs) {
        stop("After applying grid_policy='", grid_policy, "', only ", s,
             " shared specifications remain; require at least min_specs = ", min_specs, ".")
    }

    kept_placebos <- uniqueN(spec_match[unit_name %in% kept_units & unit_type == "control", unit_name])
    if (kept_placebos < min_placebos) {
        stop("After applying grid_policy='", grid_policy, "', only ", kept_placebos,
             " placebo units remain; require at least min_placebos = ", min_placebos, ".")
    }

    kept_units <- sort(unique(kept_units))
    tau_by_spec <- tau_by_spec[unit_name %in% kept_units & full_spec_id %in% final_spec_ids]

    pre_rmspe_by_spec <- NULL
    if ("pre_rmspe_percentile" %in% weighting) {
        if (!"rmse" %in% names(analysis_dt)) {
            stop("weighting='pre_rmspe_percentile' requires a finite 'rmse' column in filtered_results.")
        }
        pre_rmspe_by_spec <- analysis_dt[, .(
            n_unique_rmse = uniqueN(rmse),
            pre_rmspe = rmse[1]
        ), by = .(full_spec_id, unit_name, unit_type)]
        pre_rmspe_by_spec <- pre_rmspe_by_spec[
            unit_name %in% kept_units & full_spec_id %in% final_spec_ids
        ]

        inconsistent_rmse <- pre_rmspe_by_spec[n_unique_rmse != 1]
        if (nrow(inconsistent_rmse) > 0) {
            bad_msg <- paste(utils::head(
                paste0(inconsistent_rmse$unit_name, "@", inconsistent_rmse$full_spec_id), 10
            ), collapse = ", ")
            stop("rmse must be unique per unit/specification. Violations: ", bad_msg)
        }

        bad_pre <- pre_rmspe_by_spec[is.na(pre_rmspe) | is.nan(pre_rmspe) | is.infinite(pre_rmspe)]
        if (nrow(bad_pre) > 0) {
            bad_msg <- paste(utils::head(
                paste0(bad_pre$unit_name, "@", bad_pre$full_spec_id), 10
            ), collapse = ", ")
            stop("pre_rmspe contains NA/NaN/Inf for: ", bad_msg)
        }

        pre_spec_sets <- pre_rmspe_by_spec[, .(
            spec_count = uniqueN(full_spec_id),
            spec_ids = list(sort(unique(full_spec_id)))
        ), by = .(unit_name, unit_type)]
        bad_pre_specs <- pre_spec_sets[
            spec_count != s | !vapply(spec_ids, identical, logical(1), treated_spec_ids)
        ]
        if (nrow(bad_pre_specs) > 0) {
            bad_units <- paste(utils::head(bad_pre_specs$unit_name, 10), collapse = ", ")
            stop("pre_rmspe grid does not match treated specification grid for units: ", bad_units)
        }
    }

    units_in_order <- unit_spec_sets[unit_name %in% kept_units][order(unit_type, unit_name), unit_name]
    stats_list <- vector("list", length(units_in_order))

    for (i in seq_along(units_in_order)) {
        unit_i <- units_in_order[[i]]

        tau_i <- tau_by_spec[unit_name == unit_i][
            match(final_spec_ids, full_spec_id), tau_s
        ]
        if (anyNA(tau_i)) {
            stop("Missing tau_s values after specification alignment for unit ", unit_i, ".")
        }

        unit_type_i <- tau_by_spec[unit_name == unit_i, unit_type][1]
        pre_i <- NULL
        if ("pre_rmspe_percentile" %in% weighting) {
            pre_i <- pre_rmspe_by_spec[unit_name == unit_i][
                match(final_spec_ids, full_spec_id), pre_rmspe
            ]
            if (anyNA(pre_i)) {
                stop("Missing pre_rmspe values after specification alignment for unit ", unit_i, ".")
            }
        }

        rows <- list()
        for (stat_i in curve_stat) {
            for (w_i in weighting) {
                weights_i <- NULL
                if (w_i == "pre_rmspe_percentile") {
                    if (is.null(pre_i)) {
                        stop("Missing pre_rmspe values for weighting='pre_rmspe_percentile' in unit ", unit_i, ".")
                    }
                    weights_i <- compute_pre_rmspe_percentile_weights(pre_i)
                }
                estimate_i <- apply_curve_statistic(stat_i, tau_i, weights = weights_i)

                rows[[length(rows) + 1L]] <- data.table(
                    unit_name = unit_i,
                    unit_type = unit_type_i,
                    curve_statistic = stat_i,
                    weighting = w_i,
                    estimate = estimate_i,
                    n_specs = s
                )
            }
        }
        stats_list[[i]] <- data.table::rbindlist(rows)
    }
    stats_by_unit <- data.table::rbindlist(stats_list)

    stat_levels <- c("abs_median_tau", "median_abs_tau", "consistency_ratio")
    stat_combos <- unique(stats_by_unit[, .(curve_statistic, weighting)])
    stat_combos[, stat_rank := match(curve_statistic, stat_levels)]
    if (any(is.na(stat_combos$stat_rank))) {
        bad_stats <- unique(stat_combos[is.na(stat_rank), curve_statistic])
        stop("Unknown curve_statistic value(s): ", paste(bad_stats, collapse = ", "))
    }
    stat_combos[, weight_rank := ifelse(weighting == "none", 1L, 2L)]
    data.table::setorder(stat_combos, stat_rank, weight_rank)

    treated_summary_list <- vector("list", nrow(stat_combos))
    for (k in seq_len(nrow(stat_combos))) {
        stat_k <- stat_combos$curve_statistic[k]
        weight_k <- stat_combos$weighting[k]
        observed <- stats_by_unit[
            unit_name == treated_unit & curve_statistic == stat_k & weighting == weight_k, estimate
        ]
        if (length(observed) != 1) {
            stop("Expected exactly one treated-unit value for curve statistic '", stat_k,
                 "' with weighting '", weight_k, "'.")
        }
        placebo <- stats_by_unit[
            unit_type == "control" & curve_statistic == stat_k & weighting == weight_k, estimate
        ]
        p_info <- compute_placebo_curve_pvalue(
            observed_stat = observed,
            placebo_stats = placebo,
            two_sided = two_sided,
            expected_direction = expected_direction,
            stat_name = stat_k
        )
        treated_summary_list[[k]] <- data.table(
            curve_statistic = stat_k,
            weighting = weight_k,
            estimate = observed,
            p_value = p_info$p_value,
            n_extreme = p_info$n_extreme,
            n_placebos = p_info$n_placebos,
            two_sided = two_sided
        )
    }
    treated_summary <- data.table::rbindlist(treated_summary_list)
    treated_summary[, stat_rank := match(curve_statistic, stat_levels)]
    if (any(is.na(treated_summary$stat_rank))) {
        bad_stats <- unique(treated_summary[is.na(stat_rank), curve_statistic])
        stop("Unknown curve_statistic value(s) in treated summary: ", paste(bad_stats, collapse = ", "))
    }
    treated_summary[, weight_rank := ifelse(weighting == "none", 1L, 2L)]
    data.table::setorder(treated_summary, stat_rank, weight_rank)
    treated_summary[, c("stat_rank", "weight_rank") := NULL]

    attrition_report <- data.table::data.table(
        grid_policy = grid_policy,
        n_units_initial = units_initial,
        n_units_kept = uniqueN(kept_units),
        n_units_dropped = nrow(dropped_units),
        n_placebos_initial = placebos_initial,
        n_placebos_kept = kept_placebos,
        n_specs_treated = specs_initial,
        n_specs_final = s
    )
    dropped_spec_ids <- setdiff(treated_spec_ids, final_spec_ids)

    list(
        treated_summary = treated_summary,
        stats_by_unit = stats_by_unit,
        attrition_report = attrition_report,
        dropped_units = dropped_units,
        dropped_spec_ids = dropped_spec_ids
    )
}
