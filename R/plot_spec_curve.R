#' Plot Specification Curve Using Direct Long Format Data
#'
#' @title Plot Specification Curve with Perfect SHAP Alignment
#' @description Creates specification curve plots directly from long format data,
#' ensuring perfect alignment between SHAP values and specifications through
#' shared full_spec_id identifiers.
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
#' @param shap_values Data.table or NULL. SHAP values from run_catboost_shap_analysis().
#'   Should contain full_spec_id column for perfect alignment. Default is NULL.
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
#'   - "mad": Crop to median ± 3*MAD (robust outlier detection)
#'   - c(ymin, ymax): Manual y-axis limits
#' @param sort_by Character. Method for sorting specifications on x-axis.
#'   Options: "tau" (default), "pvalue", "rmspe_ratio".
#'   - "tau": Sort by treatment effect magnitude (original behavior)
#'   - "pvalue": Sort by statistical significance (most significant first)
#'   - "rmspe_ratio": Sort by post/pre RMSPE ratio
#' @param calculate_shap_pvalues Logical. Whether to automatically calculate SHAP significance when
#'   multi-unit SHAP data is available. When TRUE (default), if SHAP values for control units are
#'   found, calculates p-values for treated unit's feature importance and displays significance via
#'   point transparency. Set to FALSE to disable SHAP significance testing.
#' @param shap_pvalue_type Character. Method for calculating SHAP p-values.
#'   Options: "absolute" (default) or "signed".
#'   - "absolute": Calculates p-values based on the absolute magnitude of SHAP values.
#'   - "signed": Calculates p-values based on the signed magnitude, considering positive and negative effects separately.
#'
#' @return A ggplot object representing the specification curve.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # New recommended workflow with perfect alignment
#'
#' # 1. Generate long format data with bootstrap inference
#' spec_results_long <- run_spec_curve_analysis(
#'   dataset,
#'   params = params,
#'   inference_type = "bootstrap"
#' )
#'
#' # 2. Run SHAP analysis
#' shap_results <- run_catboost_shap_analysis(spec_results_long, config)
#'
#' # 3. Plot using placebo null with continuous p-values (default - RMSE ratio)
#' plot_rmse_ratio <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_results$shapley,
#'   test_statistic = "rmse_ratio"
#' )
#'
#' # 4. Plot using treatment effect test statistic p-values
#' plot_treatment_effect <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_results$shapley,
#'   test_statistic = "treatment_effect"
#' )
#'
#' # 5. Plot using normalized treatment effect test statistic p-values
#' plot_normalized_te <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_results$shapley,
#'   test_statistic = "normalized_te"
#' )
#'
#' # 6. Plot with SHAP significance testing (requires multi-unit SHAP data)
#' # Run CatBoost analysis on all units
#' config_all_units <- create_catboost_config(
#'   dataset_name = "example",
#'   treated_unit_name = "TREATED_ID",
#'   treated_unit_only = FALSE   # Analyze all units for significance testing
#' )
#' shap_all_units <- run_catboost_shap_analysis(spec_results_long, config_all_units)
#'
#' # Plot with automatic SHAP significance testing (absolute method, default)
#' plot_with_shap_significance_abs <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_all_units$shapley,   # Multi-unit SHAP data
#'   calculate_shap_pvalues = TRUE   # Enable significance testing (default)
#' )
#' # Displays SHAP values with transparency indicating statistical significance (absolute)
#'
#' # Plot with SHAP significance testing using signed p-values
#' plot_with_shap_significance_signed <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_all_units$shapley,
#'   calculate_shap_pvalues = TRUE,
#'   shap_pvalue_type = "signed" # New option
#' )
#'
#' # 7. Disable SHAP significance testing if not wanted
#' plot_no_shap_significance <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_all_units$shapley,
#'   calculate_shap_pvalues = FALSE   # Disable significance testing
#' )
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
    show_pvalues = TRUE,
    p_threshold = 0.05,
    prefer_bootstrap_pvalues = FALSE,
    test_statistic = "rmse_ratio",
    null_distribution = "placebo",
    crop_outliers = "none",
    sort_by = "tau",
    calculate_shap_pvalues = TRUE,
    shap_pvalue_type = "absolute" # NEW PARAMETER
) {

    # Libraries imported via NAMESPACE

    # Input validation
    valid_test_statistics <- c("rmse_ratio", "treatment_effect", "normalized_te")
    if (!test_statistic %in% valid_test_statistics) {
        stop("test_statistic must be one of: ", paste(valid_test_statistics, collapse = ", "))
    }
    valid_shap_pvalue_types <- c("absolute", "signed") # NEW VALIDATION
    if (!shap_pvalue_type %in% valid_shap_pvalue_types) {
        stop("shap_pvalue_type must be one of: ", paste(valid_shap_pvalue_types, collapse = ", "))
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

    # Filter by outcomes (if specified)
    if (!is.null(outcomes) && "outcome" %in% names(main_data)) {
        sc_results_df <- main_data[outcome %in% outcomes]
    } else {
        sc_results_df <- main_data
    }

    # Filter by RMSE if threshold provided
    if ("rmse" %in% names(sc_results_df) && rmse_threshold < Inf) {
        sc_results_df <- sc_results_df[rmse < rmse_threshold]
    }

    if (nrow(sc_results_df) == 0) {
        stop("No data remaining after filtering")
    }

    # Calculate average effects for plotting (post-treatment period only)
    # Check if post_pre_ratio exists before trying to aggregate it
    if ("post_pre_ratio" %in% names(sc_results_df)) {
        average_effect_df <- sc_results_df[
            post_period == TRUE,
            .(tau = mean(tau, na.rm = TRUE),
              rmse = mean(rmse, na.rm = TRUE),
              post_pre_ratio = mean(post_pre_ratio, na.rm = TRUE)),
            by = c('unit_name', 'unit_type', 'full_spec_id', 'outcome',
            'outcome_model', 'const', 'fw', 'feat', 'data_sample')
        ]
    } else {
        average_effect_df <- sc_results_df[
            post_period == TRUE,
            .(tau = mean(tau, na.rm = TRUE),
              rmse = mean(rmse, na.rm = TRUE)),
            by = c('unit_name', 'unit_type', 'full_spec_id', 'outcome',
            'outcome_model', 'const', 'fw', 'feat', 'data_sample')
        ]
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
            # Median ± 3 * Median Absolute Deviation
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
    panel_b_data <- melt(average_effect_df[unit_name == name_treated_unit],
                         id.vars = c('unit_name', 'full_spec_id', 'tau', 'rmse', 'unit_type'),
                         variable.name = 'feature_group', value.name = 'feature')

    # Apply the same specification numbering
    panel_b_data <- merge(panel_b_data, spec_mapping, by = "full_spec_id", all.x = TRUE)

    # Merge SHAP values if provided
    if (!is.null(shap_values)) {
        # Validate SHAP data structure
        required_shap_cols <- c('unit', 'full_spec_id', 'feature_group', 'feature', 'shapley_value')
        missing_shap_cols <- setdiff(required_shap_cols, names(shap_values))
        if (length(missing_shap_cols) > 0) {
            warning(paste("Missing SHAP columns:", paste(missing_shap_cols, collapse = ", "),
                          "- SHAP values will not be displayed"))
        } else {
            # Check for multi-unit SHAP data and calculate significance if requested
            unique_units <- unique(shap_values$unit)
            has_control_shaps <- length(unique_units) > 1 && any(unique_units != name_treated_unit)

            if (has_control_shaps && calculate_shap_pvalues) {
                message("Found SHAP values for ", length(unique_units), " units (",
                        sum(unique_units != name_treated_unit), " controls). Calculating SHAP significance using '", shap_pvalue_type, "' method...")
                message("To disable SHAP significance testing, set calculate_shap_pvalues = FALSE")

                # Calculate SHAP significance on filtered data - PASS NEW PARAM
                shap_significance <- calculate_shap_significance(shap_values, name_treated_unit, pvalue_method = shap_pvalue_type)
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
            } else if (has_control_shaps && !calculate_shap_pvalues) {
                message("Multi-unit SHAP data found but significance testing disabled by calculate_shap_pvalues = FALSE")
            }

            panel_b_data <- merge(panel_b_data, shap_values,
                                  by.x = c('unit_name', 'full_spec_id', 'feature_group', 'feature'),
                                  by.y = c('unit', 'full_spec_id', 'feature_group', 'feature'),
                                  all.x = TRUE)
        }
    }

    setnames(panel_b_data, c('tau', 'unit_name', 'rmse'),
             c('Estimate', 'Unit Name', 'RMSE'))

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
    
    # Transform constraint names for display (preserve descriptive names in Features group)
    panel_b_data[feature == 'simplex' & feature_group=='Weight\nMethod', feature := "Original"]
    panel_b_data[feature == 'lasso' & feature_group=='Weight\nMethod', feature := "Original + Penalty Lasso"]
    panel_b_data[feature == 'ridge' & feature_group=='Weight\nMethod', feature := "Original + Penalty Ridge"]
    panel_b_data[feature == 'ols' & feature_group=='Weight\nMethod', feature := "OLS Weights"]
    
    # DO NOT transform feature names in the 'Features' group - preserve descriptive names as-is

    # Filter out feature_groups with only one unique feature value
    feature_counts <- panel_b_data[, .(n_unique = uniqueN(feature)), by = feature_group]
    groups_to_keep <- feature_counts[n_unique > 1, feature_group]
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
                        barwidth = unit(4, "cm"),
                        barheight = unit(0.5, "cm")
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

    # Add common elements to both plot types
    p1 <- p1 +

        guides(
            color = guide_legend(title.position = "top", override.aes = list(size = 4))
        ) +

        theme_minimal() +
        theme(
            legend.position = "left",
            legend.box = "vertical",
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.text = element_text(colour = "black")
        ) +
        labs(x = NULL, y = y_label)

    # Calculate specification curve p-values on filtered data
    # This ensures p-values reflect the actual data being plotted
    spec_curve_pvals <- NULL
    if (is.list(long_data)) {
        # Get expected direction from long_data or default to negative
        expected_direction <- if (!is.null(long_data$expected_direction)) long_data$expected_direction else "negative"

        # Get abadie inference if available
        abadie_inference <- if (!is.null(long_data$abadie_inference)) long_data$abadie_inference else NULL

        # Calculate p-values on the filtered data (sc_results_df after all filtering)
        spec_curve_pvals <- tryCatch({
            calculate_spec_curve_pvalues_filtered(
                filtered_results = sc_results_df,
                abadie_inference = abadie_inference,
                expected_direction = expected_direction,
                name_treated_unit = name_treated_unit
            )
        }, error = function(e) {
            warning("Could not calculate specification curve p-values: ", e$message)
            NULL
        })
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

            # Add annotation to lower-left of Panel A (treatment effects plot)
            p1 <- p1 +
                annotate("text",
                        x = -Inf, y = -Inf,
                        label = spec_annotation_text,
                        hjust = -0.1, vjust = -0.5,   # Left-align and move up slightly
                        size = 3, color = "#000000",
                        fontface = "plain")
        }
    }

    # Apply y-axis limits if outlier cropping is requested
    if (!is.null(y_limits)) {
        p1 <- p1 + ylim(y_limits[1], y_limits[2])
    }

    # Now you can print p1 and combine it with p2 using plot_grid as before.
    # p2: Specification Choices and Shapley Values (Bottom)
    # Panel B: SHAP visualization with optional significance transparency
    has_shap_significance <- "shap_pvalue" %in% names(plot_data_p2) && "shap_significance_level" %in% names(plot_data_p2)

    if (has_shap_significance) {
        # Create p-value-based alpha mapping for continuous significance scale
        # Lower p-value = more significant = more solid (higher alpha)
        # Map p-values: p=0 -> alpha=1.0 (solid), p=1 -> alpha=0.2 (transparent)
        plot_data_p2[, shap_alpha_value := pmax(0.2, 1.0 - pmin(shap_pvalue, 1.0))]
    

        p2 <- ggplot(plot_data_p2, aes(x = Specification, y = feature)) +
            geom_point(aes(color = shapley_value, alpha = shap_alpha_value), shape = 15, size = 2.5) +
            scale_color_gradient2(
                name = "Shapley Value",
                low = "#CA0020",    # Red (better contrast)
                mid = "#969696",    # Medium gray instead of light gray
                high = "#0571B0",   # Blue (better contrast)
                midpoint = 0
            ) +
            scale_alpha_continuous(
                name = "SHAP Significance\n(Solid = Low p-value)",
                range = c(0.2, 1.0),
                breaks = c(0.2, 0.5, 0.95, 0.99, 0.999),
                labels = c("p~1.0", "p~0.5", "p~0.05", "p~0.01", "p~0.0")
            )
    } else {
        p2 <- ggplot(plot_data_p2, aes(x = Specification, y = feature)) +
            geom_point(aes(color = shapley_value), shape = 15, size = 2.5) +
            scale_color_gradient2(
                name = "Shapley Value",
                low = "#CA0020",    # Red (better contrast)
                mid = "#969696",    # Medium gray instead of light gray
                high = "#0571B0",   # Blue (better contrast)
                midpoint = 0
            )
    }

    # Add remaining plot elements
    p2 <- p2 +
        facet_grid(feature_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
        theme_minimal() +
        theme(
            axis.line.x = element_line(color = "black", linewidth = 0.5),
            axis.text.y = element_text(colour = "black", size = 8),
            strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
            panel.spacing = unit(.25, "lines"),
            legend.position = "bottom",
            legend.key.width = unit(1.5, "cm"),
            legend.box = "horizontal"   # Place legends side by side
        ) +
        labs(x = "Specification Number", y = "")


    # --- 3. Combine Plots (No changes) ---
    final_plot <- plot_grid(
        p1,
        p2,
        ncol = 1,
        align = 'v',
        axis = "lr",
        rel_heights = c(0.4, 0.6)
    )

    # Save plot if file path is provided
    if (!is.na(file_path_save)) {
        ggplot2::ggsave(file_path_save, plot = final_plot, width = width, height = height)
    }

    return(list(final_plot=final_plot, plot_data_p1=plot_data_p1, plot_data_p2=plot_data_p2))
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
#' For each specification × feature_group combination:
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

    # Calculate significance for each specification × feature_group combination
    significance_results <- shap_data[, {

        treated_shap <- shapley_value[unit == treated_unit_name]

        if (length(treated_shap) == 0) {
            # Treated unit not found for this spec × feature_group
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