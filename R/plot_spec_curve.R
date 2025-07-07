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
#' @param outcomes Character vector. Names of the outcome variables to plot.
#' @param normalize_outcomes Character. Method for normalizing treatment effects for comparability.
#'   Options: "none" (default), "percent_of_mean", "percent_of_preperiod", "standardized", "rmspe_ratio".
#'   - "percent_of_mean": (tau / mean(pre-period treated)) * 100
#'   - "percent_of_preperiod": (tau / final_preperiod_value) * 100  
#'   - "standardized": tau / sd(control_effects)
#'   - "rmspe_ratio": post_RMSPE / pre_RMSPE (Abadie-style)
#' @param rmse_threshold Numeric. Threshold for root mean square error to filter results. Default is Inf.
#' @param shap_values Data.table or NULL. SHAP values from run_xgboost_shap_analysis().
#'   Should contain full_spec_id column for perfect alignment. Default is NULL.
#' @param file_path_save Character or NA. File path to save the plot. If NA, plot is not saved. Default is NA.
#' @param width Numeric. Width of the saved plot in inches. Default is 6.
#' @param height Numeric. Height of the saved plot in inches. Default is 10.
#' @param show_pvalues Logical. Whether to include significance coloring for treatment effects.
#'   Uses Abadie p-values by default, or bootstrap p-values if available. Default is TRUE.
#' @param p_threshold Numeric. P-value threshold for significance coloring. Default is 0.05.
#' @param pvalue_style Character. Style for p-value visualization. Options: "continuous" (default) or "categorical".
#'   - "continuous": Uses -log10(p-value) color scale with continuous gradation
#'   - "categorical": Uses discrete significance categories (original behavior)
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
#'   test_statistic = "rmse_ratio",
#'   pvalue_style = "continuous"
#' )
#' 
#' # 4. Plot using treatment effect test statistic p-values
#' plot_treatment_effect <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_results$shapley,
#'   test_statistic = "treatment_effect",
#'   pvalue_style = "continuous"
#' )
#' 
#' # 5. Plot using normalized treatment effect test statistic p-values
#' plot_normalized_te <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_results$shapley,
#'   test_statistic = "normalized_te",
#'   pvalue_style = "categorical"
#' )
#' }
plot_spec_curve <- function(
    long_data,
    name_treated_unit,
    outcomes,
    normalize_outcomes = "none",
    rmse_threshold = Inf,
    shap_values = NULL,
    file_path_save = NA,
    width = 6,
    height = 10,
    show_pvalues = TRUE,
    p_threshold = 0.05,
    pvalue_style = "continuous",
    prefer_bootstrap_pvalues = FALSE,
    test_statistic = "rmse_ratio",
    null_distribution = "placebo",
    crop_outliers = "none",
    sort_by = "tau"
) {

  # Libraries imported via NAMESPACE 

  # Input validation
  valid_test_statistics <- c("rmse_ratio", "treatment_effect", "normalized_te")
  if (!test_statistic %in% valid_test_statistics) {
    stop("test_statistic must be one of: ", paste(valid_test_statistics, collapse = ", "))
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
    
    long_data <- main_data
  } else if (!data.table::is.data.table(long_data)) {
    long_data <- data.table::as.data.table(long_data)
  }
  
  # Check required columns
  required_cols <- c("unit_name", "tau", "post_period", "full_spec_id")
  missing_cols <- setdiff(required_cols, names(long_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in long_data:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check if p-value data is available
  has_abadie_pvalues <- "p_value" %in% names(long_data)
  has_bootstrap_pvalues <- "p_value_two_tailed" %in% names(long_data)
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
  
  # Filter by outcomes
  if ("outcome" %in% names(long_data)) {
    sc_results_df <- long_data[outcome %in% outcomes]
  } else {
    sc_results_df <- long_data
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
  panel_b_data[feature == 'simplex' & feature_group=='Weight\nMethod', feature := "Original"]
  panel_b_data[feature == 'lasso' & feature_group=='Weight\nMethod', feature := "Original + Penalty Lasso"]
  panel_b_data[feature == 'ridge' & feature_group=='Weight
Method', feature := "Original + Penalty Ridge"]
  panel_b_data[feature == 'ols' & feature_group=='Weight
Method', feature := "OLS Weights"]

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
          midpoint = -log10(p_threshold),  # Threshold at p=0.05
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
        # p_value < 0.01, "p < 0.01)",
        # p_value < 0.05, "Significant (p < 0.05)", 
        # p_value < 0.10, "Marginally Significant (p < 0.10)",
        # p_value >= 0.10, "Not Significant (p >= 0.10)",
        # default = "Unknown"
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
          # name = "Statistical Significance",
          # values = c("p < 0.01)" = "#08519c", "Significant (p < 0.05)" = "#3182bd",
          #            "Marginally Significant (p < 0.10)" = "#fd8d3c", "Not Significant (p >= 0.10)" = "#d94701",
          #            "Unknown" = "#999999"),
          # breaks = c("p < 0.01)", "Significant (p < 0.05)",
          #            "Marginally Significant (p < 0.10)", "Not Significant (p >= 0.10)")
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

# Apply y-axis limits if outlier cropping is requested
if (!is.null(y_limits)) {
  p1 <- p1 + ylim(y_limits[1], y_limits[2])
}

# Now you can print p1 and combine it with p2 using plot_grid as before.
# print(p1)
# p2: Specification Choices and Shapley Values (Bottom) 
# Panel B: Standard SHAP visualization (no significance strip)
p2 <- ggplot(plot_data_p2, aes(x = Specification, y = feature)) +
  geom_point(aes(color = shapley_value), shape = 15, size = 2.5) +
  scale_color_gradient2(
    name = "Shapley Value",
    low = "#762a83",    # Purple (better contrast)
    mid = "#f7f7f7",    # Light gray  
    high = "#1b7837",   # Green (better contrast)
    midpoint = 0
  )

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
    legend.box = "horizontal"  # Place legends side by side
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

return(final_plot)
}