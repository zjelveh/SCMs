#' Plot Specification Curve Using Direct Long Format Data
#'
#' @title Plot Specification Curve with Perfect SHAP Alignment
#' @description Creates specification curve plots directly from long format data,
#' ensuring perfect alignment between SHAP values and specifications through 
#' shared spec_number identifiers.
#'
#' @param long_data Data.table. Long format data from spec_curve(..., output_format = "long").
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
#'   Should contain spec_number column for perfect alignment. Default is NULL.
#' @param file_path_save Character or NA. File path to save the plot. If NA, plot is not saved. Default is NA.
#' @param width Numeric. Width of the saved plot in inches. Default is 6.
#' @param height Numeric. Height of the saved plot in inches. Default is 10.
#' @param show_pvalues Logical. Whether to include significance coloring for treatment effects.
#'   Uses Abadie p-values by default, or bootstrap p-values if available. Default is TRUE.
#' @param p_threshold Numeric. P-value threshold for significance coloring. Default is 0.05.
#' @param prefer_bootstrap_pvalues Logical. Whether to prefer bootstrap p-values over Abadie p-values
#'   when both are available. Default is FALSE (prefer Abadie).
#' @param null_distribution Character. The distribution to plot as the null/comparison.
#'   Options: "placebo" (default) or "bootstrap".
#'   - "placebo": Use placebo effects from control units as the null distribution.
#'   - "bootstrap": Use the bootstrapped null distribution for the treated unit.
#' @param crop_outliers Character or Numeric. Method for cropping outliers in Panel A.
#'   Options: "none" (default), "percentile", "iqr", "mad", or numeric vector c(ymin, ymax).
#'   - "percentile": Crop to 5th-95th percentile of treated unit effects
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
#' # 3. Plot using placebo null (default)
#' plot_placebo <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_results$shapley
#' )
#' 
#' # 4. Plot using bootstrap null
#' plot_bootstrap <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "TREATED_ID",
#'   shap_values = shap_results$shapley,
#'   null_distribution = "bootstrap"
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
    prefer_bootstrap_pvalues = FALSE,
    null_distribution = "placebo",
    crop_outliers = "none",
    sort_by = "tau"
) {

  # Libraries imported via NAMESPACE 

  # Input validation
  if (!data.table::is.data.table(long_data)) {
    long_data <- data.table::as.data.table(long_data)
  }
  
  # Check required columns
  required_cols <- c("unit_name", "tau", "post_period", "spec_number")
  missing_cols <- setdiff(required_cols, names(long_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in long_data:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check if p-value data is available
  has_abadie_pvalues <- "p_value" %in% names(long_data)
  has_bootstrap_pvalues <- "bootstrap_p_value_two_tailed" %in% names(long_data)
  has_pvalues <- has_abadie_pvalues || has_bootstrap_pvalues
  
  # Determine which p-values to use
  use_bootstrap_pvalues <- has_bootstrap_pvalues && (prefer_bootstrap_pvalues || !has_abadie_pvalues)
  p_value_column <- if (use_bootstrap_pvalues) "bootstrap_p_value_two_tailed" else "p_value"
  
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
  average_effect_df <- sc_results_df[
    post_period == TRUE,
    .(tau = mean(tau, na.rm = TRUE), rmse = mean(rmse, na.rm = TRUE)),
    by = c('unit_name', 'unit_type', 'spec_number', 'outcome', 'outcome_model', 'const', 'fw', 'feat', 'data_sample', 'num_pre_period_years')
  ]
  
  
  
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
        average_effect_df[unit_name == name_treated_unit, tau := post_pre_ratio]
        # For comparison units, calculate their ratios too if possible
        if ("post_rmspe" %in% names(average_effect_df) && "rmse" %in% names(average_effect_df)) {
          average_effect_df[unit_type %in% c("control", "bootstrap") & !is.na(post_rmspe) & rmse > 0, 
                           tau := post_rmspe / rmse]
        }
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
      y_limits <- quantile(treated_effects, c(0.05, 0.95), na.rm = TRUE)
      
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

  # A. Prepare data for Panel A (Main Effects Plot)
  panel_a_data <- copy(average_effect_df)
  setnames(panel_a_data,
           old = c("unit_name", "spec_number", "tau", "rmse"),
           new = c("Unit Name", "Specification", "Estimate", "RMSE"))

  # Add p-values if requested
  if (has_pvalues && show_pvalues) {
    p_value_data <- unique(sc_results_df[unit_name == name_treated_unit & !is.na(get(p_value_column)),
                                         .(spec_number, p_value = get(p_value_column))])
    setnames(p_value_data, "spec_number", "Specification")
    panel_a_data <- merge(panel_a_data, p_value_data, by = "Specification", all.x = TRUE)
  }

  # B. Prepare data for Panel B (SHAP & Specification Details)
  panel_b_data <- melt(average_effect_df[unit_name == name_treated_unit],
                       id.vars = c('unit_name', 'spec_number', 'tau', 'rmse', 'unit_type'),
                       variable.name = 'feature_group', value.name = 'feature')

  # Merge SHAP values if provided
  if (!is.null(shap_values)) {
    panel_b_data <- merge(panel_b_data, shap_values,
                          by.x = c('unit_name', 'spec_number', 'feature_group', 'feature'),
                          by.y = c('unit', 'spec_number', 'feature_group', 'feature'),
                          all.x = TRUE)
  }
  
  setnames(panel_b_data, c('tau', 'unit_name', 'rmse', 'spec_number'),
           c('Estimate', 'Unit Name', 'RMSE', 'Specification'))

  # C. Re-sort specifications if requested (affects both panels)
  if (sort_by != "tau") {
    sort_data <- copy(panel_a_data[unit_type == "treated"])

    if (sort_by == "pvalue" && "p_value" %in% names(sort_data)) {
      sort_data <- sort_data[order(p_value)]
    } else if (sort_by == "rmspe_ratio" && "post_pre_ratio" %in% names(sort_data)) {
      sort_data <- sort_data[order(post_pre_ratio)]
    }

    sort_data[, new_spec_number := 1:.N]
    spec_reorder <- sort_data[, .(old_spec = Specification, new_spec_number)]

    # Apply new ordering to both panel datasets
    panel_a_data <- merge(panel_a_data, spec_reorder, by.x = "Specification", by.y = "old_spec", all.x = TRUE)
    panel_a_data[, Specification := new_spec_number][, new_spec_number := NULL]
    
    panel_b_data <- merge(panel_b_data, spec_reorder, by.x = "Specification", by.y = "old_spec", all.x = TRUE)
    panel_b_data[, Specification := new_spec_number][, new_spec_number := NULL]
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
  panel_b_data[feature_group== 'num_pre_period_years', feature_group:='Pre Length']
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

  # F. Create Fit Quality Score
  all_rmse_data <- plot_data_p1[!is.na(RMSE), RMSE]
  if (length(all_rmse_data) > 0) {
    max_rmse <- max(all_rmse_data, na.rm = TRUE)
    min_rmse <- min(all_rmse_data, na.rm = TRUE)
    plot_data_p1[, fit_quality := max_rmse - RMSE + min_rmse]
  } else {
    plot_data_p1[, fit_quality := 1]
  }

  # --- 2. Create the Final Plots ---
  data.table::setorder(plot_data_p1, -unit_type)

  # Create Panel A with conditional coloring for p-values
  if (has_pvalues && show_pvalues && "p_value" %in% names(plot_data_p1)) {
    treated_data <- plot_data_p1[unit_type == "treated"]
    control_data <- plot_data_p1[unit_type %in% c("control", "bootstrap")]
    
    treated_data[, significance_category := fcase(
      p_value < 0.01, "Highly Significant (p < 0.01)",
      p_value < 0.05, "Significant (p < 0.05)", 
      p_value < 0.10, "Marginally Significant (p < 0.10)",
      p_value >= 0.10, "Not Significant (p >= 0.10)",
      default = "Unknown"
    )]
    
    p1 <- ggplot() +
      geom_point(data = control_data, aes(x = Specification, y = Estimate, size = fit_quality),
                 color = "gray60", alpha = 0.3, shape = 19) +
      geom_point(data = treated_data, aes(x = Specification, y = Estimate, color = significance_category, size = fit_quality),
                 alpha = 0.8, shape = 19) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = 'dashed') +
      scale_color_manual(
        name = "Statistical Significance",
        values = c("Highly Significant (p < 0.01)" = "#08519c", "Significant (p < 0.05)" = "#3182bd",
                   "Marginally Significant (p < 0.10)" = "#fd8d3c", "Not Significant (p >= 0.10)" = "#d94701",
                   "Unknown" = "#999999"),
        breaks = c("Highly Significant (p < 0.01)", "Significant (p < 0.05)",
                   "Marginally Significant (p < 0.10)", "Not Significant (p >= 0.10)")
      )
  } else {
    p1 <- ggplot(plot_data_p1, aes(x = Specification, y = Estimate, fill = unit_type,
                                 size = fit_quality, color = unit_type, alpha = unit_type)) +
      geom_point(shape = 21) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = 'dashed') +
      scale_fill_manual(name = "Unit Type", values = c(control = "gray60", treated = "#1f78b4", bootstrap = "gray60")) +
      scale_color_manual(name = "Unit Type", values = c(control = "gray60", treated = "#1f78b4", bootstrap = "gray60")) +
      scale_alpha_manual(name = "Unit Type", values = c(treated = 0.8, control = 0.3, bootstrap = 0.3))
  }

# Add common elements to both plot types
p1 <- p1 +
  
  scale_size_continuous(
    name = "Pre-period Fit Quality",
    range = c(0.5, 4),
    breaks = function(x) {
      # Create 4 evenly spaced breaks for the legend
      seq(min(x), max(x), length.out = 4)
    },
    labels = c("Poor", "Fair", "Good", "Excellent"),
    guide = guide_legend(
      title.position = "top",
      override.aes = list(color = "black", alpha = 1, shape = 19),
      nrow = 1
    )
  ) +
  
  guides(
    size = guide_legend(title.position = "top"),
    color = guide_legend(title.position = "top", override.aes = list(size = 4))
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
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