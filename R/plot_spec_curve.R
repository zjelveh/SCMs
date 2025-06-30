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
#' @param normalize_outcomes Logical. Whether to normalize outcomes. Default is FALSE.
#' @param rmse_threshold Numeric. Threshold for root mean square error to filter results. Default is Inf.
#' @param shap_values Data.table or NULL. SHAP values from run_xgboost_shap_analysis().
#'   Should contain spec_number column for perfect alignment. Default is NULL.
#' @param file_path_save Character or NA. File path to save the plot. If NA, plot is not saved. Default is NA.
#' @param width Numeric. Width of the saved plot in inches. Default is 6.
#' @param height Numeric. Height of the saved plot in inches. Default is 10.
#'
#' @return A ggplot object representing the specification curve.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # New recommended workflow with perfect alignment
#' 
#' # 1. Generate long format data
#' spec_results_long <- spec_curve(
#'   dataset = your_data,
#'   outcomes = "num_homicide",
#'   name_treated_unit = "PAPEP0000",
#'   output_format = "long"
#' )
#' 
#' # 2. Run SHAP analysis with perfect alignment
#' config <- create_catboost_config(
#'   dataset_name = "homicide_study",
#'   treated_unit_name = "PAPEP0000"
#' )
#' shap_results <- run_catboost_shap_analysis(spec_results_long, config)
#' 
#' # 3. Plot with guaranteed perfect alignment
#' plot <- plot_spec_curve(
#'   long_data = spec_results_long,
#'   name_treated_unit = "PAPEP0000",
#'   outcomes = "num_homicide",
#'   shap_values = shap_results$shapley
#' )
#' }
plot_spec_curve <- function(
    long_data,
    name_treated_unit,
    outcomes,
    normalize_outcomes = FALSE,
    rmse_threshold = Inf,
    shap_values = NULL,
    file_path_save = NA,
    width = 6,
    height = 10
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
    by = c('unit_name', 'spec_number', 'outcome', 'outcome_model', 'const', 'fw', 'feat', 'data_sample', 'num_pre_period_years')
  ]
  
  # Add unit type information
  average_effect_df[, unit_type := ifelse(unit_name == name_treated_unit, "treated", "control")]
  
  # Normalize outcomes if requested
  if (normalize_outcomes) {
    # Calculate standard deviation from control units
    sd_outcome <- sd(average_effect_df[unit_type == "control"]$tau, na.rm = TRUE)
    if (sd_outcome > 0) {
      average_effect_df[, tau := tau / sd_outcome]
    }
  }

average_effect_df_long = melt(average_effect_df, c('unit_name', 'spec_number', 'tau', 'rmse', 'unit_type'), 
                              variable.name = 'feature_group', value.name = 'feature')
average_effect_df_long = average_effect_df_long[unit_name==name_treated_unit & feature_group%in%shap_values$feature_group]

average_effect_df_long_shap = merge(average_effect_df_long, 
                                    shap_values, by.x=c('unit_name', 'spec_number',
                                                        'feature_group', 'feature'),
                                    by.y=c('unit', 'spec_number', 'feature_group', 'feature'))


setnames(average_effect_df_long_shap, c('tau', 'unit_name', 'rmse', 'spec_number'),
         c('Estimate', 'Unit Name', 'RMSE', 'Specification'))
# Rename columns to match original plotting expectations
average_effect_df_long_shap[feature_group=='const', feature_group:= 'Weight\nMethod']
average_effect_df_long_shap[feature_group=='outcome', feature_group:=  'Outcome']
average_effect_df_long_shap[feature_group=='outcome_model', feature_group:= 'Outcome\nModel']
average_effect_df_long_shap[feature_group== 'fw', feature_group:= 'V Weights']
average_effect_df_long_shap[feature_group== 'feat', feature_group:= 'Features']
average_effect_df_long_shap[feature_group== 'data_sample', feature_group:='Donor Pool']
average_effect_df_long_shap[feature_group== 'num_pre_period_years', feature_group:='Pre Length']

# Transform constraint names to match original format
average_effect_df_long_shap[feature == 'simplex' & feature_group=='Weight\nMethod', feature := "Original"]
average_effect_df_long_shap[feature == 'lasso' & feature_group=='Weight\nMethod', feature := "Original + Penalty Lasso"]
average_effect_df_long_shap[feature == 'ridge' & feature_group=='Weight\nMethod', feature := "Original + Penalty Ridge"]
average_effect_df_long_shap[feature == 'ols' & feature_group=='Weight\nMethod', feature := "OLS Weights"]
average_effect_df_long_shap[, Type := 'Average']

control_effects = average_effect_df[unit_type=='control']


# --- 1. Data Preparation ---

# Assume 'average_effect_df_long_shap' and 'control_effects' are your starting data.tables

# Prepare data for the top plot (p1)
# Get unique estimates for the treated unit
treated_estimates <- unique(average_effect_df_long_shap[, .(
  Specification,
  Estimate,
  RMSE,
  unit_type
)])

# Standardize column names for control unit data
setnames(control_effects,
         old = c("unit_name", "spec_number", "tau", "rmse"),
         new = c("Unit Name", "Specification", "Estimate", "RMSE"))


# Combine treated and control data for the main effect plot
plot_data_p1 <- rbindlist(list(
  treated_estimates,
  control_effects[, .(Specification, Estimate, RMSE, unit_type)]
), use.names = TRUE, fill = TRUE)

# Ensure 'unit_type' is a factor to control plotting order and legend
plot_data_p1[, unit_type := factor(unit_type, levels = c("treated", "control"))]

# Prepare data for the bottom plot (p2)
# The data is already in the correct "long" format for the tick marks
plot_data_p2 <- average_effect_df_long_shap

# FIX: Calculate RMSE quartiles across ALL units, not just treated.
# This uses all non-missing RMSE values to define the scale.
all_rmse_data <- plot_data_p1 %>%
  filter(!is.na(RMSE)) %>%
  pull(RMSE)

# Ensure there are no issues with non-finite values for quantile calculation
if (length(all_rmse_data) > 0) {
  rmse_quartiles <- quantile(all_rmse_data, probs = c(0, 1/3, 2/3, 1.0), na.rm = TRUE)
} else {
  # Provide a fallback if there's no RMSE data, to prevent crashes
  rmse_quartiles <- c(0, 1)
}

plot_data_p1[, bin_q:= cut(RMSE, rmse_quartiles)]

plot_data_p1[, circle_size:=mean(RMSE) , by=.(bin_q)]

# --- 2. Create the Final Plots ---
data.table::setorder(plot_data_p1, -unit_type)

p1 <- plot_data_p1 %>%
  # Set up all aesthetics in the main ggplot() call
  ggplot(aes(
    x = Specification,
    y = Estimate,
    fill = unit_type,
    size = RMSE,
    color = unit_type,
    alpha = unit_type
  )) +
  
  # Draw all points using a single layer
  geom_point(shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.5, linetype = 'dashed') +
  
  # --- Define the Scales for Each Mapped Aesthetic ---
  
  scale_fill_manual( 
    name = "Unit Type",
    values = c(control = "firebrick3", treated = "deepskyblue3")
  ) +

  scale_color_manual( 
    name = "Unit Type",
    values = c(control = "firebrick3", treated = "deepskyblue3")
  ) +
  
  scale_size(
    name = "Pre-period Fit (RMSE)",
    range = c(0.5, 5), # Maps the data domain to a visual size range of 5 (large) to 0.5 (small)
    breaks = sort(unique(plot_data_p1[!is.na(bin_q)]$circle_size)),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  
  scale_alpha_manual(
    name = "Unit Type", # Use the same name as `scale_fill_manual` to merge legends
    values = c(treated = 0.9, control = 0.2)
  ) +
  
  # --- Final Touches ---
  
  # A simpler guides call to ensure the legends look correct
  guides(
    alpha = "none", # Don't draw a separate legend for alpha
    size = guide_legend(title.position = "top"),
    fill = guide_legend(override.aes = list(alpha = 1, size = 4)) # Make legend glyphs opaque and larger
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(colour = "black")
  ) +
  labs(x = NULL, y = 'Average Treatment Effect') +
  ylim(-200, 200)# + geom_hline(yintercept= -1.6)
  # xlim(-10, 10)

# Now you can print p1 and combine it with p2 using plot_grid as before.
# print(p1)
# p2: Specification Choices and Shapley Values (Bottom) - No changes
# ... (p2 code remains the same as the previous version) ...
p2 <- ggplot(plot_data_p2, aes(x = Specification, y = feature)) +
  geom_point(aes(color = shapley_value), shape = 15, size = 2.5) +
  scale_color_gradient2(
    name = "Shapley Value",
    low = "#7b3294",
    mid = "grey85",
    high = "#008837",
    midpoint = 0
  ) +
  facet_grid(feature_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_minimal() +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.text.y = element_text(colour = "black", size = 8),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    panel.spacing = unit(.25, "lines"),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
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