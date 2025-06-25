#' Complete Specification Curve Analysis Pipeline
#' 
#' This example demonstrates the complete workflow for running specification
#' curve analysis followed by XGBoost + SHAP analysis to understand what
#' drives variation in treatment effect estimates.

# Load required packages
library(SCMs)
library(data.table)
library(ggplot2)
library(xgboost)
library(fastshap)

# =============================================================================
# STEP 1: SPECIFICATION CURVE ANALYSIS
# =============================================================================

# Load and process data (replace with your data loading function)
# dataset <- process_homicide_data("path/to/your/data.csv")

# Create analysis configuration
analysis_config <- create_analysis_config("homicide", custom_params = list(
  name_treated_unit = "PAPEP0000",
  treated_period = 2015,
  min_period = 2010,
  end_period = 2019,
  cores = 4  # Adjust based on your system
))

# Run specification curve analysis
cat("Running specification curve analysis...\n")
spec_results <- run_spec_curve_analysis(
  dataset = dataset,
  params = analysis_config,
  output_dir = "results/"
)

# Extract and save results
compiled_results <- extract_spec_curve_results(spec_results)
save_spec_curve_results(compiled_results, "results/spec_curve_results")

# Create basic specification curve plot
plot_basic <- plot_spec_curve(spec_results, 
                             outcomes = "num_homicide", 
                             name_treated_unit = "PAPEP0000")

# =============================================================================
# STEP 2: XGBOOST + SHAP ANALYSIS 
# =============================================================================

# Create XGBoost analysis configuration
xgb_config <- create_xgboost_config(
  dataset_name = "homicide_analysis",
  file_path = "results/spec_curve_results.csv",
  treated_unit_name = "PAPEP0000",
  outcome_filter = "num_homicide",
  spec_train_test_split = TRUE,
  spec_split_ratio = 0.7
)

# Run XGBoost + SHAP analysis
cat("Running XGBoost + SHAP analysis...\n")
xgb_results <- run_xgboost_shap_analysis(xgb_config)

# =============================================================================
# STEP 3: ADVANCED VISUALIZATIONS
# =============================================================================

# Plot SHAP value distributions
shap_plot <- plot_shapley_distributions(xgb_results, metric = "mean_abs_shap")

# Calculate specification interactions
interaction_results <- calculate_specification_interactions(
  xgb_results, 
  dataset_name = "homicide_analysis"
)

# Plot interaction heatmap
interaction_plot <- plot_interaction_heatmap(interaction_results)

# =============================================================================
# STEP 4: SAVE ALL RESULTS
# =============================================================================

# Save plots
plots <- list(
  spec_curve = plot_basic,
  shap_distributions = shap_plot,
  interactions = interaction_plot
)

save_visualization_plots(plots, "results/analysis_plots")

# Print summary
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Specification curve results saved to: results/spec_curve_results.csv\n")
cat("XGBoost + SHAP results available in xgb_results object\n")
cat("Plots saved to: results/analysis_plots_*.pdf\n")

# =============================================================================
# STEP 5: INTERPRETATION HELPERS
# =============================================================================

# Show most important features for treated unit
if (nrow(xgb_results$results) > 0) {
  treated_features <- xgb_results$results[unit == "PAPEP0000"]
  if (nrow(treated_features) > 0) {
    cat("\nMost important specification features for treated unit:\n")
    top_features <- treated_features[order(-mean_abs_shap)][1:5]
    print(top_features[, .(feature, mean_abs_shap, train_correlation, test_correlation)])
  }
}

# Show variance explained by individual features
if (!is.null(interaction_results$individual_summary)) {
  cat("\nVariance in treatment effects explained by specification choices:\n")
  print(interaction_results$individual_summary[order(-variance_explained_pct)])
}

cat("\n=== END PIPELINE ===\n")