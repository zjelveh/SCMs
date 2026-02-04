#' German Reunification Analysis - Complete SCMs Workflow Example
#' 
#' This example demonstrates the full SCMs package workflow using the German 
#' reunification case study from Abadie, Diamond, and Hainmueller (2015).
#' The analysis examines the economic impact of German reunification on 
#' West Germany's GDP per capita.
#' 
#' @author Zubin Jelveh
#' @references Abadie, Alberto, Alexis Diamond, and Jens Hainmueller. 
#'   "Comparative politics and the synthetic control method." American Journal 
#'   of Political Science 59.2 (2015): 495-510.

library(SCMs)
library(data.table)

# ===== DATA LOADING =====

# Load the German reunification dataset
# This dataset contains GDP per capita and covariates for West Germany and 
# potential donor countries from 1960-2003
dataset <- fread(system.file("extdata/scpi_germany.csv", package = "SCMs"))

# Examine the data structure
cat("Dataset dimensions:", nrow(dataset), "rows x", ncol(dataset), "columns\n")
cat("Countries included:", length(unique(dataset$country)), "\n")
cat("Time period:", min(dataset$year), "-", max(dataset$year), "\n")
cat("Treated unit: West Germany\n")
cat("Treatment period: 1991 (German reunification)\n\n")

# ===== SPECIFICATION CURVE ANALYSIS SETUP =====

# Define comprehensive specification parameters
# This creates a large specification space to test robustness
params <- list(
  cores = 4,  # Adjust based on your system
  outcomes = "gdp",
  col_name_unit_name = "country",
  name_treated_unit = "West Germany",
  
  # Multiple covariate specifications following different approaches
  covagg = list(
    # Specification 1: ADH-style with all variables
    adh_full = list(
      label = 'ADH Full',
      gdp_each = list(var = "gdp", each = TRUE),
      invest_final = list(var = "invest80", periods = 1980),
      schooling_census = list(var = "schooling", every_n = 5),
      trade_avg = list(var = "trade", average = "full_pre"),
      industry_final = list(var = "industry", periods = 1980)
    ),
    
    # Specification 2: Parsimonious approach
    parsimonious = list(
      label = 'GDP + Trade Average',
      gdp_avg = list(var = "gdp", average = "full_pre"),
      trade_ave = list(var = 'trade', average = 'full_pre')
    ),
    
    # Specification 3: GDP-only matching
    gdp_only = list(
      label = 'GDP Only',
      gdp_each = list(var = "gdp", each = TRUE)
    )
  ),
  
  # Treatment parameters
  treated_period = 1991,  # German reunification
  min_period = 1960,
  end_period = 2003,
  col_name_period = "year",
  
  # Specification curve dimensions
  feature_weights = c("uniform", "optimize"),
  donor_sample = c("all", "most_similar"),
  outcome_models = c("none", "augsynth", "lasso", "ridge", "ols"),
  constraints = list(
    list(name = "simplex"),
    list(name = "lasso"),
    list(name = "ridge")
  ),
  
  # Inference configuration
  inference_type = 'placebo',
  inference_config = list(
    bootstrap_n_replications = 100,
    verbose = TRUE
  )
)

cat("Specification space size:", 
    length(params$outcomes) *
    length(params$covagg) *
    length(params$feature_weights) *
    length(params$donor_sample) *
    length(params$outcome_models) *
    length(params$constraints), "specifications\n\n")

# ===== RUN SPECIFICATION CURVE ANALYSIS =====

cat("Running specification curve analysis...\n")
cat("This may take several minutes depending on cores and specification space size.\n\n")

# Run the complete specification curve analysis
results_germany <- run_spec_curve_analysis(
  dataset = dataset, 
  params = params
)

cat("Specification curve analysis complete!\n")
cat("Results structure:\n")
cat("- Main results:", nrow(results_germany$results), "rows\n")
cat("- Placebo inference available:", !is.null(results_germany$abadie_inference), "\n")
cat("- Test statistics:", if(!is.null(results_germany$abadie_inference)) 
    names(results_germany$abadie_inference) else "None", "\n\n")

# ===== XGBOOST + SHAP ANALYSIS =====

cat("Running XGBoost + SHAP analysis for interpretability...\n")

# Create XGBoost analysis configuration
xgb_config <- create_xgboost_config(
  dataset_name = "german_reunification_analysis",
  treated_unit_name = "West Germany",
  outcome_filter = "gdp",
  spec_features = c('feat', 'outcome_model', 'const', 'data_sample', 'fw'),
  treated_unit_only = TRUE
)

# Run XGBoost model with SHAP values
germany_shap <- run_xgboost_shap_analysis(results_germany$results, xgb_config)

cat("XGBoost analysis complete!\n")
cat("SHAP results available for", nrow(germany_shap$shapley), "specifications\n\n")

# ===== VISUALIZATION =====

cat("Creating specification curve visualizations...\n")

# Plot 1: RMSE ratio test statistic with SHAP coloring
plot_rmse <- plot_spec_curve(
  pvalue_style = 'continuous',
  long_data = results_germany,
  show_pvalues = TRUE, 
  prefer_bootstrap_pvalues = FALSE,
  name_treated_unit = 'West Germany', 
  outcomes = 'gdp', 
  normalize_outcomes = 'none',
  shap_values = germany_shap$shapley,
  null_distribution = 'placebo',
  test_statistic = 'rmse_ratio'  # Traditional Abadie approach
)

print(plot_rmse)

# Plot 2: Treatment effect test statistic  
plot_te <- plot_spec_curve(
  pvalue_style = 'continuous',
  long_data = results_germany,
  show_pvalues = TRUE, 
  name_treated_unit = 'West Germany', 
  outcomes = 'gdp', 
  normalize_outcomes = 'none',
  shap_values = germany_shap$shapley,
  null_distribution = 'placebo',
  test_statistic = 'treatment_effect'  # Direct treatment effect
)

print(plot_te)

# Plot 3: Normalized treatment effect test statistic
plot_norm_te <- plot_spec_curve(
  pvalue_style = 'categorical',
  long_data = results_germany,
  show_pvalues = TRUE, 
  name_treated_unit = 'West Germany', 
  outcomes = 'gdp', 
  normalize_outcomes = 'none',
  shap_values = germany_shap$shapley,
  null_distribution = 'placebo',
  test_statistic = 'normalized_te'  # Treatment effect / pre-period SD
)

print(plot_norm_te)

# ===== RESULTS SUMMARY =====

cat("\n===== ANALYSIS SUMMARY =====\n")

# Extract key statistics
treated_specs <- results_germany$results[unit_name == "West Germany" & post_period == TRUE]
avg_effects <- treated_specs[, .(avg_tau = mean(tau, na.rm = TRUE)), by = full_spec_id]

cat("Treatment effect statistics:\n")
cat("- Mean effect across specifications:", round(mean(avg_effects$avg_tau, na.rm = TRUE), 3), "\n")
cat("- Median effect:", round(median(avg_effects$avg_tau, na.rm = TRUE), 3), "\n")
cat("- Standard deviation:", round(sd(avg_effects$avg_tau, na.rm = TRUE), 3), "\n")
cat("- Min effect:", round(min(avg_effects$avg_tau, na.rm = TRUE), 3), "\n")
cat("- Max effect:", round(max(avg_effects$avg_tau, na.rm = TRUE), 3), "\n")

# P-value summary across test statistics
if (!is.null(results_germany$abadie_inference)) {
  test_stats <- c("rmse_ratio", "treatment_effect", "normalized_te")
  
  cat("\nSignificance summary (p < 0.05):\n")
  for (stat in test_stats) {
    p_key <- paste0("p_values_", stat)
    if (!is.null(results_germany$abadie_inference[[p_key]])) {
      p_vals <- results_germany$abadie_inference[[p_key]]
      sig_rate <- mean(p_vals$p_value < 0.05, na.rm = TRUE)
      cat("- ", stat, ": ", round(sig_rate * 100, 1), "% of specifications significant\n")
    }
  }
}

cat("\nAnalysis complete! The specification curve shows the robustness of\n")
cat("the estimated impact of German reunification across different modeling choices.\n")
cat("SHAP values indicate which specification choices are most important\n")
cat("for predicting the magnitude of treatment effects.\n")