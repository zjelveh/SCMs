#' Basque Country Terrorism Analysis - SCMs Package Example
#' 
#' This example replicates and extends the analysis from Abadie & Gardeazabal (2003)
#' examining the economic costs of terrorism in Spain's Basque Country using
#' specification curve analysis to demonstrate software package sensitivity.
#' 
#' @author Zubin Jelveh
#' @references Abadie, Alberto, and Javier Gardeazabal. 
#'   "The economic costs of conflict: A case study of the Basque Country." 
#'   American Economic Review 93.1 (2003): 113-132.

library(SCMs)
library(data.table)

# ===== DATA LOADING =====

# Load the Basque Country dataset (included with package)
dataset <- fread(system.file("extdata/basque.csv", package = "SCMs"))

# Examine the data structure
cat("Dataset dimensions:", nrow(dataset), "rows x", ncol(dataset), "columns\n")
cat("Regions included:", length(unique(dataset$regionname)), "\n")
cat("Time period:", min(dataset$year), "-", max(dataset$year), "\n")
cat("Treated unit: Basque Country (region 17)\n")
cat("Treatment period: 1970 (onset of terrorism)\n\n")

# ===== BASIC SYNTHETIC CONTROL ESTIMATION =====

# Run basic SCM estimation replicating Abadie & Gardeazabal (2003)
basic_results <- scest(
  dataset, 
  treated_unit = "Basque Country (Pais Vasco)",
  treated_period = 1970,
  outcome = "gdpcap",
  unit_col = "regionname",
  time_col = "year",
  # Use key covariates from original study
  covariates = c("popdens", "sec.agriculture", "sec.energy", 
                 "sec.industry", "sec.construction", "sec.services.venta")
)

# View basic results
summary(basic_results)
plot(basic_results)

# Run placebo inference
inference_results <- inference_sc(
  basic_results, 
  dataset, 
  inference_type = "placebo"
)

cat("Basic SCM Results:\n")
cat("Treatment Effect:", round(basic_results$treatment_effect, 3), "\n")
cat("P-value (placebo test):", round(inference_results$abadie_pvalue, 3), "\n\n")

# ===== SPECIFICATION CURVE ANALYSIS =====

# Define comprehensive covariate specifications
# This demonstrates how different variable construction approaches affect results
covariate_specs <- list(
  # Original Abadie & Gardeazabal approach: time-averaged
  original_avg = list(
    label = 'Time-averaged (Original)',
    gdp_baseline = list(var = "gdpcap", periods = c(1960, 1965)),
    popdens_avg = list(var = "popdens", average = "full_pre"),
    agriculture_avg = list(var = "sec.agriculture", average = "full_pre"),
    industry_avg = list(var = "sec.industry", average = "full_pre")
  ),
  
  # Per-period approach: each year separately  
  per_period = list(
    label = 'All per period',
    gdp_each = list(var = "gdpcap", each = TRUE),
    popdens_each = list(var = "popdens", each = TRUE),
    agriculture_each = list(var = "sec.agriculture", each = TRUE),
    industry_each = list(var = "sec.industry", each = TRUE)
  ),
  
  # Hybrid approach: outcome per period, covariates averaged
  hybrid = list(
    label = 'Outcome per period, covariates time-averaged',
    gdp_each = list(var = "gdpcap", each = TRUE),
    popdens_avg = list(var = "popdens", average = "full_pre"),
    agriculture_avg = list(var = "sec.agriculture", average = "full_pre"),
    industry_avg = list(var = "sec.industry", average = "full_pre")
  ),
  
  # Sectoral focus: emphasize economic structure
  sectoral_focus = list(
    label = 'Sectoral composition focus',
    gdp_baseline = list(var = "gdpcap", first = 5),
    all_sectors = list(var = "sec.agriculture", average = "full_pre"),
    energy = list(var = "sec.energy", average = "full_pre"), 
    industry = list(var = "sec.industry", average = "full_pre"),
    construction = list(var = "sec.construction", average = "full_pre"),
    services = list(var = "sec.services.venta", average = "full_pre")
  )
)

# Run comprehensive specification curve analysis
spec_results <- spec_curve(
  dataset = dataset,
  outcomes = "gdpcap",
  col_name_unit_name = "regionname",
  name_treated_unit = "Basque Country (Pais Vasco)",
  covagg = covariate_specs,
  treated_period = 1970,
  min_period = 1955,
  end_period = 1997,
  col_name_period = "year",
  
  # Test different methodological approaches
  feature_weights = c("uniform", "optimize"),
  outcome_models = c("none", "augsynth", "ridge", "lasso"),
  constraints = list(
    list(name = "simplex"),
    list(name = "lasso", Q = 0.1),
    list(name = "ridge", Q = 0.1)
  ),
  
  # Inference settings
  inference_type = "placebo",
  cores = 4
)

cat("Specification Curve Analysis completed with", 
    nrow(spec_results$results), "specifications\n")

# ===== SHAP INTERPRETABILITY ANALYSIS =====

# Configure machine learning analysis to understand what drives variation
catboost_config <- create_catboost_config(
  dataset_name = "basque_terrorism",
  treated_unit_name = "Basque Country (Pais Vasco)",
  outcome_filter = "gdpcap",
  spec_features = c("feat", "outcome_model", "const", "fw")
)

# Generate SHAP values to explain specification sensitivity
shap_results <- run_catboost_shap_analysis(
  spec_results$results, 
  catboost_config
)

# ===== VISUALIZATION AND RESULTS =====

# Create comprehensive specification curve with SHAP coloring
plot_spec_curve(
  spec_results,
  name_treated_unit = "Basque Country (Pais Vasco)",
  outcomes = "gdpcap",
  shap_values = shap_results$shapley,
  test_statistic = "treatment_effect",
  title = "Basque Country: Specification Sensitivity Analysis",
  show_shap = TRUE,
  file_path_save = "basque_spec_curve_analysis.pdf"
)

# Summary statistics
cat("\n===== SPECIFICATION CURVE SUMMARY =====\n")
results_summary <- spec_results$results[, .(
  min_effect = min(treatment_effect, na.rm = TRUE),
  max_effect = max(treatment_effect, na.rm = TRUE),
  median_effect = median(treatment_effect, na.rm = TRUE),
  mean_effect = mean(treatment_effect, na.rm = TRUE),
  sd_effect = sd(treatment_effect, na.rm = TRUE),
  pct_negative = mean(treatment_effect < 0, na.rm = TRUE) * 100,
  pct_significant = mean(pvalue_placebo < 0.05, na.rm = TRUE) * 100
)]

print(results_summary)

cat("\nKey Findings:\n")
cat("• Effect estimates range from", round(results_summary$min_effect, 2), 
    "to", round(results_summary$max_effect, 2), "\n")
cat("•", round(results_summary$pct_negative, 1), "% of specifications show negative effects\n")
cat("•", round(results_summary$pct_significant, 1), "% of specifications are statistically significant\n")

# Show which specification choices matter most
if (!is.null(shap_results$feature_importance)) {
  cat("\nMost important specification dimensions (by SHAP importance):\n")
  print(shap_results$feature_importance[order(-abs(mean_shap))][1:5])
}

cat("\n===== POLICY IMPLICATIONS =====\n")
cat("This analysis demonstrates that:\n")
cat("1. Single SCM estimates can be misleading due to software implementation choices\n")
cat("2. Specification curve analysis reveals the full range of defensible estimates\n") 
cat("3. SHAP analysis identifies which methodological choices drive result variation\n")
cat("4. Transparent reporting requires showing specification robustness\n")