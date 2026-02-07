# Specification Curve + SHAP Example
#
# Compact example aligned with the README workflow.

library(SCMs)

set.seed(42)
d <- expand.grid(
  unit = c("treated", "c1", "c2", "c3", "c4"),
  year = 2001:2010,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Keep IDs simple/alphanumeric for robust matching across internal reshaping paths.
d$unit_id <- gsub("[^A-Za-z0-9_]", "_", d$unit)

d$outcome <-
  100 + 0.5 * (d$year - 2001) +
  ifelse(d$unit == "treated" & d$year >= 2007, 2, 0) +
  rnorm(nrow(d), 0, 0.3)

d$population <-
  50 + 0.1 * (d$year - 2001) +
  rnorm(nrow(d), 0, 0.2)

spec_results <- spec_curve(
  dataset = d,
  outcomes = "outcome",
  col_name_unit_name = "unit_id",
  name_treated_unit = "treated",
  covagg = list(
    "Outcome Per Period" = list(
      label = "Outcome Per Period",
      per_period = "outcome_var"
    ),
    "Outcome Per Period + Pop Mean" = list(
      label = "Outcome Per Period + Pop Mean",
      per_period = "outcome_var",
      pre_period_mean = "population"
    )
  ),
  treated_period = 2007,
  min_period = 2001,
  end_period = 2010,
  col_name_period = "year",
  feature_weights = c("uniform", "optimize"),
  outcome_models = c("none", "ridge"),
  donor_sample = "all",
  constraints = list(list(name = "simplex")),
  constants = FALSE,
  cores = 1,
  verbose = FALSE,
  inference_type = "placebo"
)

xgb_cfg <- create_xgboost_config(
  dataset_name = "toy",
  treated_unit_name = "treated",
  outcome_filter = "outcome",
  spec_features = c("feat", "outcome_model", "const", "fw", "data_sample")
)

shap_results <- run_xgboost_shap_analysis(
  spec_results$results,
  xgb_cfg,
  compute_loo = FALSE
)

plot_obj <- plot_spec_curve(
  long_data = spec_results,
  name_treated_unit = "treated",
  outcomes = "outcome",
  shap_values = shap_results$shapley,
  show_pvalues = TRUE,
  test_statistic = "treatment_effect"
)

# Display in interactive sessions:
# print(plot_obj$final_plot)
