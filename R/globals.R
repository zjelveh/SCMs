# Global variable declarations to avoid CRAN check notes
# These are mostly data.table column names and intermediate variables used in data operations

utils::globalVariables(c(
  # Data.table column bindings and data operations
  ".", ".data", "..available_cols", "..available_features", "..features_to_use",
  
  # Common data columns 
  "ID", "Time", "Unit Name", "Estimate", "Predicted", "Specification",
  "unit_name", "unit_type", "unit_numbers", "full_spec_id", "spec_number",
  "outcome", "outcome_model", "period", "post_period", "tau", "value",
  
  # Specification features
  "const", "constant", "feat", "feature", "fw", "data_sample",
  "feature_group", "feature_display", "feature_with_shap",
  "spec_combination",
  
  # SHAP and ML variables
  "shapley_value", "shap_summary", "shap_color", 
  "shap_color_value", "mean_shap", "mean_abs_shap",
  "abs_mean_interaction_contribution", "mean_abs_shap_interaction",
  "total_abs_shap_value", "categorical_levels", "category", "component",
  "interaction_label", "interaction_partner", "interaction_percent", 
  "interaction_ratio", "is_categorical", "main_effect_percent",
  "main_plus_interactions", "rnk_abs_mean",
  
  # Statistical inference variables
  "p_value", "p_value_one_sided", "p_value_two_tailed", "significance_category",
  "is_significant", "is_treated", "pval_rank_med", "pval_rank_z", 
  "rank_med", "rank_z", "zscore", "test_statistic_value",
  
  # Treatment effect variables
  "avg_tau", "ave_tau", "avg_treatment_effect", "median_tau", 
  "post_pre_ratio", "pre_rmse", "rmse", "predicted_loo",
  "diff_with_treated", "treated_unit", "mean_z", "unit",
  
  # Outcome and covariate processing
  "magnitude", "label", "new_specification", "original_feature",
  "original_unit_name", "ybar_control", "n_unique",
  
  # Dataset-specific variables (these may vary by application)
  "hr_rate", "num_homicide", "officers", "pop", "population", "year",
  "mdate", "modate", "month", "mh", "ofp", "ois", "oiso", "oisp", "ori9",
  "trt", "unarmed", "result", "read_xlsx",
  
  # Legacy function references
  "b.est.cvxr"
))
