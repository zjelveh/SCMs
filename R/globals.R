#' Internal Imports and NSE Globals
#'
#' Centralizes imports used via non-qualified calls and declares
#' non-standard-evaluation symbols to satisfy R CMD check.
#'
#' @name scm-globals
#' @keywords internal
#' @noRd
#' @importFrom grDevices dev.off pdf
#' @importFrom stats model.matrix
#' @importFrom utils getS3method
#' @importFrom R6 R6Class
NULL

utils::globalVariables(
  c(
    ".", ".data", ".fitted", ".resid", "..available_cols", "..features_to_use",
    "ID", "Legend", "Specification", "Time", "Unit Name", "actual", "avg_tau",
    "avg_treatment_effect", "categorical_levels", "const", "constant",
    "curve_statistic", "data_sample", "diff_with_treated", "estimate", "Estimate",
    "feat", "feature", "feature_display", "feature_display_chr",
    "feature_display_ordered", "feature_group", "feature_level",
    "feature_with_shap", "full_spec_id", "fw", "grid_matches_treated",
    "is_categorical", "is_significant", "mean_abs_shap", "n_unique",
    "n_unique_rmse", "new_specification", "onehot_col", "original_unit_name",
    "outcome", "outcome_model", "p.value", "p_value", "p_value_two_tailed",
    "post_period", "post_pre_ratio", "pre_rmse", "pre_rmspe", "Predicted",
    "predicted_loo", "rmse", "rnk_abs_mean", "shap_color", "shap_color_value",
    "shap_summary", "shapley_value", "significance_category", "significant",
    "spec_count", "spec_ids", "spec_number", "stat_rank", "tau", "tau_s",
    "test_statistic_value", "trt", "unit", "unit_name", "unit_numbers",
    "unit_type", "weight_rank", "x", "y", "y_numeric", "ybar_control"
  )
)
