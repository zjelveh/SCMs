library(testthat)

create_plot_spec_curve_results <- function() {
  set.seed(42)
  d <- expand.grid(
    unit = c("treated", "c1", "c2", "c3", "c4"),
    year = 2001:2010,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  d$outcome <- 100 + 0.5 * (d$year - 2001) +
    ifelse(d$unit == "treated" & d$year >= 2007, 2, 0) +
    rnorm(nrow(d), 0, 0.3)
  spec_curve(
    dataset = d,
    outcomes = "outcome",
    col_name_unit_name = "unit",
    name_treated_unit = "treated",
    covagg = list(
      "Outcome Per Period" = list(
        label = "Outcome Per Period",
        operations = list(
          list(var = "outcome_var", partition_periods = list(type = "by_period"))
        )
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
}

test_that("plot_spec_curve fails hard when SHAP values do not cover the treated grid", {
  skip_if_not_installed("CVXR")

  spec_results <- suppressMessages(create_plot_spec_curve_results())
  treated_ids <- unique(spec_results$results[unit_name == "treated", full_spec_id])
  expect_gt(length(treated_ids), 1)

  plot_with_shap <- function(shap) {
    plot_spec_curve(
      long_data = spec_results,
      name_treated_unit = "treated",
      outcomes = "outcome",
      show_shap = TRUE,
      shap_values = shap,
      test_statistic = "treatment_effect"
    )
  }

  partial_shap <- data.table::data.table(
    unit = "treated",
    full_spec_id = treated_ids[1],
    feature_group = "fw",
    feature = "uniform",
    shapley_value = 0.1
  )
  expect_error(suppressMessages(plot_with_shap(partial_shap)), "did not merge onto Panel B")

  other_unit_shap <- data.table::copy(partial_shap)[, unit := "c1"]
  expect_error(suppressMessages(plot_with_shap(other_unit_shap)), "No SHAP values found for treated unit")
})
