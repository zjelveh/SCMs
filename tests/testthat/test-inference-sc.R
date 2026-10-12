library(testthat)

create_inference_mock_data <- function() {
  set.seed(321)
  units <- c("treated", "c1", "c2", "c3")
  years <- 2001:2010
  data <- expand.grid(unit = units, year = years, stringsAsFactors = FALSE)

  data$outcome <- 50 + 0.4 * (data$year - 2001) +
    ifelse(data$unit == "treated" & data$year >= 2007, 2.5, 0) +
    rnorm(nrow(data), 0, 0.3)
  data$cov1 <- 20 + 0.2 * (data$year - 2001) + rnorm(nrow(data), 0, 0.2)
  data
}

test_that("inference_sc returns placebo inference from estimate_sc result", {
  skip_if_not_installed("CVXR")

  mock_data <- create_inference_mock_data()
  covagg <- list(
    list(var = "outcome_var", partition_periods = list(type = "by_period")),
    list(var = "cov1", compute = "mean")
  )

  sc_result <- estimate_sc(
    dataset = mock_data,
    outcome = "outcome",
    covagg = covagg,
    col_name_unit_name = "unit",
    name_treated_unit = "treated",
    col_name_period = "year",
    treated_period = 2007,
    min_period = 2001,
    end_period = 2010,
    outcome_models = "none",
    feature_weights = "uniform",
    w.constr = list(name = "simplex")
  )

  inf <- inference_sc(
    sc.pred = sc_result,
    dataset = mock_data,
    cores = 1,
    verbose = FALSE
  )

  expect_s3_class(inf, "scm_inference")
  expect_true(all(c("inference.results", "rmse", "abadie_significance") %in% names(inf)))
  expect_gt(nrow(inf$inference.results), 0)
})

test_that("inference_placebo returns expected structure and test-stat blocks", {
  skip_if_not_installed("CVXR")

  mock_data <- create_inference_mock_data()
  covagg <- list(
    list(var = "outcome_var", partition_periods = list(type = "by_period")),
    list(var = "cov1", compute = "mean")
  )

  sc_result <- estimate_sc(
    dataset = mock_data,
    outcome = "outcome",
    covagg = covagg,
    col_name_unit_name = "unit",
    name_treated_unit = "treated",
    col_name_period = "year",
    treated_period = 2007,
    min_period = 2001,
    end_period = 2010,
    outcome_models = "none",
    feature_weights = "uniform",
    w.constr = list(name = "simplex")
  )

  placebo <- inference_placebo(
    sc.pred = sc_result,
    dataset = mock_data,
    cores = 1,
    verbose = FALSE,
    expected_direction = "two_sided"
  )

  expect_true(all(c("taus", "rmse", "abadie_significance") %in% names(placebo)))
  expect_true(data.table::is.data.table(placebo$taus))
  expect_true(data.table::is.data.table(placebo$rmse))
  expect_true("treated" %in% placebo$taus$unit_name)
  expect_true("treated" %in% placebo$rmse$unit_name)

  abadie_blocks <- c("rmse_ratio", "treatment_effect", "normalized_te")
  expect_true(all(abadie_blocks %in% names(placebo$abadie_significance)))
  expect_true(all(vapply(
    abadie_blocks,
    function(b) data.table::is.data.table(placebo$abadie_significance[[b]]$p_values),
    logical(1)
  )))
})

test_that("inference_placebo fails hard on missing metadata and bad dataset column", {
  skip_if_not_installed("CVXR")

  mock_data <- create_inference_mock_data()
  covagg <- list(
    list(var = "outcome_var", partition_periods = list(type = "by_period")),
    list(var = "cov1", compute = "mean")
  )

  sc_result <- estimate_sc(
    dataset = mock_data,
    outcome = "outcome",
    covagg = covagg,
    col_name_unit_name = "unit",
    name_treated_unit = "treated",
    col_name_period = "year",
    treated_period = 2007,
    min_period = 2001,
    end_period = 2010,
    outcome_models = "none",
    feature_weights = "uniform",
    w.constr = list(name = "simplex")
  )

  missing_meta <- sc_result
  missing_meta$col_name_unit_name <- NULL
  expect_error(
    inference_placebo(
      sc.pred = missing_meta,
      dataset = mock_data,
      cores = 1,
      verbose = FALSE
    ),
    "Missing col_name_unit_name or name_treated_unit"
  )

  bad_dataset <- data.table::copy(data.table::as.data.table(mock_data))
  bad_dataset[, unit := NULL]
  expect_error(
    inference_placebo(
      sc.pred = sc_result,
      dataset = bad_dataset,
      cores = 1,
      verbose = FALSE
    ),
    "Column unit not found in dataset"
  )
})

test_that("inference_sc fails hard on scest objects missing metadata", {
  skip_if_not_installed("CVXR")

  mock_data <- create_inference_mock_data()
  scm_data <- scdata(
    df = mock_data,
    id.var = "unit",
    time.var = "year",
    outcome.var = "outcome",
    period.pre = 2001:2006,
    period.post = 2007:2010,
    unit.tr = "treated",
    unit.co = c("c1", "c2", "c3"),
    covagg = list(list(var = "outcome_var", partition_periods = list(type = "by_period")))
  )
  fit <- scest(scm_data, w.constr = list(name = "simplex"))

  expect_error(
    inference_sc(sc.pred = fit, dataset = mock_data, cores = 1, verbose = FALSE),
    "missing required metadata"
  )
})

make_abadie_inputs <- function(treated_effect) {
  units <- c("treated", "c1", "c2", "c3")
  post_effects <- c(treated_effect, -0.5, 0.4, -1)
  taus <- data.table::data.table(
    unit_name = rep(units, each = 3),
    unit_type = rep(c("treated", "control", "control", "control"), each = 3),
    outcome_model = "none",
    post_period = rep(c(FALSE, FALSE, TRUE), times = 4),
    tau = as.vector(rbind(0.1, -0.1, post_effects))
  )
  rmse <- data.table::data.table(
    unit_name = units,
    unit_type = c("treated", "control", "control", "control"),
    outcome_model = "none",
    pre_rmse = 0.1
  )
  list(taus = taus, rmse = rmse)
}

test_that("two-sided placebo p-values rank units by absolute effect", {
  for (effect in c(3, -3)) {
    inputs <- make_abadie_inputs(effect)
    sig <- calculate_abadie_significance(
      inputs$taus, inputs$rmse, expected_direction = "two_sided"
    )
    p <- sig$treatment_effect$p_values
    expect_equal(p[unit_name == "treated", p_value], 0.25)
    expect_equal(p[unit_name == "c3", p_value], 0.50)
    expect_equal(p[unit_name == "c1", p_value], 0.75)
    expect_equal(p[unit_name == "c2", p_value], 1.00)
  }
})

test_that("calculate_abadie_significance rejects an unknown expected_direction", {
  inputs <- make_abadie_inputs(3)
  expect_error(
    calculate_abadie_significance(
      inputs$taus, inputs$rmse, expected_direction = "two-sided"
    ),
    "expected_direction must be one of"
  )
})
