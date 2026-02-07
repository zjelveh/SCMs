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
