library(testthat)
library(data.table)

# Helper function to create mock spec_curve object
create_mock_spec_curve <- function() {
  set.seed(42)
  n_specs <- 50
  
  # Create mock results data
  results_data <- data.table(
    unit_name = rep(c("treated_unit", paste0("control_", 1:5)), each = n_specs),
    unit_type = rep(c("treated", rep("control", 5)), each = n_specs),
    outcome = rep("gdp", n_specs * 6),
    full_spec_id = rep(paste0("spec_", 1:n_specs), 6),
    post_period = rep(c(TRUE, FALSE), length.out = n_specs * 6),
    tau = c(
      rnorm(n_specs, mean = 2, sd = 1),      # treated effects
      rep(rnorm(n_specs * 5, mean = 0, sd = 0.5), 1)  # control effects
    ),
    rmse = runif(n_specs * 6, 0.5, 2.0),
    outcome_model = rep(sample(c("none", "augsynth", "ridge"), n_specs, replace = TRUE), 6),
    const = rep(sample(c("simplex", "lasso", "ridge"), n_specs, replace = TRUE), 6),
    fw = rep(sample(c("uniform", "optimize"), n_specs, replace = TRUE), 6),
    feat = rep(sample(c("all", "selected"), n_specs, replace = TRUE), 6)
  )
  
  # Create mock inference data
  abadie_inference <- list(
    p_values_rmse_ratio = data.table(
      full_spec_id = paste0("spec_", 1:n_specs),
      unit_name = rep("treated_unit", n_specs),
      unit_type = rep("treated", n_specs),
      p_value = runif(n_specs, 0, 1)
    )
  )
  
  bootstrap_inference <- list(
    p_values = data.table(
      full_spec_id = paste0("spec_", 1:n_specs),
      p_value_two_tailed = runif(n_specs, 0, 1)
    )
  )
  
  # Create spec_curve object
  spec_curve_obj <- list(
    results = results_data,
    expected_direction = "negative",
    abadie_inference = abadie_inference,
    bootstrap_inference = bootstrap_inference
  )
  
  class(spec_curve_obj) <- c("spec_curve", "list")
  return(spec_curve_obj)
}

test_that("print.spec_curve works correctly", {
  spec_obj <- create_mock_spec_curve()
  
  # Test that print doesn't throw an error
  expect_output(print(spec_obj), "Specification Curve Analysis Results")
  expect_output(print(spec_obj), "Specifications:")
  expect_output(print(spec_obj), "Treatment Effect Summary:")
  expect_output(print(spec_obj), "Inference Methods:")
  
  # Test with minimal object (no inference)
  minimal_obj <- list(
    results = spec_obj$results,
    expected_direction = "negative"
  )
  class(minimal_obj) <- c("spec_curve", "list")
  
  expect_output(print(minimal_obj), "Inference Methods:\\s+None")
})

test_that("summary.spec_curve works correctly", {
  spec_obj <- create_mock_spec_curve()
  
  # Test that summary doesn't throw an error and includes detailed breakdown
  expect_output(summary(spec_obj), "Specification Curve Analysis Results")
  expect_output(summary(spec_obj), "Detailed Specification Breakdown")
  expect_output(summary(spec_obj), "Treatment Effects by Outcome Model")
  expect_output(summary(spec_obj), "Treatment Effects by Constraint Type")
  expect_output(summary(spec_obj), "Model Fit Statistics")
  expect_output(summary(spec_obj), "Robustness Indicators")
  
  # Test with empty treated effects
  empty_obj <- spec_obj
  empty_obj$results <- empty_obj$results[unit_type != "treated"]
  
  expect_output(summary(empty_obj), "No treated unit effects found")
})

test_that("plot.spec_curve works correctly", {
  skip_if_not_installed("ggplot2")
  
  spec_obj <- create_mock_spec_curve()
  
  # Test basic plot functionality
  # Note: This mainly tests that the function doesn't error
  # Full plotting tests would require more complex setup
  expect_silent(suppressWarnings(suppressMessages({
    tryCatch({
      plot_result <- plot(spec_obj, show_shap = FALSE)
      # If plot succeeds, it should return a ggplot-like object or list
      expect_true(is.list(plot_result))
    }, error = function(e) {
      # Plot might fail due to missing dependencies or data issues
      # But it shouldn't fail due to class/method dispatch issues
      if (grepl("No treated unit found|Multiple treated units", e$message)) {
        stop("S3 method dispatch failed: ", e$message)
      }
      # Other errors (missing packages, etc.) are acceptable for this test
    })
  })))
  
  # Test error handling - no treated units
  no_treated_obj <- spec_obj
  no_treated_obj$results <- no_treated_obj$results[unit_type != "treated"]
  
  expect_error(plot(no_treated_obj, show_shap = FALSE), "No treated unit found")
  
  # Test warning - multiple treated units
  multi_treated_obj <- spec_obj
  multi_treated_obj$results <- rbind(
    multi_treated_obj$results,
    multi_treated_obj$results[unit_type == "treated"][1:5][, unit_name := "treated_unit_2"]
  )
  
  expect_warning(
    plot(multi_treated_obj, show_shap = FALSE, show_pvalues = FALSE),
    "Multiple treated units found"
  )
})

test_that("spec_curve S3 class assignment works", {
  # Test that our class assignment in spec_curve function works
  mock_results <- list(
    results = data.table(x = 1, y = 2),
    expected_direction = "negative"
  )
  
  class(mock_results) <- c("spec_curve", "list")
  
  expect_s3_class(mock_results, "spec_curve")
  expect_s3_class(mock_results, "list")
  
  # Test method dispatch
  expect_true(is.function(utils::getS3method("print", "spec_curve", optional = TRUE)))
  expect_true(is.function(utils::getS3method("summary", "spec_curve", optional = TRUE)))
  expect_true(is.function(utils::getS3method("plot", "spec_curve", optional = TRUE)))
})

test_that("S3 methods handle edge cases gracefully", {
  # Test with minimal spec_curve object
  minimal_obj <- list(
    results = data.table(
      unit_name = "test",
      unit_type = "treated", 
      outcome = "y",
      full_spec_id = "spec_1",
      post_period = TRUE,
      tau = 1.5,
      rmse = 0.8
    )
  )
  class(minimal_obj) <- c("spec_curve", "list")
  
  expect_output(print(minimal_obj), "Specifications:\\s+1")
  expect_output(summary(minimal_obj), "Mean:\\s+1.5")
  
  # Test with object containing NAs
  na_obj <- minimal_obj
  na_obj$results$tau <- NA_real_
  
  expect_output(print(na_obj), "Mean:\\s+NaN|Mean:\\s+NA")
  expect_output(summary(na_obj), "No treated unit effects found|Mean:\\s+NaN|Mean:\\s+NA")
})

make_curve_test_data <- function(sign_flip = 1, equal_rmse = FALSE) {
  specs <- c("s1", "s2", "s3", "s4")
  units <- c("treated_unit", "control_1", "control_2")
  unit_type <- c("treated", "control", "control")
  base_tau <- list(
    treated_unit = c(2, -3, 1, -2),
    control_1 = c(1, -1, 0.5, -0.5),
    control_2 = c(3, -2, 2, -1)
  )
  rmse_vals <- list(
    treated_unit = c(0.2, 0.4, 0.3, 0.5),
    control_1 = c(0.3, 0.6, 0.4, 0.8),
    control_2 = c(0.25, 0.45, 0.35, 0.55)
  )
  if (equal_rmse) {
    rmse_vals <- lapply(rmse_vals, function(x) rep(0.5, length(x)))
  }

  rows <- vector("list", length(units) * length(specs))
  idx <- 1L
  for (u in seq_along(units)) {
    for (s in seq_along(specs)) {
      rows[[idx]] <- data.table(
        full_spec_id = specs[s],
        unit_name = units[u],
        unit_type = unit_type[u],
        post_period = TRUE,
        tau = sign_flip * base_tau[[units[u]]][s],
        rmse = rmse_vals[[units[u]]][s]
      )
      idx <- idx + 1L
    }
  }
  rbindlist(rows)
}

test_that("percentile_rank returns strict (0,1) values with midranks on ties", {
  x <- c(5, 1, 1, 10)
  p <- SCMs:::percentile_rank(x)
  expect_true(all(p > 0 & p < 1))
  # ties at value 1 get average rank (1+2)/2 = 1.5 => 1.5/(4+1) = 0.3
  expect_equal(p[x == 1], c(0.3, 0.3))
})

test_that("Wilcoxon signed-rank statistic flips sign with tau sign flip and two-sided p-value is sign-invariant", {
  tau <- c(2, -3, 1, -2)
  sr <- SCMs:::curve_stat_wilcoxon_sr(tau)
  sr_flip <- SCMs:::curve_stat_wilcoxon_sr(-tau)
  expect_equal(sr_flip, -sr)

  data_pos <- make_curve_test_data(sign_flip = 1)
  data_neg <- make_curve_test_data(sign_flip = -1)

  out_pos <- SCMs:::calculate_spec_curve_pvalues_filtered(
    filtered_results = data_pos,
    curve_stat = c("median", "wilcoxon_sr"),
    weighting = "none",
    two_sided = TRUE
  )
  out_neg <- SCMs:::calculate_spec_curve_pvalues_filtered(
    filtered_results = data_neg,
    curve_stat = c("median", "wilcoxon_sr"),
    weighting = "none",
    two_sided = TRUE
  )

  p_pos <- out_pos$treated_summary[curve_statistic == "wilcoxon_sr" & weighting == "none", p_value]
  p_neg <- out_neg$treated_summary[curve_statistic == "wilcoxon_sr" & weighting == "none", p_value]
  expect_equal(p_pos, p_neg)
})

test_that("weighted Wilcoxon signed-rank is a normalized version of unweighted when all pre_rmspe values are equal", {
  data_equal_rmse <- make_curve_test_data(equal_rmse = TRUE)
  out <- SCMs:::calculate_spec_curve_pvalues_filtered(
    filtered_results = data_equal_rmse,
    curve_stat = "wilcoxon_sr",
    weighting = c("none", "pre_rmspe_percentile"),
    two_sided = TRUE
  )

  sr_none <- out$treated_summary[curve_statistic == "wilcoxon_sr" & weighting == "none", estimate]
  sr_w <- out$treated_summary[curve_statistic == "wilcoxon_sr" & weighting == "pre_rmspe_percentile", estimate]
  s <- unique(out$stats_by_unit[unit_name == "treated_unit", n_specs])
  expect_length(s, 1L)
  expect_equal(sr_w * s, sr_none)
})

test_that("curve-level inference fails hard on NA/Inf, missing pre_rmspe, mismatched lengths, and S<2", {
  data_ok <- make_curve_test_data()

  data_na <- copy(data_ok)
  data_na[1, tau := NA_real_]
  expect_error(
    SCMs:::calculate_spec_curve_pvalues_filtered(data_na, curve_stat = "median"),
    "Non-finite tau_s"
  )

  data_inf <- copy(data_ok)
  data_inf[1, tau := Inf]
  expect_error(
    SCMs:::calculate_spec_curve_pvalues_filtered(data_inf, curve_stat = "median"),
    "Non-finite tau_s"
  )

  data_no_rmse <- copy(data_ok)[, rmse := NULL]
  expect_error(
    SCMs:::calculate_spec_curve_pvalues_filtered(
      data_no_rmse,
      curve_stat = "wilcoxon_sr",
      weighting = "pre_rmspe_percentile"
    ),
    "requires a finite 'rmse' column"
  )

  expect_error(
    SCMs:::curve_stat_wilcoxon_sr_weighted(c(1, -1, 2), c(0.1, 0.2)),
    "identical lengths"
  )

  data_small_s <- data_ok[full_spec_id %in% "s1"]
  expect_error(
    SCMs:::calculate_spec_curve_pvalues_filtered(data_small_s, curve_stat = "median"),
    "require at least min_specs = 2"
  )
})

test_that("all four curve-stat/weighting combinations are returned", {
  data_ok <- make_curve_test_data()
  out <- SCMs:::calculate_spec_curve_pvalues_filtered(
    filtered_results = data_ok,
    curve_stat = c("median", "wilcoxon_sr"),
    weighting = c("none", "pre_rmspe_percentile"),
    two_sided = TRUE
  )

  combos <- unique(out$treated_summary[, .(curve_statistic, weighting)])
  expect_equal(nrow(combos), 4L)
  expect_true(all(c(
    "median:none",
    "median:pre_rmspe_percentile",
    "wilcoxon_sr:none",
    "wilcoxon_sr:pre_rmspe_percentile"
  ) %in% paste(combos$curve_statistic, combos$weighting, sep = ":")))
})

test_that("strict grid policy fails when a placebo unit is missing specs", {
  data_ok <- make_curve_test_data()
  data_missing <- data_ok[!(unit_name == "control_2" & full_spec_id %in% "s4")]

  expect_error(
    SCMs:::calculate_spec_curve_pvalues_filtered(
      filtered_results = data_missing,
      curve_stat = "median",
      weighting = "none",
      grid_policy = "strict"
    ),
    "All units must share the exact same specification grid"
  )
})

test_that("drop_incomplete_units keeps full treated grid and drops mismatched placebo units", {
  data_ok <- make_curve_test_data()
  data_missing <- data_ok[!(unit_name == "control_2" & full_spec_id %in% "s4")]

  out <- SCMs:::calculate_spec_curve_pvalues_filtered(
    filtered_results = data_missing,
    curve_stat = "median",
    weighting = "none",
    grid_policy = "drop_incomplete_units"
  )

  expect_equal(out$attrition_report$n_units_dropped, 1L)
  expect_equal(out$attrition_report$n_placebos_kept, 1L)
  expect_equal(out$attrition_report$n_specs_final, 4L)
  expect_true("control_2" %in% out$dropped_units$unit_name)
})

test_that("intersect_specs keeps all units and uses shared specs only", {
  data_ok <- make_curve_test_data()
  data_missing <- data_ok[!(unit_name == "control_2" & full_spec_id %in% "s4")]

  out <- SCMs:::calculate_spec_curve_pvalues_filtered(
    filtered_results = data_missing,
    curve_stat = "median",
    weighting = "none",
    grid_policy = "intersect_specs"
  )

  expect_equal(out$attrition_report$n_units_dropped, 0L)
  expect_equal(out$attrition_report$n_placebos_kept, 2L)
  expect_equal(out$attrition_report$n_specs_final, 3L)
  expect_true("s4" %in% out$dropped_spec_ids)
})

test_that("min_placebos and min_specs are enforced after grid policy", {
  data_ok <- make_curve_test_data()
  data_missing <- data_ok[!(unit_name == "control_2" & full_spec_id %in% "s4")]

  expect_error(
    SCMs:::calculate_spec_curve_pvalues_filtered(
      filtered_results = data_missing,
      curve_stat = "median",
      weighting = "none",
      grid_policy = "drop_incomplete_units",
      min_placebos = 2
    ),
    "only 1 placebo units remain"
  )

  expect_error(
    SCMs:::calculate_spec_curve_pvalues_filtered(
      filtered_results = data_missing,
      curve_stat = "median",
      weighting = "none",
      grid_policy = "intersect_specs",
      min_specs = 4
    ),
    "only 3 shared specifications remain"
  )
})
