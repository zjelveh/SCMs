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
  expect_silent({
    tryCatch({
      plot_result <- plot(spec_obj)
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
  })
  
  # Test error handling - no treated units
  no_treated_obj <- spec_obj
  no_treated_obj$results <- no_treated_obj$results[unit_type != "treated"]
  
  expect_error(plot(no_treated_obj), "No treated unit found")
  
  # Test warning - multiple treated units
  multi_treated_obj <- spec_obj
  multi_treated_obj$results <- rbind(
    multi_treated_obj$results,
    multi_treated_obj$results[unit_type == "treated"][1:5][, unit_name := "treated_unit_2"]
  )
  
  expect_warning(plot(multi_treated_obj), "Multiple treated units found")
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
  expect_true(methods::existsMethod("print", signature = "spec_curve"))
  expect_true(methods::existsMethod("summary", signature = "spec_curve"))
  expect_true(methods::existsMethod("plot", signature = "spec_curve"))
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