library(testthat)
library(data.table)

# Helper function to create mock scest object for testing
create_mock_scest_object <- function(include_bootstrap = FALSE) {
  set.seed(123)
  
  n_donors <- 4
  n_pre <- 8  
  n_post <- 5
  
  # Create mock weights
  weights <- c(0.4, 0.3, 0.2, 0.1)
  names(weights) <- paste0("donor_", 1:n_donors)
  
  # Create mock fitted values
  Y_pre_fit <- rnorm(n_pre, mean = 50, sd = 5)
  Y_post_fit <- rnorm(n_post, mean = 45, sd = 5)
  
  # Create mock actual values
  Y_pre <- Y_pre_fit + rnorm(n_pre, 0, 2)
  Y_post <- Y_post_fit + rnorm(n_post, 0, 2)
  
  # Create mock constraint object
  w_constr <- list(name = "simplex")
  
  # Create mock V matrix
  V_matrix <- diag(n_pre)
  
  # Build estimation results
  est_results <- list(
    w = weights,
    Y.pre.fit = Y_pre_fit,
    Y.post.fit = Y_post_fit,
    w.constr = w_constr,
    V = V_matrix
  )
  
  # Add bootstrap results if requested
  if (include_bootstrap) {
    # Mock bootstrap weights (n_bootstrap x n_donors)
    n_bootstrap <- 100
    bootstrap_weights <- matrix(
      abs(rnorm(n_bootstrap * n_donors, rep(weights, each = n_bootstrap), 0.1)),
      nrow = n_bootstrap, ncol = n_donors
    )
    # Normalize to sum to 1
    bootstrap_weights <- bootstrap_weights / rowSums(bootstrap_weights)
    colnames(bootstrap_weights) <- names(weights)
    
    est_results$bootstrap_weights <- bootstrap_weights
  }
  
  # Create data component
  data_comp <- list(
    Y.pre = Y_pre,
    Y.post = Y_post,
    B = matrix(rnorm(n_pre * n_donors), nrow = n_pre, ncol = n_donors)
  )
  colnames(data_comp$B) <- names(weights)
  
  # Create scest object
  scest_obj <- list(
    est.results = est_results,
    data = data_comp
  )
  
  class(scest_obj) <- "scest"
  return(scest_obj)
}

# Helper function to create mock spec_curve object
create_mock_spec_curve_object <- function(include_inference = FALSE) {
  set.seed(456)
  
  n_specs <- 20
  
  # Create mock results data
  results <- data.table(
    full_spec_id = paste0("spec_", 1:n_specs),
    unit_name = rep("treated_unit", n_specs),
    unit_type = rep("treated", n_specs),
    outcome = rep("gdp", n_specs),
    post_period = rep(TRUE, n_specs),
    tau = rnorm(n_specs, -5, 2),
    rmse = runif(n_specs, 0.5, 2.0),
    outcome_model = sample(c("none", "ridge", "lasso"), n_specs, replace = TRUE),
    const = sample(c("simplex", "lasso", "ridge"), n_specs, replace = TRUE),
    fw = sample(c("uniform", "optimize"), n_specs, replace = TRUE)
  )
  
  spec_curve_obj <- list(
    results = results
  )
  
  # Add inference results if requested
  if (include_inference) {
    # Abadie inference
    abadie_inference <- list(
      p_values_rmse_ratio = data.table(
        full_spec_id = paste0("spec_", 1:n_specs),
        unit_name = rep("treated_unit", n_specs),
        unit_type = rep("treated", n_specs),
        p_value = runif(n_specs, 0, 1)
      )
    )
    
    # Bootstrap inference
    bootstrap_inference <- list(
      p_values = data.table(
        full_spec_id = paste0("spec_", 1:n_specs),
        p_value_two_tailed = runif(n_specs, 0, 1)
      )
    )
    
    spec_curve_obj$abadie_inference <- abadie_inference
    spec_curve_obj$bootstrap_inference <- bootstrap_inference
  }
  
  class(spec_curve_obj) <- c("spec_curve", "list")
  return(spec_curve_obj)
}

test_that("tidy.scest works correctly", {
  mock_scest <- create_mock_scest_object()
  
  # Basic tidy call
  tidy_result <- tidy(mock_scest)
  
  # Check structure
  expect_s3_class(tidy_result, "data.table")
  expect_true(all(c("term", "estimate") %in% names(tidy_result)))
  expect_equal(nrow(tidy_result), 4)  # Four donors
  
  # Check derived columns
  expect_true("donor_rank" %in% names(tidy_result))
  expect_true("weight_percentage" %in% names(tidy_result))
  expect_true("is_active" %in% names(tidy_result))
  expect_true("constraint_type" %in% names(tidy_result))
  
  # Check values make sense
  expect_equal(sum(tidy_result$estimate), 1)  # Weights sum to 1 for simplex
  expect_equal(sum(tidy_result$weight_percentage), 100)  # Percentages sum to 100
  expect_true(all(tidy_result$is_active))  # All weights should be active
  expect_equal(unique(tidy_result$constraint_type), "simplex")
  
  # Check ordering (highest weight first)
  expect_true(all(diff(tidy_result$estimate) <= 0))
})

test_that("tidy.scest with confidence intervals works", {
  mock_scest <- create_mock_scest_object(include_bootstrap = TRUE)
  
  # Tidy with confidence intervals
  tidy_result <- tidy(mock_scest, conf.int = TRUE)
  
  # Check additional columns exist
  expect_true("std.error" %in% names(tidy_result))
  expect_true("conf.low" %in% names(tidy_result))
  expect_true("conf.high" %in% names(tidy_result))
  
  # Check confidence intervals make sense
  expect_true(all(tidy_result$conf.low <= tidy_result$estimate))
  expect_true(all(tidy_result$estimate <= tidy_result$conf.high))
  expect_true(all(tidy_result$std.error >= 0))
})

test_that("tidy.scest handles missing bootstrap data gracefully", {
  mock_scest <- create_mock_scest_object(include_bootstrap = FALSE)
  
  # Should warn but not error when requesting confidence intervals
  expect_warning(
    tidy_result <- tidy(mock_scest, conf.int = TRUE),
    "Confidence intervals requested but no bootstrap results available"
  )
  
  # Should still return the data with NA columns
  expect_true("std.error" %in% names(tidy_result))
  expect_true(all(is.na(tidy_result$std.error)))
})

test_that("glance.scest works correctly", {
  mock_scest <- create_mock_scest_object()
  
  # Basic glance call
  glance_result <- glance(mock_scest)
  
  # Check structure
  expect_s3_class(glance_result, "data.table")
  expect_equal(nrow(glance_result), 1)  # Single row summary
  
  # Check essential columns
  essential_cols <- c("n_donors", "n_pre_periods", "n_post_periods", 
                      "rmse_pre", "rmse_post", "constraint_type",
                      "n_active_donors", "max_weight", "herfindahl_index")
  expect_true(all(essential_cols %in% names(glance_result)))
  
  # Check values make sense
  expect_equal(glance_result$n_donors, 4)
  expect_equal(glance_result$n_pre_periods, 8)  
  expect_equal(glance_result$n_post_periods, 5)
  expect_equal(glance_result$constraint_type, "simplex")
  expect_equal(glance_result$max_weight, 0.4)  # Highest weight from mock data
  expect_true(glance_result$herfindahl_index > 0 && glance_result$herfindahl_index <= 1)
})

test_that("tidy.spec_curve works correctly", {
  mock_spec_curve <- create_mock_spec_curve_object()
  
  # Basic tidy call
  tidy_result <- tidy(mock_spec_curve)
  
  # Check structure
  expect_s3_class(tidy_result, "data.table")
  expect_true(all(c("spec_id", "estimate", "outcome", "unit_name") %in% names(tidy_result)))
  expect_equal(nrow(tidy_result), 20)  # Twenty specifications
  
  # Check specification characteristics are included
  expect_true("outcome_model" %in% names(tidy_result))
  expect_true("constraint_type" %in% names(tidy_result))
  expect_true("feature_weights" %in% names(tidy_result))
  
  # Check ordering (by effect size)
  expect_true(all(diff(tidy_result$estimate) >= 0))
})

test_that("tidy.spec_curve with inference works", {
  mock_spec_curve <- create_mock_spec_curve_object(include_inference = TRUE)
  
  # Tidy with inference
  tidy_result <- tidy(mock_spec_curve, conf.int = TRUE)
  
  # Check inference columns exist
  expect_true("p.value" %in% names(tidy_result))
  expect_true("method" %in% names(tidy_result))
  
  # Should also have bootstrap p-values
  expect_true("p.value.bootstrap" %in% names(tidy_result))
  
  # Check significance indicator
  if ("significant" %in% names(tidy_result)) {
    expect_true(all(is.logical(tidy_result$significant)))
  }
})

test_that("augment.scest works correctly", {
  mock_scest <- create_mock_scest_object()
  
  # Basic augment call
  augmented_result <- augment(mock_scest)
  
  # Check structure
  expect_s3_class(augmented_result, "data.table")
  expect_true(all(c(".fitted", ".period_type", ".is_treated") %in% names(augmented_result)))
  
  # Check period types
  expect_true(all(unique(augmented_result$.period_type) %in% c("pre", "post")))
  expect_true(all(augmented_result$.is_treated))  # All treated unit data
  
  # Check residuals if actual values present
  if (".resid" %in% names(augmented_result)) {
    expect_true(all(is.finite(augmented_result$.resid)))
  }
})

test_that("tidy methods handle invalid input gracefully", {
  # Test with wrong class
  fake_object <- list(some_data = 1:10)
  class(fake_object) <- "fake_class"
  
  # These should error with method dispatch or our validation errors
  expect_error(tidy(fake_object), "no applicable method|Object must be of class")
  expect_error(glance(fake_object), "no applicable method|Object must be of class") 
  expect_error(augment(fake_object), "no applicable method|Object must be of class")
})

test_that("tidy methods work with broom ecosystem", {
  # Test that our methods work with broom's approach
  mock_scest <- create_mock_scest_object()
  
  # Should be able to call the methods through generic dispatch
  expect_no_error(tidy_result <- tidy(mock_scest))
  expect_no_error(glance_result <- glance(mock_scest))
  expect_no_error(augmented_result <- augment(mock_scest))
  
  # Check that results are data.table (our enhancement over standard broom)
  expect_s3_class(tidy_result, "data.table")
  expect_s3_class(glance_result, "data.table") 
  expect_s3_class(augmented_result, "data.table")
})

test_that("tidy methods preserve data.table efficiency", {
  mock_scest <- create_mock_scest_object()
  
  # Test that our results can be used in data.table operations
  tidy_result <- tidy(mock_scest)
  
  # Should work with data.table syntax
  expect_no_error(
    filtered <- tidy_result[estimate > 0.2]
  )
  
  expect_no_error(
    summarized <- tidy_result[, .(mean_weight = mean(estimate)), by = constraint_type]
  )
  
  # Results should maintain data.table class
  expect_s3_class(filtered, "data.table")
  expect_s3_class(summarized, "data.table")
})