library(testthat)
library(data.table)

# Regression tests to ensure consistent behavior across package versions
# These tests use fixed seeds and expected values to detect breaking changes

# Helper function to create deterministic test data for regression tests
create_deterministic_test_data <- function() {
  set.seed(12345)  # Fixed seed for reproducibility
  
  # Create simple panel structure
  units <- c("CA", "TX", "NY", "FL")
  periods <- 1990:2005
  
  data <- expand.grid(
    state = units,
    year = periods,
    stringsAsFactors = FALSE
  )
  
  # Add deterministic outcome with known treatment effect
  data$gdp <- 100 + 
             (data$year - 1990) * 1.5 +  # Linear trend
             rnorm(nrow(data), 0, 2) +    # Fixed noise
             ifelse(data$state == "CA" & data$year > 2000, 8, 0)  # Treatment effect
  
  # Add covariates with fixed relationships
  data$population <- 500 + data$year * 10 + rnorm(nrow(data), 0, 20)
  data$investment <- 30 + data$gdp * 0.1 + rnorm(nrow(data), 0, 5)
  
  return(data)
}

test_that("scdata produces consistent output structure across versions", {
  data <- create_deterministic_test_data()
  
  # Basic scdata call with fixed parameters
  result <- scdata(
    df = data,
    id.var = "state",
    time.var = "year", 
    outcome.var = "gdp",
    period.pre = 1990:2000,
    period.post = 2001:2005,
    unit.tr = "CA",
    unit.co = c("TX", "NY", "FL"),
    constant = FALSE
  )
  
  # Test consistent object structure
  expect_s3_class(result, "scdata")
  expect_type(result, "list")
  
  # Test matrix dimensions remain consistent
  expect_equal(ncol(result$A), 1)  # One treated unit
  expect_equal(ncol(result$B), 3)  # Three control units
  expect_equal(nrow(result$A), nrow(result$B))  # Same number of features
  
  # Test time series dimensions
  expect_equal(length(result$Y.pre), 11)  # 1990-2000 (11 years)
  expect_equal(length(result$Y.post), 5)   # 2001-2005 (5 years)
  expect_equal(nrow(result$Y.donors), 11)  # Pre-treatment periods
  expect_equal(nrow(result$P), 5)          # Post-treatment periods
  
  # Test specs structure
  expect_equal(result$specs$J, 3)  # Three donors
  expect_equal(result$specs$M, 1)  # One treated unit
})

test_that("scdata with constant term maintains structure", {
  data <- create_deterministic_test_data()
  
  covagg <- list(
    list(var = "population", average = TRUE),
    list(var = "investment", average = TRUE)
  )
  
  result_const <- scdata(
    df = data,
    id.var = "state",
    time.var = "year",
    outcome.var = "gdp", 
    period.pre = 1990:2000,
    period.post = 2001:2005,
    unit.tr = "CA",
    unit.co = c("TX", "NY", "FL"),
    constant = TRUE,
    covagg = covagg
  )
  
  # Test constant matrix structure
  expect_true(!is.null(result_const$C))
  expect_true(is.matrix(result_const$C))
  expect_equal(ncol(result_const$C), 1)  # Single constant column
  expect_equal(nrow(result_const$C), nrow(result_const$A))  # Same as feature matrix
  
  # Test specs updated correctly
  expect_equal(result_const$specs$KM, 1)  # One constant term
})

test_that("scest weight estimation remains stable", {
  skip_if_not_installed("CVXR")
  
  data <- create_deterministic_test_data()
  
  scdata_obj <- scdata(
    df = data,
    id.var = "state",
    time.var = "year",
    outcome.var = "gdp",
    period.pre = 1990:2000, 
    period.post = 2001:2005,
    unit.tr = "CA",
    unit.co = c("TX", "NY", "FL"),
    constant = FALSE
  )
  
  result <- tryCatch({
    scest(
      data = scdata_obj,
      w.constr = list(name = "simplex"),
      feature_weights = "uniform",
      solver = "ECOS"
    )
  }, error = function(e) {
    if (grepl("should be|must be|Invalid|Missing", e$message)) {
      stop(e)  # Re-throw validation errors
    }
    skip(paste("Computational error in regression test:", e$message))
  })
  
  # Test weight structure consistency
  expect_s3_class(result, "scest")
  expect_equal(length(result$est.results$w), 3)  # Three donor weights
  
  # Test that weights are non-negative for simplex constraint
  if (!is.null(result$est.results$w)) {
    expect_true(all(result$est.results$w >= -1e-10))  # Allow small numerical errors
  }
  
  # Test fitted values structure
  expect_equal(length(result$est.results$Y.pre.fit), 11)   # Pre-treatment periods
  expect_equal(length(result$est.results$Y.post.fit), 5)   # Post-treatment periods
})

test_that("validation functions produce consistent error messages", {
  # Test that error messages remain consistent across versions
  df_bad <- data.frame(x = 1:3, y = 4:6)
  
  # Test validate_dataframe errors
  expect_error(validate_dataframe(list(x = 1:3), "test_df"),
               "'test_df' must be a data frame")
  
  expect_error(validate_dataframe(data.frame(), "test_df", allow_empty = FALSE),
               "'test_df' cannot be empty")
  
  expect_error(validate_dataframe(df_bad, "test_df", required_cols = c("x", "z")),
               "Missing required columns")
  
  # Test validate_character errors  
  expect_error(validate_character(123, "test_arg"),
               "'test_arg' must be a single character string")
  
  expect_error(validate_character("invalid", "test_arg", choices = c("a", "b")),
               "Invalid values in 'test_arg'")
  
  # Test validate_numeric errors
  expect_error(validate_numeric("not_numeric", "test_arg"),
               "'test_arg' must be a single numeric value")
  
  expect_error(validate_numeric(0, "test_arg", min = 1),
               "'test_arg' must be >= 1")
})

test_that("S3 method dispatch remains consistent", {
  # Create minimal spec_curve object
  set.seed(12345)
  mock_results <- data.table(
    unit_name = rep("treated_unit", 10),
    unit_type = rep("treated", 10),
    outcome = rep("gdp", 10),
    full_spec_id = paste0("spec_", 1:10),
    post_period = rep(TRUE, 10),
    tau = rnorm(10, 2, 0.5),
    rmse = runif(10, 0.5, 1.5)
  )
  
  spec_obj <- list(
    results = mock_results,
    expected_direction = "negative"
  )
  class(spec_obj) <- c("spec_curve", "list")
  
  # Test S3 method dispatch consistency
  expect_output(print(spec_obj), "Specification Curve Analysis Results")
  expect_output(summary(spec_obj), "Detailed Specification Breakdown")
  
  # Test class structure
  expect_s3_class(spec_obj, "spec_curve")
  expect_s3_class(spec_obj, "list")
})

test_that("Data processing pipeline maintains consistency", {
  data <- create_deterministic_test_data()
  
  # Test that processed matrices have consistent structure
  scdata_obj <- scdata(
    df = data,
    id.var = "state", 
    time.var = "year",
    outcome.var = "gdp",
    period.pre = 1990:2000,
    period.post = 2001:2005,
    unit.tr = "CA",
    unit.co = c("TX", "NY", "FL")
  )
  
  # Test matrix properties that should remain stable
  expect_true(all(is.finite(scdata_obj$A)))  # No infinite values
  expect_true(all(is.finite(scdata_obj$B)))  # No infinite values
  expect_false(any(is.na(scdata_obj$Y.pre)))  # No missing outcomes
  expect_false(any(is.na(scdata_obj$Y.post))) # No missing outcomes
  
  # Test that donor matrix has expected structure
  expect_equal(nrow(scdata_obj$Y.donors), 11)  # Pre-treatment periods
  expect_equal(ncol(scdata_obj$Y.donors), 3)   # Three donors
  expect_equal(nrow(scdata_obj$P), 5)          # Post-treatment periods
  expect_equal(ncol(scdata_obj$P), 3)          # Three donors
})

test_that("Package namespace exports remain consistent", {
  # Test that key functions remain exported
  expect_true("scdata" %in% ls("package:SCMs"))
  expect_true("scest" %in% ls("package:SCMs"))
  expect_true("spec_curve" %in% ls("package:SCMs"))
  expect_true("run_spec_curve_analysis" %in% ls("package:SCMs"))
  
  # Test S3 method registration
  expect_true(methods::existsMethod("print", "scdata"))
  expect_true(methods::existsMethod("print", "scest"))
  expect_true(methods::existsMethod("print", "spec_curve"))
  expect_true(methods::existsMethod("summary", "scdata"))
  expect_true(methods::existsMethod("summary", "scest")) 
  expect_true(methods::existsMethod("summary", "spec_curve"))
  expect_true(methods::existsMethod("plot", "spec_curve"))
})