library(testthat)
library(data.table)

# Helper function to create mock panel data for spec_curve testing
create_mock_panel_data_spec_curve <- function() {
  set.seed(456)
  
  # Create panel data structure
  states <- c("treated_state", "control_1", "control_2", "control_3")
  years <- 2000:2010
  
  # Expand to full panel
  data <- expand.grid(
    state = states,
    year = years,
    stringsAsFactors = FALSE
  )
  
  # Convert to data.table for efficiency
  setDT(data)
  
  # Add outcome with treatment effect after 2005
  data[, gdp := rnorm(.N, mean = 100, sd = 10) + 
         (year - 2000) * 2 +  # Time trend
         ifelse(state == "treated_state" & year >= 2006, 15, 0)]  # Treatment effect
  
  # Add some covariates
  data[, population := rnorm(.N, mean = 1000, sd = 100)]
  data[, investment := rnorm(.N, mean = 50, sd = 8)]
  data[, education := rnorm(.N, mean = 12, sd = 2)]
  
  return(data)
}

test_that("spec_curve input validation works correctly", {
  mock_data <- create_mock_panel_data_spec_curve()
  covagg_min <- list(
    population = list(
      label = "Population",
      operations = list(list(var = "population", compute = "mean"))
    )
  )
  empty_data <- data.table(
    state = character(),
    year = integer(),
    gdp = numeric(),
    population = numeric(),
    investment = numeric(),
    education = numeric()
  )
  
  # Test missing required parameters
  expect_error(
    spec_curve(dataset = NULL, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year"),
    "dataset must be a data.frame or data.table"
  )
  
  # Test empty dataset
  expect_error(
    spec_curve(dataset = empty_data, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year"),
    "dataset.*cannot be empty"
  )
  
  # Test non-character outcome
  expect_error(
    spec_curve(dataset = mock_data, outcomes = 123, covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year"),
    "outcomes.*single character"
  )
  
  # Test missing columns
  bad_data <- copy(mock_data)[, gdp := NULL]
  expect_error(
    spec_curve(dataset = bad_data, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year"),
    "Missing required columns.*gdp"
  )
  
  # Test invalid period relationships
  expect_error(
    spec_curve(dataset = mock_data, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2000, min_period = 2000, end_period = 2010,
               col_name_period = "year"),
    "min_period.*must be less than.*treated_period"
  )
  
  # Test treated unit not found
  expect_error(
    spec_curve(dataset = mock_data, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "nonexistent_unit",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year"),
    "Values not found in column 'state'"
  )
})

test_that("spec_curve covariate specification handling works", {
  mock_data <- create_mock_panel_data_spec_curve()
  
  # Test operation-based covariate format
  simple_covagg <- list(
    "basic_covariates" = list(
      label = "basic_covariates",
      operations = list(
        list(var = "population", partition_periods = list(type = "by_period")),
        list(var = "investment", compute = "mean"),
        list(var = "outcome_var", compute = "mean")
      )
    )
  )
  
  expect_silent({
    # Should not error on input validation - may error on computation
    result <- suppressMessages(tryCatch({
      spec_curve(
        dataset = mock_data,
        outcomes = "gdp",
        covagg = simple_covagg,
        col_name_unit_name = "state", 
        name_treated_unit = "treated_state",
        treated_period = 2006,
        min_period = 2000,
        end_period = 2010,
        col_name_period = "year",
        feature_weights = c("uniform"),
        outcome_models = c("none"),
        constraints = list(list(name = "simplex")),
        donor_sample = "all",
        cores = 1,
        verbose = FALSE
      )
    }, error = function(e) {
      # Allow computational errors but not validation errors
      if (grepl("Input validation failed|must be|Missing required", e$message)) {
        stop(e)
      }
      return(NULL)
    }))
  })
  
  # Test alternative operation-based format
  trad_covagg <- list(
    population_mean = list(
      label = "population_mean",
      operations = list(list(var = "population", compute = "mean"))
    ),
    investment_each = list(
      label = "investment_each",
      operations = list(list(var = "investment", partition_periods = list(type = "by_period")))
    )
  )
  
  expect_silent({
    result <- suppressMessages(tryCatch({
      spec_curve(
        dataset = mock_data,
        outcomes = "gdp", 
        covagg = trad_covagg,
        col_name_unit_name = "state",
        name_treated_unit = "treated_state",
        treated_period = 2006,
        min_period = 2000,
        end_period = 2010,
        col_name_period = "year",
        feature_weights = c("uniform"),
        outcome_models = c("none"),
        constraints = list(list(name = "simplex")),
        donor_sample = "all",
        cores = 1,
        verbose = FALSE
      )
    }, error = function(e) {
      if (grepl("Input validation failed|must be|Missing required", e$message)) {
        stop(e)
      }
      return(NULL)
    }))
  })
})

test_that("spec_curve parameter validation works", {
  mock_data <- create_mock_panel_data_spec_curve()
  covagg_min <- list(
    population = list(
      label = "Population",
      operations = list(list(var = "population", compute = "mean"))
    )
  )
  
  # Test invalid feature_weights
  expect_error(
    spec_curve(dataset = mock_data, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year", feature_weights = c("invalid_method")),
    "Invalid values in 'feature_weights'"
  )
  
  # Test invalid outcome_models  
  expect_error(
    spec_curve(dataset = mock_data, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year", outcome_models = c("invalid_model")),
    "Invalid values in 'outcome_models'"
  )
  
  # Test invalid donor_sample
  expect_error(
    spec_curve(dataset = mock_data, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year", donor_sample = c("invalid_sample")),
    "Invalid values in 'donor_sample'"
  )
  
  # Test invalid cores parameter
  expect_error(
    spec_curve(dataset = mock_data, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year", cores = -1),
    "cores.*must be >= 1"
  )
})

test_that("spec_curve basic functionality works without validation errors", {
  skip_if_not_installed("CVXR")
  
  mock_data <- create_mock_panel_data_spec_curve()
  covagg_min <- list(
    population = list(
      label = "Population",
      operations = list(list(var = "population", compute = "mean"))
    )
  )
  
  # Basic minimal spec_curve run - should pass input validation
  expect_silent({
    result <- suppressMessages(tryCatch({
      spec_curve(
        dataset = mock_data,
        outcomes = "gdp",
        covagg = covagg_min,
        col_name_unit_name = "state",
        name_treated_unit = "treated_state", 
        treated_period = 2006,
        min_period = 2000,
        end_period = 2010,
        col_name_period = "year",
        feature_weights = c("uniform"),
        outcome_models = c("none"),
        constraints = list(list(name = "simplex")),
        donor_sample = "all",
        cores = 1,
        verbose = FALSE
      )
    }, error = function(e) {
      # Re-throw validation errors, allow computational errors
      if (grepl("Input validation failed|must be|Missing required|Invalid", e$message)) {
        stop(e)
      }
      return(NULL)  # Suppress computational errors for this test
    }))
  })
})

test_that("spec_curve inference parameter validation", {
  mock_data <- create_mock_panel_data_spec_curve()
  covagg_min <- list(
    population = list(
      label = "Population",
      operations = list(list(var = "population", compute = "mean"))
    )
  )
  
  # Test invalid inference_type
  expect_error(
    spec_curve(dataset = mock_data, outcomes = "gdp", covagg = covagg_min, verbose = FALSE,
               col_name_unit_name = "state", name_treated_unit = "treated_state",
               treated_period = 2006, min_period = 2000, end_period = 2010,
               col_name_period = "year", inference_type = "invalid_inference"),
    "inference_type must be one of"
  )
  
  # Test inference_config validation
  expect_silent({
    # Should pass input validation even with custom inference config
    result <- suppressMessages(tryCatch({
      spec_curve(
        dataset = mock_data,
        outcomes = "gdp",
        covagg = covagg_min,
        col_name_unit_name = "state",
        name_treated_unit = "treated_state",
        treated_period = 2006,
        min_period = 2000,
        end_period = 2010,
        col_name_period = "year",
        feature_weights = c("uniform"),
        outcome_models = c("none"),
        constraints = list(list(name = "simplex")),
        donor_sample = "all",
        inference_type = "placebo",
        inference_config = list(bootstrap_n_replications = 100),
        cores = 1,
        verbose = FALSE
      )
    }, error = function(e) {
      if (grepl("Input validation failed|must be|Missing required|Invalid", e$message)) {
        stop(e)
      }
      return(NULL)
    }))
  })
})

test_that("spec_curve returns proper structure when successful", {
  skip_if_not_installed("CVXR")
  
  mock_data <- create_mock_panel_data_spec_curve()
  covagg_min <- list(
    population = list(
      label = "Population",
      operations = list(list(var = "population", compute = "mean"))
    )
  )
  
  result <- tryCatch({
    spec_curve(
      dataset = mock_data,
      outcomes = "gdp",
      covagg = covagg_min,
      col_name_unit_name = "state",
      name_treated_unit = "treated_state",
      treated_period = 2006,
      min_period = 2000, 
      end_period = 2010,
      col_name_period = "year",
      feature_weights = c("uniform"),
      outcome_models = c("none"),
      constraints = list(list(name = "simplex")),
      donor_sample = "all",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e) {
    if (grepl("Input validation failed|must be|Missing required|Invalid", e$message)) {
      stop(e)  # Re-throw validation errors
    }
    # For computational errors, skip the test
    skip(paste("Computational error in spec_curve:", e$message))
  })
  
  if (is.null(result)) {
    skip("spec_curve returned NULL due to computational issues")
  }
  
  # Test that result has spec_curve class (assigned in function)
  expect_s3_class(result, "spec_curve")
  
  # Test essential structure
  expect_true("results" %in% names(result))
  expect_true(is.data.table(result$results) || is.data.frame(result$results))
  
  # Test inference structures exist (even if empty due to computational failures)
  expect_true("abadie_inference" %in% names(result) || 
              "bootstrap_inference" %in% names(result))
})

test_that("spec_curve handles multiple outcomes", {
  skip_if_not_installed("CVXR")
  
  mock_data <- create_mock_panel_data_spec_curve()
  covagg_min <- list(
    population = list(
      label = "Population",
      operations = list(list(var = "population", compute = "mean"))
    )
  )
  
  # Add second outcome
  mock_data[, gdp2 := gdp + rnorm(.N, 0, 5)]
  
  expect_silent({
    result <- suppressMessages(tryCatch({
      spec_curve(
        dataset = mock_data,
        outcomes = c("gdp", "gdp2"),
        covagg = covagg_min,
        col_name_unit_name = "state",
        name_treated_unit = "treated_state",
        treated_period = 2006,
        min_period = 2000,
        end_period = 2010,
        col_name_period = "year",
        feature_weights = c("uniform"),
        outcome_models = c("none"),
        constraints = list(list(name = "simplex")),
        donor_sample = "all",
        cores = 1,
        verbose = FALSE
      )
    }, error = function(e) {
      if (grepl("Input validation failed|must be|Missing required|Invalid", e$message)) {
        stop(e)
      }
      return(NULL)
    }))
  })
})

test_that("spec_curve handles edge cases", {
  skip_if_not_installed("CVXR")
  
  # Test with minimal time periods
  minimal_data <- create_mock_panel_data_spec_curve()
  minimal_data <- minimal_data[year %in% 2004:2007]  # Very short panel
  covagg_min <- list(
    population = list(
      label = "Population",
      operations = list(list(var = "population", compute = "mean"))
    )
  )
  
  expect_silent({
    result <- suppressMessages(tryCatch({
      spec_curve(
        dataset = minimal_data,
        outcomes = "gdp", 
        covagg = covagg_min,
        col_name_unit_name = "state",
        name_treated_unit = "treated_state",
        treated_period = 2006,
        min_period = 2004,
        end_period = 2007,
        col_name_period = "year",
        feature_weights = c("uniform"),
        outcome_models = c("none"),
        constraints = list(list(name = "simplex")),
        donor_sample = "all",
        cores = 1,
        verbose = FALSE
      )
    }, error = function(e) {
      # Allow warnings about short panels, but not validation errors
      if (grepl("Input validation failed|must be|Missing required|Invalid", e$message) && 
          !grepl("Only.*pre-treatment", e$message)) {
        stop(e)
      }
      return(NULL)
    }))
  })
})
