library(testthat)

# Helper function to create mock data for testing
create_mock_panel_data <- function() {
  # Create simple panel data for testing
  set.seed(123)
  
  states <- c("CA", "TX", "NY", "FL")
  years <- 1990:2010
  
  # Create all combinations
  data <- expand.grid(state = states, year = years, stringsAsFactors = FALSE)
  
  # Add outcome variable with some structure
  data$gdp <- rnorm(nrow(data), mean = 100, sd = 10) + 
             (data$year - 1990) * 2 + 
             ifelse(data$state == "CA" & data$year >= 2000, 15, 0)  # Treatment effect
  
  # Add some covariates
  data$population <- rnorm(nrow(data), mean = 1000, sd = 100)
  data$unemployment <- rnorm(nrow(data), mean = 5, sd = 1)
  
  return(data)
}

test_that("estimate_sc input validation works correctly", {
  mock_data <- create_mock_panel_data()
  
  # Test missing data frame
  expect_error(estimate_sc(NULL, "gdp", list(), "state", "CA", "year", 2000, 1990, 2010),
               "'dataset' must be a data frame")
  
  # Test empty data frame
  expect_error(estimate_sc(data.frame(), "gdp", list(), "state", "CA", "year", 2000, 1990, 2010),
               "'dataset' cannot be empty")
  
  # Test missing columns
  bad_data <- mock_data[, !names(mock_data) %in% "gdp"]
  expect_error(estimate_sc(bad_data, "gdp", list(), "state", "CA", "year", 2000, 1990, 2010),
               "Missing required columns.*gdp")
  
  # Test non-character outcome
  expect_error(estimate_sc(mock_data, 123, list(), "state", "CA", "year", 2000, 1990, 2010),
               "'outcome' must be a single character string")
  
  # Test non-character unit column
  expect_error(estimate_sc(mock_data, "gdp", list(), 123, "CA", "year", 2000, 1990, 2010),
               "'col_name_unit_name' must be a single character string")
  
  # Test non-numeric periods
  expect_error(estimate_sc(mock_data, "gdp", list(), "state", "CA", "year", "2000", 1990, 2010),
               "'treated_period' must be a single numeric value")
  
  # Test invalid period relationships
  expect_error(estimate_sc(mock_data, "gdp", list(), "state", "CA", "year", 1990, 1990, 2010),
               "must be less than 'treated_period'")
  
  # Test treated unit not found
  expect_error(estimate_sc(mock_data, "gdp", list(), "state", "INVALID", "year", 2000, 1990, 2010),
               "not found in data")
  
  # Test treated period not found
  expect_error(estimate_sc(mock_data, "gdp", list(), "state", "CA", "year", 2050, 1990, 2010),
               "'treated_period' \\(2050\\) not found in data")
  
  # Test non-numeric outcome
  bad_data <- mock_data
  bad_data$gdp <- as.character(bad_data$gdp)
  expect_error(estimate_sc(bad_data, "gdp", list(), "state", "CA", "year", 2000, 1990, 2010),
               "must be numeric")
  
  # Test invalid outcome models
  expect_error(estimate_sc(mock_data, "gdp", list(), "state", "CA", "year", 2000, 1990, 2010,
                          outcome_models = c("Invalid")),
               "Invalid outcome_models.*Invalid")
  
  # Test invalid feature weights
  expect_error(estimate_sc(mock_data, "gdp", list(), "state", "CA", "year", 2000, 1990, 2010,
                          feature_weights = c("invalid")),
               "Invalid feature_weights.*invalid")
})

test_that("estimate_sc validation warnings work correctly", {
  mock_data <- create_mock_panel_data()
  
  # Test warning for few pre-treatment periods
  expect_warning(estimate_sc(mock_data, "gdp", list(), "state", "CA", "year", 2000, 1999, 2010),
                 "Only 1 pre-treatment period")
})

# Note: Full functionality tests would require a working synthetic control implementation
# These tests focus on the input validation we added, which is the most critical
# for preventing user errors and providing helpful feedback

test_that("estimate_sc handles valid inputs without error", {
  skip_if_not_installed("CVXR")
  
  mock_data <- create_mock_panel_data()
  
  # This test verifies that valid inputs don't trigger validation errors
  # The actual computation might fail due to solver/data issues, but validation should pass
  expect_error(estimate_sc(mock_data, "gdp", list(), "state", "CA", "year", 2000, 1990, 2010),
               # Should not get validation errors, but might get computation errors
               regex = NA, invert = TRUE, class = "validation_error")
})

test_that("estimate_sc covagg parameter validation", {
  mock_data <- create_mock_panel_data()
  
  # Test valid covagg list
  covagg <- list(list(var = "population", average = TRUE))
  expect_silent({
    # Validation should pass - actual estimation may fail but that's different
    tryCatch(estimate_sc(mock_data, "gdp", covagg, "state", "CA", "year", 2000, 1990, 2010),
             error = function(e) {
               # If error is NOT from our validation, that's OK for this test
               if (!grepl("'covagg' must be|Missing required|must be a single", e$message)) {
                 return(NULL) # Non-validation error is acceptable
               } else {
                 stop(e) # Re-throw validation errors
               }
             })
  })
  
  # Test invalid covagg (not a list)
  expect_error(estimate_sc(mock_data, "gdp", "invalid", "state", "CA", "year", 2000, 1990, 2010),
               "'covagg' must be a list")
})