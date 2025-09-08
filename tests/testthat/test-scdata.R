library(testthat)

# Helper function to create test data
create_test_data <- function() {
  set.seed(123)
  
  # Create panel data
  units <- c("CA", "TX", "NY", "FL", "WA")
  periods <- 1990:2005
  
  data <- expand.grid(
    state = units,
    year = periods,
    stringsAsFactors = FALSE
  )
  
  # Add outcome with treatment effect for CA after 2000
  data$gdp <- 100 + (data$year - 1990) * 2 + rnorm(nrow(data), 0, 5) +
             ifelse(data$state == "CA" & data$year > 2000, 10, 0)
  
  # Add covariates
  data$population <- 1000 + rnorm(nrow(data), 0, 100)
  data$investment <- 50 + rnorm(nrow(data), 0, 10)
  
  return(data)
}

test_that("scdata basic functionality works", {
  test_data <- create_test_data()
  
  # Test basic scdata creation
  expect_silent({
    result <- scdata(
      df = test_data,
      id.var = "state", 
      time.var = "year",
      outcome.var = "gdp",
      period.pre = 1990:2000,
      period.post = 2001:2005,
      unit.tr = "CA",
      unit.co = c("TX", "NY", "FL", "WA")
    )
  })
  
  # Check that result has scdata class
  result <- scdata(
    df = test_data,
    id.var = "state", 
    time.var = "year", 
    outcome.var = "gdp",
    period.pre = 1990:2000,
    period.post = 2001:2005,
    unit.tr = "CA",
    unit.co = c("TX", "NY", "FL", "WA")
  )
  
  expect_s3_class(result, "scdata")
  expect_type(result, "list")
  
  # Check essential components exist
  expect_true("A" %in% names(result))  # Pre-treatment features matrix
  expect_true("B" %in% names(result))  # Control units features
  expect_true("Y.pre" %in% names(result))  # Pre-treatment outcome for treated
  expect_true("Y.post" %in% names(result))  # Post-treatment outcome for treated
  expect_true("Y.donors" %in% names(result))  # Control units outcomes
  expect_true("specs" %in% names(result))  # Specifications
})

test_that("scdata input validation works", {
  test_data <- create_test_data()
  
  # Test non-dataframe input
  expect_error(
    scdata(df = list(x = 1), id.var = "state", time.var = "year", outcome.var = "gdp",
           period.pre = 1990:2000, period.post = 2001:2005, unit.tr = "CA", unit.co = "TX"),
    "Data input should be a dataframe object"
  )
  
  # Test missing id.var
  expect_error(
    scdata(df = test_data, id.var = 123, time.var = "year", outcome.var = "gdp", 
           period.pre = 1990:2000, period.post = 2001:2005, unit.tr = "CA", unit.co = "TX"),
    "You should specify the name of id.var as a character"
  )
  
  # Test missing outcome.var
  expect_error(
    scdata(df = test_data, id.var = "state", time.var = "year", outcome.var = 123,
           period.pre = 1990:2000, period.post = 2001:2005, unit.tr = "CA", unit.co = "TX"),
    "You should specify the name of outcome.var as a character"
  )
  
  # Test missing time.var  
  expect_error(
    scdata(df = test_data, id.var = "state", time.var = 123, outcome.var = "gdp",
           period.pre = 1990:2000, period.post = 2001:2005, unit.tr = "CA", unit.co = "TX"),
    "You should specify the name of time.var as a character"
  )
  
  # Test missing column
  expect_error(
    scdata(df = test_data, id.var = "missing_col", time.var = "year", outcome.var = "gdp",
           period.pre = 1990:2000, period.post = 2001:2005, unit.tr = "CA", unit.co = "TX"),
    "ID variable \\(id.var\\) not found in the input dataframe"
  )
  
  # Test non-numeric outcome
  bad_data <- test_data
  bad_data$gdp <- as.character(bad_data$gdp)
  expect_error(
    scdata(df = bad_data, id.var = "state", time.var = "year", outcome.var = "gdp",
           period.pre = 1990:2000, period.post = 2001:2005, unit.tr = "CA", unit.co = "TX"),
    "Outcome variable \\(outcome.var\\) must be numeric"
  )
})

test_that("scdata handles covariate specifications", {
  test_data <- create_test_data()
  
  # Test with basic covariate aggregation
  covagg_spec <- list(
    list(var = "population", average = TRUE),
    list(var = "investment", each = TRUE)
  )
  
  expect_silent({
    result <- scdata(
      df = test_data,
      id.var = "state",
      time.var = "year", 
      outcome.var = "gdp",
      period.pre = 1990:2000,
      period.post = 2001:2005,
      unit.tr = "CA",
      unit.co = c("TX", "NY", "FL", "WA"),
      covagg = covagg_spec
    )
  })
  
  # Check that matrices have appropriate dimensions
  result <- scdata(
    df = test_data,
    id.var = "state",
    time.var = "year",
    outcome.var = "gdp", 
    period.pre = 1990:2000,
    period.post = 2001:2005,
    unit.tr = "CA",
    unit.co = c("TX", "NY", "FL", "WA"),
    covagg = covagg_spec
  )
  
  # A matrix should have features as rows, control units as columns
  expect_true(is.matrix(result$A))
  expect_true(nrow(result$A) > 0)  # Should have features
  expect_equal(ncol(result$A), 1)  # One treated unit
  
  # B matrix should have same number of rows as A, columns = number of donors
  expect_true(is.matrix(result$B))
  expect_equal(nrow(result$B), nrow(result$A))
  expect_equal(ncol(result$B), 4)  # Four control units
})

test_that("scdata handles constant terms correctly", {
  test_data <- create_test_data()
  
  covagg_spec <- list(
    list(var = "population", average = TRUE),
    list(var = "investment", average = TRUE)
  )
  
  # Test without constant
  result_no_const <- scdata(
    df = test_data,
    id.var = "state",
    time.var = "year",
    outcome.var = "gdp",
    period.pre = 1990:2000, 
    period.post = 2001:2005,
    unit.tr = "CA",
    unit.co = c("TX", "NY", "FL", "WA"),
    covagg = covagg_spec,
    constant = FALSE
  )
  
  expect_null(result_no_const$C)
  
  # Test with constant  
  result_with_const <- scdata(
    df = test_data,
    id.var = "state",
    time.var = "year", 
    outcome.var = "gdp",
    period.pre = 1990:2000,
    period.post = 2001:2005,
    unit.tr = "CA",
    unit.co = c("TX", "NY", "FL", "WA"),
    covagg = covagg_spec,
    constant = TRUE
  )
  
  expect_true(!is.null(result_with_const$C))
  expect_true(is.matrix(result_with_const$C))
})

test_that("scdata print and summary methods work", {
  test_data <- create_test_data()
  
  result <- scdata(
    df = test_data,
    id.var = "state",
    time.var = "year",
    outcome.var = "gdp", 
    period.pre = 1990:2000,
    period.post = 2001:2005,
    unit.tr = "CA",
    unit.co = c("TX", "NY", "FL", "WA")
  )
  
  # Test print method
  expect_output(print(result), "Synthetic Control Data")
  
  # Test summary method
  expect_output(summary(result), "Synthetic Control Data")
})

test_that("scdata handles edge cases", {
  test_data <- create_test_data()
  
  # Test with minimal time periods
  expect_warning({
    result <- scdata(
      df = test_data,
      id.var = "state", 
      time.var = "year",
      outcome.var = "gdp",
      period.pre = c(1990, 1991),  # Very few pre-periods
      period.post = c(2001, 2002),
      unit.tr = "CA", 
      unit.co = "TX"
    )
  })
  
  # Test with single control unit
  expect_silent({
    result <- scdata(
      df = test_data,
      id.var = "state",
      time.var = "year",
      outcome.var = "gdp",
      period.pre = 1990:2000,
      period.post = 2001:2005, 
      unit.tr = "CA",
      unit.co = "TX"  # Only one control
    )
  })
  
  expect_equal(ncol(result$B), 1)
})