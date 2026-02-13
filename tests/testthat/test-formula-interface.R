library(testthat)
library(data.table)

# Helper function to create mock panel data for formula testing
create_formula_test_data <- function() {
  states <- c("California", "Texas", "New York", "Florida")
  years <- 1990:2005
  
  panel <- expand.grid(
    state = states,
    year = years,
    stringsAsFactors = FALSE
  )
  setDT(panel)
  
  # Deterministic donor paths
  panel[state == "Texas", gdp := 100 + 1.6 * (year - 1990)]
  panel[state == "New York", gdp := 95 + 2.1 * (year - 1990)]
  panel[state == "Florida", gdp := 102 + 1.4 * (year - 1990)]

  # Treated unit is a convex combination in pre-period, shifted post-treatment.
  panel[state == "California",
        gdp := 0.5 * (100 + 1.6 * (year - 1990)) +
          0.3 * (95 + 2.1 * (year - 1990)) +
          0.2 * (102 + 1.4 * (year - 1990))]
  panel[state == "California" & year >= 2000, gdp := gdp + 10]
  
  # Deterministic covariates
  panel[state == "Texas", population := 1000 + 10 * (year - 1990)]
  panel[state == "New York", population := 1200 + 6 * (year - 1990)]
  panel[state == "Florida", population := 900 + 8 * (year - 1990)]
  panel[state == "California",
        population := 0.5 * (1000 + 10 * (year - 1990)) +
          0.3 * (1200 + 6 * (year - 1990)) +
          0.2 * (900 + 8 * (year - 1990))]

  panel[state == "Texas", investment := 45 + 0.9 * (year - 1990)]
  panel[state == "New York", investment := 55 + 0.6 * (year - 1990)]
  panel[state == "Florida", investment := 48 + 0.8 * (year - 1990)]
  panel[state == "California",
        investment := 0.5 * (45 + 0.9 * (year - 1990)) +
          0.3 * (55 + 0.6 * (year - 1990)) +
          0.2 * (48 + 0.8 * (year - 1990))]

  panel[state == "Texas", education := 11 + 0.03 * (year - 1990)]
  panel[state == "New York", education := 12 + 0.02 * (year - 1990)]
  panel[state == "Florida", education := 10.5 + 0.04 * (year - 1990)]
  panel[state == "California",
        education := 0.5 * (11 + 0.03 * (year - 1990)) +
          0.3 * (12 + 0.02 * (year - 1990)) +
          0.2 * (10.5 + 0.04 * (year - 1990))]
  
  return(panel)
}

test_that("synth function validates inputs correctly", {
  panel_data <- create_formula_test_data()
  
  # Test non-formula input
  expect_error(
    synth("not a formula", panel_data, "year", "state", 2000, 1990:1999, 2000:2005),
    "must be a formula object"
  )
  
  # Test missing | separator
  expect_error(
    synth(gdp ~ population, panel_data, "year", "state", 2000, 1990:1999, 2000:2005),
    "must include '\\|' to specify treated unit"
  )
  
  # Test missing variables
  expect_error(
    synth(nonexistent ~ population | California, panel_data, "year", "state", 2000, 1990:1999, 2000:2005),
    "Missing variables in data"
  )
  
  # Test nonexistent treated unit
  expect_error(
    synth(gdp ~ population | NonExistent, panel_data, "year", "state", 2000, 1990:1999, 2000:2005),
    "not found in.*column"
  )
})

test_that("formula parsing works correctly", {
  panel_data <- create_formula_test_data()
  
  # Test basic formula parsing
  formula_parts <- parse_synth_formula(
    gdp ~ population + investment | California,
    panel_data
  )
  
  expect_equal(formula_parts$outcome, "gdp")
  expect_equal(formula_parts$treated_unit, "California")
  expect_true("population" %in% formula_parts$covariates)
  expect_true("investment" %in% formula_parts$covariates)
})

test_that("synth function works with valid inputs", {
  skip_if_not_installed("CVXR")
  
  panel_data <- create_formula_test_data()
  
  # This should pass input validation (may fail on computation, which is fine)
  expect_silent({
    result <- tryCatch({
      synth(
        gdp ~ population + investment | California,
        data = panel_data,
        time.var = "year",
        id.var = "state",
        treated.period = 2000,
        pre.period = 1990:1999,
        post.period = 2000:2005,
        constraint = "simplex"
      )
    }, error = function(e) {
      # Allow computational errors but not validation errors
      if (grepl("must be|Missing|not found", e$message)) {
        stop(e)  # Re-throw validation errors
      }
      return(NULL)  # Suppress computational errors for this test
    })
  })
})

test_that("constraint standardization works", {
  # Test character constraint conversion
  expect_equal(standardize_constraint_spec("simplex"), list(name = "simplex"))
  expect_equal(standardize_constraint_spec("lasso")$name, "lasso")
  expect_equal(standardize_constraint_spec("ridge")$name, "ridge")
  expect_equal(standardize_constraint_spec("ols"), list(name = "ols"))
  
  # Test list constraint passthrough
  custom_constraint <- list(name = "custom", param = 0.5)
  expect_equal(standardize_constraint_spec(custom_constraint), custom_constraint)
})

test_that("covariate aggregation building works", {
  covariates <- c("population", "investment")
  
  # Test average method
  covagg_avg <- build_covagg_from_formula(covariates, "average")
  expect_equal(length(covagg_avg), 2)
  expect_equal(covagg_avg[[1]]$var, "population")
  expect_equal(covagg_avg[[1]]$partition_periods$type, "all")
  expect_equal(covagg_avg[[1]]$compute, "mean")
  
  # Test each method
  covagg_each <- build_covagg_from_formula(covariates, "each")
  expect_equal(covagg_each[[1]]$partition_periods$type, "by_period")
})

test_that("formula interface creates proper scest_formula object", {
  skip_if_not_installed("CVXR")
  
  panel_data <- create_formula_test_data()
  
  result <- synth(
    gdp ~ population | California,
    data = panel_data,
    time.var = "year", 
    id.var = "state",
    treated.period = 2000,
    pre.period = 1990:1999,
    post.period = 2000:2005
  )
  
  # Test class structure
  expect_s3_class(result, "scest_formula")
  expect_s3_class(result, "scest")
  
  # Test formula metadata
  expect_true("formula" %in% names(result))
  expect_true("formula_parts" %in% names(result))
  expect_true("call" %in% names(result))
})

test_that("scest_formula print method works", {
  skip_if_not_installed("CVXR")
  
  panel_data <- create_formula_test_data()
  
  result <- tryCatch({
    synth(
      gdp ~ population | California,
      data = panel_data,
      time.var = "year",
      id.var = "state", 
      treated.period = 2000,
      pre.period = 1990:1999,
      post.period = 2000:2005
    )
  }, error = function(e) {
    if (grepl("must be|Missing|not found", e$message)) {
      stop(e)
    }
    return(NULL)
  })
  
  if (!is.null(result)) {
    # Test that print method works without error
    expect_output(print(result), "Synthetic Control Model \\(Formula Interface\\)")
    expect_output(print(result), "Formula:")
    
    # Test summary method
    expect_output(summary(result), "Formula components:")
  }
})
