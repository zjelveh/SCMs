library(testthat)

test_that("validate_dataframe works correctly", {
  # Valid data frame
  df <- data.frame(x = 1:3, y = 4:6)
  expect_silent(validate_dataframe(df, "test_df"))
  
  # Test with required columns
  expect_silent(validate_dataframe(df, "test_df", required_cols = c("x", "y")))
  
  # Error cases
  expect_error(validate_dataframe(list(x = 1:3), "test_df"), 
               "'test_df' must be a data frame")
  expect_error(validate_dataframe(data.frame(), "test_df", allow_empty = FALSE),
               "'test_df' cannot be empty")
  expect_error(validate_dataframe(df, "test_df", required_cols = c("x", "z")),
               "Missing required columns")
})

test_that("validate_character works correctly", {
  # Valid single character
  expect_silent(validate_character("test", "arg"))
  
  # Valid character vector
  expect_silent(validate_character(c("a", "b"), "arg", length = 2))
  
  # Valid choices
  expect_silent(validate_character("a", "arg", choices = c("a", "b", "c")))
  
  # Error cases
  expect_error(validate_character(123, "arg"),
               "'arg' must be a single character string")
  expect_error(validate_character(c("a", "b"), "arg", length = 1),
               "'arg' must be a single character string")
  expect_error(validate_character("d", "arg", choices = c("a", "b", "c")),
               "Invalid values in 'arg'")
})

test_that("validate_numeric works correctly", {
  # Valid numbers
  expect_silent(validate_numeric(5, "arg"))
  expect_silent(validate_numeric(c(1, 2, 3), "arg", length = 3))
  
  # With bounds
  expect_silent(validate_numeric(5, "arg", min = 1, max = 10))
  
  # Error cases
  expect_error(validate_numeric("5", "arg"),
               "'arg' must be a single numeric value")
  expect_error(validate_numeric(c(1, 2), "arg", length = 1),
               "'arg' must be a single numeric value")
  expect_error(validate_numeric(0, "arg", min = 1),
               "'arg' must be >= 1")
  expect_error(validate_numeric(11, "arg", max = 10),
               "'arg' must be <= 10")
})

test_that("validate_logical works correctly", {
  # Valid logical
  expect_silent(validate_logical(TRUE, "arg"))
  expect_silent(validate_logical(c(TRUE, FALSE), "arg", length = 2))
  
  # Error cases
  expect_error(validate_logical("TRUE", "arg"),
               "'arg' must be TRUE or FALSE")
  expect_error(validate_logical(c(TRUE, FALSE), "arg", length = 1),
               "'arg' must be TRUE or FALSE")
})

test_that("validate_period_relationships works correctly", {
  # Valid periods
  expect_silent(validate_period_relationships(1990, 2000, 2010))
  
  # Should give warning for few pre-periods
  expect_warning(validate_period_relationships(1999, 2000, 2010),
                 "Only 1 pre-treatment period")
  
  # Error cases
  expect_error(validate_period_relationships(2000, 2000, 2010),
               "must be less than 'treated_period'")
  expect_error(validate_period_relationships(1990, 2010, 2010),
               "must be less than 'end_period'")
})

test_that("validate_values_in_column works correctly", {
  df <- data.frame(state = c("CA", "TX", "NY"), year = c(2000, 2001, 2002))
  
  # Valid values
  expect_silent(validate_values_in_column(df, "state", "CA", "treated_unit"))
  
  # Error case
  expect_error(validate_values_in_column(df, "state", "FL", "treated_unit"),
               "Values not found in column")
})
