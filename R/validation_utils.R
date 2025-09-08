#' Input Validation Utility Functions
#'
#' @description
#' This file contains utility functions for input validation to ensure
#' consistent error handling across the SCMs package. These functions
#' provide informative error messages that help users identify and fix
#' input problems quickly.
#'
#' @name validation-utils
#' @keywords internal

#' Validate Data Frame Input
#'
#' @param data Object to validate as data frame
#' @param arg_name Character. Name of argument (for error messages)
#' @param allow_empty Logical. Whether to allow empty data frames
#' @param required_cols Character vector. Required column names
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_dataframe <- function(data, arg_name = "data", allow_empty = FALSE, required_cols = NULL) {
  if (!is.data.frame(data)) {
    stop("'", arg_name, "' must be a data frame. Received: ", class(data)[1])
  }
  
  if (!allow_empty && nrow(data) == 0) {
    stop("'", arg_name, "' cannot be empty.")
  }
  
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      stop("Missing required columns in '", arg_name, "': ", 
           paste(missing_cols, collapse = ", "), 
           ". Available columns: ", paste(names(data), collapse = ", "))
    }
  }
}

#' Validate Character Input
#'
#' @param x Object to validate
#' @param arg_name Character. Name of argument
#' @param length Expected length (default = 1)
#' @param choices Character vector of valid choices (optional)
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_character <- function(x, arg_name, length = 1, choices = NULL) {
  if (!is.character(x) || length(x) != length) {
    if (length == 1) {
      stop("'", arg_name, "' must be a single character string.")
    } else {
      stop("'", arg_name, "' must be a character vector of length ", length, ".")
    }
  }
  
  if (!is.null(choices)) {
    invalid <- setdiff(x, choices)
    if (length(invalid) > 0) {
      stop("Invalid values in '", arg_name, "': ", paste(invalid, collapse = ", "),
           ". Valid options: ", paste(choices, collapse = ", "))
    }
  }
}

#' Validate Numeric Input
#'
#' @param x Object to validate
#' @param arg_name Character. Name of argument
#' @param length Expected length (default = 1)
#' @param min Minimum allowed value (optional)
#' @param max Maximum allowed value (optional)
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_numeric <- function(x, arg_name, length = 1, min = NULL, max = NULL) {
  if (!is.numeric(x) || length(x) != length) {
    if (length == 1) {
      stop("'", arg_name, "' must be a single numeric value.")
    } else {
      stop("'", arg_name, "' must be a numeric vector of length ", length, ".")
    }
  }
  
  if (!is.null(min) && any(x < min, na.rm = TRUE)) {
    stop("'", arg_name, "' must be >= ", min, ". Received: ", paste(x, collapse = ", "))
  }
  
  if (!is.null(max) && any(x > max, na.rm = TRUE)) {
    stop("'", arg_name, "' must be <= ", max, ". Received: ", paste(x, collapse = ", "))
  }
}

#' Validate Logical Input
#'
#' @param x Object to validate
#' @param arg_name Character. Name of argument
#' @param length Expected length (default = 1)
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_logical <- function(x, arg_name, length = 1) {
  if (!is.logical(x) || length(x) != length) {
    if (length == 1) {
      stop("'", arg_name, "' must be TRUE or FALSE.")
    } else {
      stop("'", arg_name, "' must be a logical vector of length ", length, ".")
    }
  }
}

#' Validate List Input
#'
#' @param x Object to validate
#' @param arg_name Character. Name of argument
#' @param required_names Character vector of required list element names
#' @param allow_null Logical. Whether NULL is acceptable
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_list <- function(x, arg_name, required_names = NULL, allow_null = TRUE) {
  if (is.null(x) && allow_null) {
    return(invisible(NULL))
  }
  
  if (!is.list(x)) {
    stop("'", arg_name, "' must be a list.")
  }
  
  if (!is.null(required_names)) {
    missing_names <- setdiff(required_names, names(x))
    if (length(missing_names) > 0) {
      stop("Missing required elements in '", arg_name, "': ", 
           paste(missing_names, collapse = ", "))
    }
  }
}

#' Validate Column Existence in Data
#'
#' @param data Data frame
#' @param col_name Column name to check
#' @param arg_name Name of argument (for error messages)
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_column_exists <- function(data, col_name, arg_name) {
  if (!col_name %in% names(data)) {
    stop("Column '", col_name, "' (specified in '", arg_name, "') not found in data. ",
         "Available columns: ", paste(names(data), collapse = ", "))
  }
}

#' Validate Column is Numeric
#'
#' @param data Data frame
#' @param col_name Column name to check
#' @param arg_name Name of argument (for error messages)
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_column_numeric <- function(data, col_name, arg_name) {
  if (!is.numeric(data[[col_name]])) {
    stop("Column '", col_name, "' (specified in '", arg_name, "') must be numeric. ",
         "Current type: ", class(data[[col_name]])[1])
  }
}

#' Validate Values Exist in Column
#'
#' @param data Data frame
#' @param col_name Column name
#' @param values Values to check for
#' @param arg_name Name of argument (for error messages)
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_values_in_column <- function(data, col_name, values, arg_name) {
  available_values <- unique(data[[col_name]])
  missing_values <- setdiff(values, available_values)
  if (length(missing_values) > 0) {
    stop("Values not found in column '", col_name, "': ", paste(missing_values, collapse = ", "),
         " (specified in '", arg_name, "'). Available values: ", 
         paste(sort(available_values), collapse = ", "))
  }
}

#' Validate Period Relationships for Synthetic Control
#'
#' @param min_period Numeric. Minimum period
#' @param treated_period Numeric. Treatment period
#' @param end_period Numeric. End period
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_period_relationships <- function(min_period, treated_period, end_period) {
  if (min_period >= treated_period) {
    stop("'min_period' (", min_period, ") must be less than 'treated_period' (", 
         treated_period, "). There must be pre-treatment periods for synthetic control estimation.")
  }
  
  if (treated_period >= end_period) {
    stop("'treated_period' (", treated_period, ") must be less than 'end_period' (", 
         end_period, "). There must be post-treatment periods for effect evaluation.")
  }
  
  # Check for reasonable number of pre-treatment periods
  pre_periods <- treated_period - min_period
  if (pre_periods < 2) {
    warning("Only ", pre_periods, " pre-treatment period(s) available. ",
            "Consider extending the pre-treatment window for more reliable estimation.",
            call. = FALSE)
  }
}

#' Validate S3 Object Class
#'
#' @param obj Object to validate
#' @param required_class Character. Required class name
#' @param arg_name Name of argument (for error messages)
#' @return NULL (throws error if validation fails)
#' @keywords internal
validate_s3_class <- function(obj, required_class, arg_name) {
  if (!methods::is(obj, required_class)) {
    stop("'", arg_name, "' must be an object of class '", required_class, "'. ",
         "Received object of class: ", paste(class(obj), collapse = ", "))
  }
}