#' Process Covariate Aggregations for Synthetic Control (Version 2)
#' 
#' @title **PRIMARY INTERFACE** - Flexible Covariate Processing
#' @description Processes covariate aggregation specifications using a modular 4-pattern system:
#' \itemize{
#'   \item **Pattern 1**: Aggregate over arbitrary periods → Single variable
#'   \item **Pattern 2**: Aggregate over last N / first M periods → Single variable
#'   \item **Pattern 3**: Each period for specified periods → Multiple variables
#'   \item **Pattern 4**: Each period for last N / first M periods → Multiple variables
#' }
#' 
#' @param data Data frame or data.table in long format
#' @param covagg List of covariate aggregation specifications with new flexible format:
#'   \itemize{
#'     \item \code{var} - Required: variable name
#'     \item \code{periods} - Optional: specific periods (vector or range)
#'     \item \code{last} - Optional: last N periods
#'     \item \code{first} - Optional: first M periods
#'     \item \code{each} - Optional: create individual variables (default: FALSE)
#'     \item \code{agg_fun} - Optional: aggregation function (default: mean)
#'     \item \code{label} - Optional: custom variable name
#'   }
#' @param period.pre Vector of pre-treatment periods
#' @param id.var Character name of ID variable
#' @param time.var Character name of time variable
#' @return List of processed covariate matrices
#' @importFrom data.table as.data.table is.data.table copy setorderv setnames
#' @export
#' 
#' @examples
#' \dontrun{
#' # Pattern 1: Aggregate over arbitrary periods
#' covagg1 = list(gdp_custom = list(var = "gdp", periods = c(1980, 1985, 1990)))
#' 
#' # Pattern 2: Aggregate over last N periods
#' covagg2 = list(gdp_last10 = list(var = "gdp", last = 10, agg_fun = median))
#' 
#' # Pattern 3: Each period for specified periods
#' covagg3 = list(gdp_each = list(var = "gdp", periods = c(1980, 1985), each = TRUE))
#' 
#' # Pattern 4: Each period for last N periods
#' covagg4 = list(gdp_each_last = list(var = "gdp", last = 5, each = TRUE))
#' }
process_covariates <- function(data, covagg, period.pre, id.var, time.var) {
  
  # Convert to data.table if not already
  if (!data.table::is.data.table(data)) {
    dt <- data.table::as.data.table(data)
  } else {
    dt <- data.table::copy(data)
  }
  
  # Filter to pre-period only
  dt_pre <- dt[get(time.var) %in% period.pre]
  
  processed_covs <- list()
  
  for (cov_name in names(covagg)) {
    cov_spec <- covagg[[cov_name]]
    
    # Handle nested covagg structures (like those with label)
    if (is.list(cov_spec) && is.null(cov_spec$var) && !is.null(cov_spec$label)) {
      # This is a nested structure - process each sub-specification
      for (sub_name in names(cov_spec)) {
        sub_spec <- cov_spec[[sub_name]]
        
        # Skip the label element and only process actual covariate specs
        if (!is.list(sub_spec) || is.null(sub_spec$var) || sub_name == "label") {
          next
        }
        
        # Process this sub-specification with combined name
        combined_name <- paste(cov_name, sub_name, sep = "_")
        
        # VALIDATE SUB-COVARIATE SPECIFICATION
        validate_covariate_spec(sub_spec, combined_name, period.pre)
        
        var_name <- sub_spec$var
        
        # Check if variable exists
        if (!var_name %in% names(dt)) {
          warning(paste("Variable", var_name, "not found in data. Skipping."))
          next
        }
        
        # Extract relevant columns
        dt_var <- dt_pre[, .SD, .SDcols = c(id.var, time.var, var_name)]
        
        # DETERMINE PROCESSING PATTERN
        pattern <- determine_processing_pattern(sub_spec)
        
        # ROUTE TO APPROPRIATE PROCESSOR
        tryCatch({
          switch(pattern,
            "aggregate_arbitrary" = {
              result <- process_aggregate_arbitrary(dt_var, sub_spec, period.pre, id.var, time.var, var_name)
              processed_covs[[combined_name]] <- result
            },
            "aggregate_range" = {
              result <- process_aggregate_range(dt_var, sub_spec, period.pre, id.var, time.var, var_name)
              processed_covs[[combined_name]] <- result
            },
            "each_arbitrary" = {
              result <- process_each_arbitrary(dt_var, sub_spec, period.pre, id.var, time.var, var_name)
              processed_covs <- append(processed_covs, result)
            },
            "each_range" = {
              result <- process_each_range(dt_var, sub_spec, period.pre, id.var, time.var, var_name)
              processed_covs <- append(processed_covs, result)
            },
            stop(paste("Unknown processing pattern:", pattern))
          )
        }, error = function(e) {
          warning(paste("Failed to process covariate", combined_name, ":", e$message))
        })
      }
      next  # Skip to next main covagg element
    }
    
    # Skip non-list elements and specs without var (flat structure)
    if (!is.list(cov_spec) || is.null(cov_spec$var)) {
      next
    }
    
    # VALIDATE COVARIATE SPECIFICATION
    validate_covariate_spec(cov_spec, cov_name, period.pre)
    
    var_name <- cov_spec$var
    
    # Check if variable exists
    if (!var_name %in% names(dt)) {
      warning(paste("Variable", var_name, "not found in data. Skipping."))
      next
    }
    
    # Extract relevant columns
    dt_var <- dt_pre[, .SD, .SDcols = c(id.var, time.var, var_name)]
    
    # DETERMINE PROCESSING PATTERN
    pattern <- determine_processing_pattern(cov_spec)
    
    # ROUTE TO APPROPRIATE PROCESSOR
    tryCatch({
      switch(pattern,
        "aggregate_arbitrary" = {
          result <- process_aggregate_arbitrary(dt_var, cov_spec, period.pre, id.var, time.var, var_name)
          processed_covs[[cov_name]] <- result
        },
        "aggregate_range" = {
          result <- process_aggregate_range(dt_var, cov_spec, period.pre, id.var, time.var, var_name)
          processed_covs[[cov_name]] <- result
        },
        "each_arbitrary" = {
          result <- process_each_arbitrary(dt_var, cov_spec, period.pre, id.var, time.var, var_name)
          processed_covs <- append(processed_covs, result)
        },
        "each_range" = {
          result <- process_each_range(dt_var, cov_spec, period.pre, id.var, time.var, var_name)
          processed_covs <- append(processed_covs, result)
        },
        stop(paste("Unknown processing pattern:", pattern))
      )
    }, error = function(e) {
      warning(paste("Failed to process covariate", cov_name, ":", e$message))
    })
  }
  
  return(processed_covs)
}

#' **VALIDATION FUNCTIONS**

#' Validate Covariate Specification
#' @keywords internal
validate_covariate_spec <- function(cov_spec, cov_name, period.pre) {
  
  # Check for required var parameter
  if (is.null(cov_spec$var) || !is.character(cov_spec$var) || length(cov_spec$var) != 1) {
    stop(paste("Covariate specification", cov_name, "must have a single 'var' parameter"))
  }
  
  # Check for conflicting period specifications
  has_periods <- !is.null(cov_spec$periods)
  has_last <- !is.null(cov_spec$last)
  has_first <- !is.null(cov_spec$first)
  
  period_specs <- sum(has_periods, has_last, has_first)
  if (period_specs > 1) {
    stop(paste("Covariate specification", cov_name, "cannot have multiple period specifications (periods, last, first)"))
  }
  
  # Check for duplicative ranges
  if (has_first && has_last) {
    first_periods <- head(sort(period.pre), cov_spec$first)
    last_periods <- tail(sort(period.pre), cov_spec$last)
    
    if (length(intersect(first_periods, last_periods)) > 0) {
      warning(paste("Covariate specification", cov_name, "has overlapping first and last periods"))
    }
  }
  
  # Validate aggregation function if provided
  if (!is.null(cov_spec$agg_fun)) {
    validate_agg_function(cov_spec$agg_fun, cov_name)
  }
  
  # Validate period ranges
  if (has_periods) {
    validate_periods(cov_spec$periods, period.pre, cov_name)
  }
  
  if (has_last) {
    validate_range_size(cov_spec$last, length(period.pre), cov_name, "last")
  }
  
  if (has_first) {
    validate_range_size(cov_spec$first, length(period.pre), cov_name, "first")
  }
}

#' Validate Aggregation Function
#' @keywords internal
validate_agg_function <- function(agg_fun, cov_name) {
  if (!is.function(agg_fun)) {
    stop(paste("agg_fun for", cov_name, "must be a function"))
  }
  
  # Test the function with sample data
  test_data <- c(1, 2, 3, 4, 5)
  
  tryCatch({
    result <- agg_fun(test_data)
    
    if (length(result) != 1) {
      stop(paste("agg_fun for", cov_name, "must return a single value, got", length(result), "values"))
    }
    
    if (!is.numeric(result) || !is.finite(result)) {
      stop(paste("agg_fun for", cov_name, "must return a finite numeric value"))
    }
    
  }, error = function(e) {
    stop(paste("agg_fun for", cov_name, "failed test:", e$message))
  })
}

#' Validate Periods
#' @keywords internal
validate_periods <- function(periods, period.pre, cov_name) {
  if (!is.numeric(periods) && !is.integer(periods)) {
    stop(paste("periods for", cov_name, "must be numeric or integer"))
  }
  
  missing_periods <- setdiff(periods, period.pre)
  if (length(missing_periods) > 0) {
    warning(paste("Covariate", cov_name, "requests periods not in data:", paste(missing_periods, collapse = ", ")))
  }
}

#' Validate Range Size
#' @keywords internal
validate_range_size <- function(range_size, total_periods, cov_name, range_type) {
  if (!is.numeric(range_size) || length(range_size) != 1 || range_size <= 0) {
    stop(paste(range_type, "for", cov_name, "must be a positive number"))
  }
  
  if (range_size > total_periods) {
    warning(paste("Covariate", cov_name, "requests", range_size, range_type, "periods but only", total_periods, "available"))
  }
}

#' **PATTERN DETECTION**

#' Determine Processing Pattern from Covariate Specification
#' @keywords internal
determine_processing_pattern <- function(cov_spec) {
  
  has_periods <- !is.null(cov_spec$periods)
  has_range <- !is.null(cov_spec$last) || !is.null(cov_spec$first)
  has_each <- !is.null(cov_spec$each) && cov_spec$each == TRUE
  
  # Determine pattern based on combinations
  if (has_periods && has_each) {
    return("each_arbitrary")
  } else if (has_range && has_each) {
    return("each_range")
  } else if (has_periods && !has_each) {
    return("aggregate_arbitrary")
  } else if (has_range && !has_each) {
    return("aggregate_range")
  } else if (!has_periods && !has_range && !has_each) {
    # Default: aggregate over all periods
    return("aggregate_range")
  } else {
    stop("Invalid covariate specification combination")
  }
}

#' **HELPER FUNCTIONS**

#' Resolve Target Periods from Specification
#' @keywords internal
resolve_target_periods <- function(cov_spec, period.pre) {
  
  if (!is.null(cov_spec$periods)) {
    # Use specified periods, filter to available ones
    return(intersect(cov_spec$periods, period.pre))
    
  } else if (!is.null(cov_spec$last)) {
    # Use last N periods
    return(tail(sort(period.pre), min(cov_spec$last, length(period.pre))))
    
  } else if (!is.null(cov_spec$first)) {
    # Use first N periods
    return(head(sort(period.pre), min(cov_spec$first, length(period.pre))))
    
  } else {
    # Default: use all periods
    return(period.pre)
  }
}

#' Get Aggregation Function
#' @keywords internal
get_agg_function <- function(cov_spec) {
  if (!is.null(cov_spec$agg_fun)) {
    return(cov_spec$agg_fun)
  } else {
    return(mean)  # Default aggregation function
  }
}

#' Generate Variable Name
#' @keywords internal
generate_var_name <- function(cov_spec, cov_name, period = NULL) {
  
  # Use custom label if provided
  if (!is.null(cov_spec$label)) {
    base_name <- cov_spec$label
  } else {
    base_name <- cov_name
  }
  
  # Add period suffix if specified
  if (!is.null(period)) {
    return(paste0(base_name, "_", period))
  } else {
    return(base_name)
  }
}

#' **CORE PROCESSING FUNCTIONS**

#' Process Pattern 1: Aggregate Over Arbitrary Periods
#' @keywords internal
process_aggregate_arbitrary <- function(dt_var, cov_spec, period.pre, id.var, time.var, var_name) {
  
  # Get target periods
  target_periods <- resolve_target_periods(cov_spec, period.pre)
  
  if (length(target_periods) == 0) {
    warning(paste("No valid periods found for", var_name))
    return(data.table::data.table())
  }
  
  # Filter to target periods
  dt_filtered <- dt_var[get(time.var) %in% target_periods]
  
  # Get aggregation function
  agg_fun <- get_agg_function(cov_spec)
  
  # Apply aggregation
  result <- dt_filtered[!is.na(get(var_name)), 
                       .(agg_value = agg_fun(get(var_name))), 
                       by = get(id.var)]
  
  # Set column names
  final_var_name <- generate_var_name(cov_spec, var_name)
  data.table::setnames(result, c(id.var, final_var_name))
  
  return(result)
}

#' Process Pattern 2: Aggregate Over Range (Last N / First M)
#' @keywords internal
process_aggregate_range <- function(dt_var, cov_spec, period.pre, id.var, time.var, var_name) {
  
  # Get target periods
  target_periods <- resolve_target_periods(cov_spec, period.pre)
  
  if (length(target_periods) == 0) {
    warning(paste("No valid periods found for", var_name))
    return(data.table::data.table())
  }
  
  # Filter to target periods
  dt_filtered <- dt_var[get(time.var) %in% target_periods]
  
  # Get aggregation function
  agg_fun <- get_agg_function(cov_spec)
  
  # Apply aggregation
  result <- dt_filtered[!is.na(get(var_name)), 
                       .(agg_value = agg_fun(get(var_name))), 
                       by = get(id.var)]
  
  # Set column names
  final_var_name <- generate_var_name(cov_spec, var_name)
  data.table::setnames(result, c(id.var, final_var_name))
  
  return(result)
}

#' Process Pattern 3: Each Period for Specified Periods
#' @keywords internal
process_each_arbitrary <- function(dt_var, cov_spec, period.pre, id.var, time.var, var_name) {
  
  # Get target periods
  target_periods <- resolve_target_periods(cov_spec, period.pre)
  
  if (length(target_periods) == 0) {
    warning(paste("No valid periods found for", var_name))
    return(list())
  }
  
  result <- list()
  
  for (period in target_periods) {
    period_data <- dt_var[get(time.var) == period]
    
    if (nrow(period_data) > 0) {
      # Generate variable name for this period
      period_var_name <- generate_var_name(cov_spec, var_name, period)
      
      # Clean data (remove NAs)
      period_clean <- period_data[!is.na(get(var_name))]
      
      if (nrow(period_clean) > 0) {
        # Create result with proper column names
        period_result <- period_clean[, .(id = get(id.var), value = get(var_name))]
        data.table::setnames(period_result, c(id.var, period_var_name))
        result[[period_var_name]] <- period_result
      }
    }
  }
  
  return(result)
}

#' Process Pattern 4: Each Period for Range (Last N / First M)
#' @keywords internal
process_each_range <- function(dt_var, cov_spec, period.pre, id.var, time.var, var_name) {
  
  # Get target periods
  target_periods <- resolve_target_periods(cov_spec, period.pre)
  
  if (length(target_periods) == 0) {
    warning(paste("No valid periods found for", var_name))
    return(list())
  }
  
  # Check for potentially large number of variables
  if (length(target_periods) > 50) {
    warning(paste("Creating", length(target_periods), "variables for", var_name, "- this may use significant memory"))
  }
  
  result <- list()
  
  for (period in target_periods) {
    period_data <- dt_var[get(time.var) == period]
    
    if (nrow(period_data) > 0) {
      # Generate variable name for this period
      period_var_name <- generate_var_name(cov_spec, var_name, period)
      
      # Clean data (remove NAs)
      period_clean <- period_data[!is.na(get(var_name))]
      
      if (nrow(period_clean) > 0) {
        # Create result with proper column names
        period_result <- period_clean[, .(id = get(id.var), value = get(var_name))]
        data.table::setnames(period_result, c(id.var, period_var_name))
        result[[period_var_name]] <- period_result
      }
    }
  }
  
  return(result)
}



