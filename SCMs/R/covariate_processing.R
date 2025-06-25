#' Process Covariate Aggregations for Synthetic Control
#' 
#' @param data Data frame or data.table in long format
#' @param covagg List of covariate aggregation specifications
#' @param period.pre Vector of pre-treatment periods
#' @param id.var Character name of ID variable
#' @param time.var Character name of time variable
#' @return List of processed covariate matrices
#' @importFrom data.table as.data.table is.data.table copy setorderv setnames frollmean shift
#' @keywords internal
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
    
    # Skip non-list elements (like label) and only process actual covariate specs
    if (!is.list(cov_spec) || is.null(cov_spec$var)) {
      next
    }
    
    var_name <- cov_spec$var
    
    # Check if variable exists
    if (!var_name %in% names(dt)) {
      warning(paste("Variable", var_name, "not found in data. Skipping."))
      next
    }
    
    # Extract relevant columns
    dt_var <- dt_pre[, .SD, .SDcols = c(id.var, time.var, var_name)]
    
    # Process based on specification type
    if (!is.null(cov_spec$each_period) && cov_spec$each_period) {
      # EACH PERIOD SEPARATELY
      processed_covs <- append(processed_covs, 
                               process_each_period_dt(dt_var, cov_name, period.pre, id.var, time.var, var_name))
      
    } else if (!is.null(cov_spec$periods)) {
      # Specific periods
      processed_covs[[cov_name]] <- process_specific_periods_dt(dt_var, cov_spec$periods, id.var, time.var, var_name)
      
    } else if (!is.null(cov_spec$every_n)) {
      # Every N periods
      processed_covs[[cov_name]] <- process_every_n_periods_dt(dt_var, cov_spec$every_n, period.pre, id.var, time.var, var_name)
      
    } else if (!is.null(cov_spec$average)) {
      # Averages
      processed_covs[[cov_name]] <- process_averages_dt(dt_var, cov_spec$average, id.var, var_name)
      
    } else if (!is.null(cov_spec$first)) {
      # First N periods
      processed_covs[[cov_name]] <- process_first_n_dt(dt_var, cov_spec$first, period.pre, id.var, time.var, var_name)
      
    } else if (!is.null(cov_spec$last)) {
      # Last N periods  
      processed_covs[[cov_name]] <- process_last_n_dt(dt_var, cov_spec$last, period.pre, id.var, time.var, var_name)
      
    } else if (!is.null(cov_spec$rolling)) {
      # Rolling windows
      processed_covs <- append(processed_covs,
                               process_rolling_dt(dt_var, cov_spec$rolling, cov_name, id.var, time.var, var_name))
      
    } else if (!is.null(cov_spec$ranges)) {
      # Custom ranges
      processed_covs <- append(processed_covs,
                               process_ranges_dt(dt_var, cov_spec$ranges, cov_name, id.var, time.var, var_name))
      
    } else if (!is.null(cov_spec$growth)) {
      # Growth rates
      processed_covs[[cov_name]] <- process_growth_dt(dt_var, cov_spec$growth, id.var, time.var, var_name)
      
    } else if (!is.null(cov_spec$volatility)) {
      # Volatility measures
      processed_covs[[cov_name]] <- process_volatility_dt(dt_var, cov_spec$volatility, id.var, time.var, var_name)
      
    } else {
      # Default: use variable as-is (backward compatibility)
      processed_covs[[cov_name]] <- dt_var[!is.na(get(var_name))]
    }
  }
  
  return(processed_covs)
}

#' Process Each Period as Separate Matching Variable
#' @keywords internal
process_each_period_dt <- function(dt_var, cov_name, period.pre, id.var, time.var, var_name) {
  result <- list()
  
  for (period in period.pre) {
    period_data <- dt_var[get(time.var) == period]
    
    if (nrow(period_data) > 0) {
      # Create wide format for this period
      var_name_period <- paste0(cov_name, "_", period)
      
      # Remove rows with NA for this variable
      period_clean <- period_data[!is.na(get(var_name))]
      
      if (nrow(period_clean) > 0) {
        # Return as named list element
        result[[var_name_period]] <- period_clean[, .(id = get(id.var), value = get(var_name))]
        data.table::setnames(result[[var_name_period]], c(id.var, var_name_period))
      }
    }
  }
  
  return(result)
}

#' Process Specific Periods
#' @keywords internal
process_specific_periods_dt <- function(dt_var, periods, id.var, time.var, var_name) {
  result <- dt_var[get(time.var) %in% periods]
  return(result[!is.na(get(var_name))])
}

#' Process Every N Periods
#' @keywords internal
process_every_n_periods_dt <- function(dt_var, n, period.pre, id.var, time.var, var_name) {
  periods_to_use <- period.pre[seq(1, length(period.pre), by = n)]
  result <- dt_var[get(time.var) %in% periods_to_use]
  return(result[!is.na(get(var_name))])
}

#' Process Averages
#' @keywords internal
process_averages_dt <- function(dt_var, avg_type, id.var, var_name) {
  if (avg_type == "full_pre") {
    # Super efficient aggregation
    result <- dt_var[!is.na(get(var_name)), 
                     .(avg_value = mean(get(var_name), na.rm = TRUE)), 
                     by = get(id.var)]
    
    data.table::setnames(result, c(id.var, paste0(var_name, "_avg")))
    return(result)
  }
}

#' Process First N Periods
#' @keywords internal
process_first_n_dt <- function(dt_var, n, period.pre, id.var, time.var, var_name) {
  first_periods <- head(sort(period.pre), n)
  result <- dt_var[get(time.var) %in% first_periods]
  return(result[!is.na(get(var_name))])
}

#' Process Last N Periods
#' @keywords internal
process_last_n_dt <- function(dt_var, n, period.pre, id.var, time.var, var_name) {
  last_periods <- tail(sort(period.pre), n)
  result <- dt_var[get(time.var) %in% last_periods]
  return(result[!is.na(get(var_name))])
}

#' Process Rolling Windows
#' @keywords internal
process_rolling_dt <- function(dt_var, window_size, cov_name, id.var, time.var, var_name) {
  result <- list()
  
  # Sort by time within each unit
  data.table::setorderv(dt_var, c(id.var, time.var))
  
  # Calculate rolling averages using data.table's powerful grouping
  dt_var[, rolling_avg := data.table::frollmean(get(var_name), window_size, na.rm = TRUE), by = get(id.var)]
  
  # Create separate variable for each complete rolling window
  complete_windows <- dt_var[!is.na(rolling_avg)]
  
  if (nrow(complete_windows) > 0) {
    var_name_rolling <- paste0(cov_name, "_roll", window_size)
    result[[var_name_rolling]] <- complete_windows[, .(id = get(id.var), 
                                                        time = get(time.var),
                                                        value = rolling_avg)]
    data.table::setnames(result[[var_name_rolling]], c(id.var, time.var, var_name_rolling))
  }
  
  return(result)
}

#' Process Custom Ranges
#' @keywords internal
process_ranges_dt <- function(dt_var, ranges, cov_name, id.var, time.var, var_name) {
  result <- list()
  
  for (i in seq_along(ranges)) {
    range_periods <- ranges[[i]]
    range_data <- dt_var[get(time.var) %in% range_periods]
    
    if (nrow(range_data) > 0) {
      # Average over this range
      var_name_range <- paste0(cov_name, "_range", i)
      range_avg <- range_data[!is.na(get(var_name)), 
                             .(range_value = mean(get(var_name), na.rm = TRUE)), 
                             by = get(id.var)]
      
      data.table::setnames(range_avg, c(id.var, var_name_range))
      result[[var_name_range]] <- range_avg
    }
  }
  
  return(result)
}

#' Process Growth Rates
#' @keywords internal
process_growth_dt <- function(dt_var, growth_type, id.var, time.var, var_name) {
  # Sort by time within each unit
  data.table::setorderv(dt_var, c(id.var, time.var))
  
  if (growth_type == "period_over_period") {
    # Calculate period-over-period growth rates
    dt_var[, growth_rate := (get(var_name) / data.table::shift(get(var_name), 1L) - 1) * 100, 
           by = get(id.var)]
    
    result <- dt_var[!is.na(growth_rate)]
    return(result[, .(id = get(id.var), time = get(time.var), growth = growth_rate)])
  }
}

#' Process Volatility Measures
#' @keywords internal
process_volatility_dt <- function(dt_var, vol_type, id.var, time.var, var_name) {
  if (vol_type == "sd") {
    result <- dt_var[!is.na(get(var_name)), 
                     .(volatility = sd(get(var_name), na.rm = TRUE)), 
                     by = get(id.var)]
    
    data.table::setnames(result, c(id.var, paste0(var_name, "_vol")))
    return(result)
  } else if (vol_type == "cv") {
    # Coefficient of variation
    result <- dt_var[!is.na(get(var_name)), 
                     .(cv = sd(get(var_name), na.rm = TRUE) / mean(get(var_name), na.rm = TRUE)), 
                     by = get(id.var)]
    
    data.table::setnames(result, c(id.var, paste0(var_name, "_cv")))
    return(result)
  }
}