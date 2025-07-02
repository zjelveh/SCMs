#' Calculate Abadie-style Significance Measures
#'
#' @title Calculate Post/Pre RMSPE Ratios and Rank-based P-values
#' @description Calculates significance measures following Abadie et al. (2010) German reunification paper.
#' Includes post/pre RMSPE ratios, rank-based p-values, and filtering for poor pre-treatment fits.
#'
#' @param taus_result Data.table containing treatment effects for all units
#' @param rmse_result Data.table containing RMSE values for all units
#' @param rmspe_threshold Numeric. Threshold for filtering poor pre-treatment fits (default: 0)
#'
#' @return List containing Abadie significance measures
calculate_abadie_significance <- function(taus_result, rmse_result, rmspe_threshold = 0) {
  
  if (nrow(taus_result) == 0 || nrow(rmse_result) == 0) {
    return(list(
      post_pre_ratios = data.table(),
      p_values = data.table(),
      filtered_results = data.table()
    ))
  }

  
  # Ensure inputs are data.tables for efficient operations
  if (!is.data.table(taus_result)) setDT(taus_result)
  if (!is.data.table(rmse_result)) setDT(rmse_result)

  # Calculate post-period RMSPE by outcome model
  post_rmspe <- taus_result[post_period == TRUE, .(
    post_rmspe = sqrt(mean(tau^2))
  ), by = .(unit_name, unit_type, outcome_model)]
  
  # Merge pre-RMSE with post-RMSPE by outcome model
  rmspe_combined <- merge(rmse_result, post_rmspe, by = c("unit_name", "unit_type", "outcome_model"))
  
  # Calculate post/pre RMSPE ratios
  rmspe_combined[, post_pre_ratio := post_rmspe / pre_rmse]
  
  # Get treated unit ratios for each outcome model
  treated_ratios <- rmspe_combined[unit_type == "treated", .(outcome_model, treated_ratio = post_pre_ratio)]
  
  if (nrow(treated_ratios) > 0) {
    # Calculate ranks and p-values - group by outcome model
    p_values_result <- rmspe_combined[unit_type == "control", {
      if (.N > 0) {
        treated_val <- treated_ratios[outcome_model == .BY$outcome_model, treated_ratio]
        if (length(treated_val) > 0) {
          # Rank-based p-value calculation (consistent with bootstrap approach)
          sorted_control_ratios <- sort(post_pre_ratio)
          rank_val <- findInterval(treated_val, sorted_control_ratios, rightmost.closed = TRUE) + 1
          total_val <- .N + 1
          
          list(
            treated_ratio = treated_val,
            rank = rank_val,
            total_units = total_val,
            p_value = rank_val / total_val
          )
        } else {
          list(treated_ratio = NA_real_, rank = NA_integer_, total_units = NA_integer_, p_value = NA_real_)
        }
      } else {
        list(treated_ratio = NA_real_, rank = NA_integer_, total_units = NA_integer_, p_value = NA_real_)
      }
    }, by = outcome_model]
  } else {
    p_values_result <- data.table()
  }
  
  # Apply RMSPE filtering using vectorized operations
  # Get treated unit pre-RMSE values for thresholding by outcome model
  treated_pre_rmse <- rmspe_combined[unit_type == "treated", .(outcome_model, treated_pre_rmse = pre_rmse)]
  
  if (nrow(treated_pre_rmse) > 0) {
    # Calculate filtered p-values - group by outcome model
    filtered_result <- rmspe_combined[, {
        treated_rmse <- treated_pre_rmse[outcome_model == .BY$outcome_model, treated_pre_rmse]
        if (length(treated_rmse) > 0) {
          threshold_val <- rmspe_threshold * treated_rmse
          
          # Filter data: keep treated + controls with good pre-treatment fit
          filtered_data <- .SD[unit_type == "treated" | (unit_type == "control" & pre_rmse <= threshold_val)]
        
          # Calculate filtered p-value if we have both treated and control units
          if (nrow(filtered_data[unit_type == "treated"]) > 0 && nrow(filtered_data[unit_type == "control"]) > 0) {
            treated_ratio_val <- filtered_data[unit_type == "treated", post_pre_ratio]
            control_ratios_val <- filtered_data[unit_type == "control", post_pre_ratio]
            
            # Use rank-based approach consistent with unfiltered calculation
            sorted_control_ratios <- sort(control_ratios_val)
            rank_filtered_val <- findInterval(treated_ratio_val, sorted_control_ratios, rightmost.closed = TRUE) + 1
            total_filtered_val <- length(control_ratios_val) + 1
            units_excluded_val <- nrow(.SD[unit_type == "control"]) - length(control_ratios_val)
            
            list(
              treated_ratio = treated_ratio_val,
              rank_filtered = rank_filtered_val,
              total_units_filtered = total_filtered_val,
              p_value_filtered = rank_filtered_val / total_filtered_val,
              rmspe_threshold = rmspe_threshold,
              units_excluded = units_excluded_val
            )
          } else {
            list(treated_ratio = NA_real_, rank_filtered = NA_integer_, 
                 total_units_filtered = NA_integer_, p_value_filtered = NA_real_,
                 rmspe_threshold = rmspe_threshold, units_excluded = NA_integer_)
          }
        } else {
          list(treated_ratio = NA_real_, rank_filtered = NA_integer_, 
               total_units_filtered = NA_integer_, p_value_filtered = NA_real_,
               rmspe_threshold = rmspe_threshold, units_excluded = NA_integer_)
        }
      }, by = outcome_model]
  } else {
    filtered_result <- data.table()
  }
  
  return(list(
    post_pre_ratios = rmspe_combined,
    p_values = p_values_result,
    filtered_results = filtered_result
  ))
}


inference_placebo <- function(
    sc.pred,
    dataset,
    cores,
    verbose
){

  # Debug: Check critical fields immediately
  if (is.null(sc.pred)) stop("sc.pred is NULL")
  if (is.null(sc.pred$data)) stop("sc.pred$data is NULL")
  if (is.null(sc.pred$data$specs)) stop("sc.pred$data$specs is NULL")
  if (is.null(sc.pred$data$specs$donors.units)) stop("donors.units is NULL")

  # Get control units and remove treated unit from dataset
  control_units = sc.pred$data$specs$donors.units
  
  # Use the correct column name and treated unit from sc.pred
  col_name_unit <- sc.pred$col_name_unit_name
  name_treated_unit <- sc.pred$name_treated_unit
  
  # Validate inputs
  if (is.null(col_name_unit) || is.null(name_treated_unit)) {
    stop("Missing col_name_unit_name or name_treated_unit in sc.pred")
  }
  
  if (!col_name_unit %in% names(dataset)) {
    stop(paste("Column", col_name_unit, "not found in dataset"))
  }
  
  dataset = dataset[get(col_name_unit) != name_treated_unit]

  control_taus_list = list()
  control_rmse_list = list()
  treated_taus_list = list()
  treated_rmse_list = list()
  
  # Validate required fields
  if (is.null(sc.pred$est.results$outcome_model) || length(sc.pred$est.results$outcome_model) == 0) {
    stop("No outcome models found in sc.pred$est.results$outcome_model")
  }
  
  # Calculate treated unit results once up front - will be reused
  for(oc in names(sc.pred$est.results$outcome_model)){
    # Treated unit calculations
    sc_post = sc.pred$est.results$outcome_model[[oc]]$Y.post.fit 
    sc_pre  = sc.pred$est.results$Y.pre.fit
    
    # Validate dimensions
    if (is.null(sc_post) || length(sc_post) == 0) {
      stop(paste("sc_post is NULL or empty for outcome model:", oc))
    }
    if (is.null(sc_pre) || length(sc_pre) == 0) {
      stop(paste("sc_pre is NULL or empty for outcome model:", oc))
    }
    
    # Calculate treatment effects
    treated_tau = sc.pred$data$Y.post - sc_post 
    treated_pre_tau = sc.pred$data$Y.pre - sc_pre
    
    treated_all_tau = c(treated_pre_tau, treated_tau)
    
    # Calculate RMSE for pre-treatment period of treated unit
    treated_pre_rmse = sqrt(mean(treated_pre_tau^2))
    
    # Validate critical fields before using them
    if (is.null(sc.pred$name_treated_unit) || length(sc.pred$name_treated_unit) == 0) {
      stop("sc.pred$name_treated_unit is NULL or empty")
    }
    if (is.null(sc.pred$min_period) || is.null(sc.pred$end_period)) {
      stop("sc.pred min_period or end_period is NULL")
    }
    if (is.null(sc.pred$treated_period)) {
      stop("sc.pred$treated_period is NULL")
    }
    # Store RMSE for treated unit
    treated_rmse_list[[oc]] = data.table(
      unit_name = sc.pred$name_treated_unit,
      pre_rmse = treated_pre_rmse,
      unit_type = 'treated',
      outcome_model = oc
    )
    
    # Store results for treated unit
    treated_taus_list[[oc]] = data.table(
      unit_name = rep(sc.pred$name_treated_unit, length(treated_all_tau)),
      period = sc.pred$min_period:sc.pred$end_period,
      tau = treated_all_tau,
      post_period = (sc.pred$min_period:sc.pred$end_period) >= sc.pred$treated_period,
      unit_type = 'treated',
      outcome_model = oc
    )
  }

  # Iterate over all control units
  for(cu in control_units){
    # Estimate synthetic control for each control unit
    est_sc <- tryCatch({
      estimate_sc(dataset=dataset,
                  outcome=sc.pred$outcome,
                  covagg=sc.pred$covagg,
                  col_name_unit_name=sc.pred$col_name_unit_name,
                  name_treated_unit=cu,
                  col_name_period=sc.pred$col_name_period,
                  treated_period=sc.pred$treated_period,
                  min_period=sc.pred$min_period,
                  end_period=sc.pred$end_period, 
                  feature_weights = sc.pred$feature_weights,
                  outcome_models = sc.pred$outcome_models,
                  w.constr = sc.pred$w.constr)
    }, error = function(e) {
      if (verbose) {
        warning(paste("Failed to estimate SC for control unit", cu, ":", e$message))
      }
      return(NULL)
    })

    # Skip this control unit if estimation failed
    if (is.null(est_sc)) {
      next
    }

    # Calculate treatment effects for each outcome model
    for(oc in names(est_sc$est.results$outcome_model)){
      sc_post = est_sc$est.results$outcome_model[[oc]]$Y.post.fit        
      sc_pre  = est_sc$est.results$Y.pre.fit
      
      # Calculate treatment effects
      tau = est_sc$data$Y.post - sc_post 
      pre_tau = est_sc$data$Y.pre - sc_pre
      all_tau = c(pre_tau, tau)
      
      # Calculate RMSE for pre-treatment period
      pre_rmse = sqrt(mean(pre_tau^2))
      
      # Store results for control unit
      control_rmse_list[[paste(cu, oc, sep="_")]] = data.table(
        unit_name = cu,
        pre_rmse = pre_rmse,
        unit_type = 'control',
        outcome_model = oc
      )
      
      control_taus_list[[paste(cu, oc, sep="_")]] = data.table(
        unit_name = rep(cu, length(all_tau)),
        period = sc.pred$min_period:sc.pred$end_period,
        tau = all_tau,
        post_period = (sc.pred$min_period:sc.pred$end_period) >= sc.pred$treated_period,
        unit_type = 'control',
        outcome_model = oc
      )
    }
  }
 
  # Combine treated and control results
  all_taus = c(treated_taus_list, control_taus_list)
  all_rmse = c(treated_rmse_list, control_rmse_list)

  # Handle empty lists safely
  if (length(all_taus) == 0) {
    taus_result <- data.table()
  } else {
    taus_result <- rbindlist(all_taus)
  }
  
  if (length(all_rmse) == 0) {
    rmse_result <- data.table()
  } else {
    rmse_result <- rbindlist(all_rmse)
  }

  # Calculate Abadie-style significance measures
  abadie_results <- calculate_abadie_significance(taus_result, rmse_result)

  # Return both treatment effects, RMSE values, and Abadie significance measures
  return(list(
    taus = taus_result,
    rmse = rmse_result,
    abadie_significance = abadie_results
  ))
}