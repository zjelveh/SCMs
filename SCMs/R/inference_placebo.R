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
  
  # Return both treatment effects and RMSE values
  return(list(
    taus = taus_result,
    rmse = rmse_result
  ))
}