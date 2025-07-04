#' Calculate Abadie-style Significance Measures
#'
#' @title Calculate Multiple Test Statistics and Rank-based P-values
#' @description Calculates significance measures following Abadie et al. (2010) German reunification paper.
#' Computes post/pre RMSPE ratios, treatment effects, normalized treatment effects, and their rank-based two-sided p-values.
#'
#' @param taus_result Data.table containing treatment effects for all units
#' @param rmse_result Data.table containing RMSE values for all units
#' @param rmspe_threshold Numeric. Threshold for filtering poor pre-treatment fits (default: 0)
#'
#' @return List containing Abadie significance measures for all test statistics
calculate_abadie_significance <- function(taus_result, rmse_result, rmspe_threshold = 0) {
  
  if (nrow(taus_result) == 0 || nrow(rmse_result) == 0) {
    return(list(
      rmse_ratio = list(test_statistics = data.table(), p_values = data.table(), filtered_results = data.table()),
      treatment_effect = list(test_statistics = data.table(), p_values = data.table()),
      normalized_te = list(test_statistics = data.table(), p_values = data.table())
    ))
  }

  # Ensure inputs are data.tables for efficient operations
  if (!is.data.table(taus_result)) setDT(taus_result)
  if (!is.data.table(rmse_result)) setDT(rmse_result)

  # Calculate all test statistics
  
  # 1. RMSE Ratio (post/pre RMSPE)
  post_rmspe <- taus_result[post_period == TRUE, .(
    post_rmspe = sqrt(mean(tau^2))
  ), by = .(unit_name, unit_type, outcome_model)]
  
  rmse_ratio_stats <- merge(rmse_result, post_rmspe, by = c("unit_name", "unit_type", "outcome_model"))
  rmse_ratio_stats[, test_statistic_value := post_rmspe / pre_rmse]
  
  # 2. Treatment Effect (average post-period treatment effect)
  te_stats <- taus_result[post_period == TRUE, .(
    test_statistic_value = mean(tau)
  ), by = .(unit_name, unit_type, outcome_model)]
  
  # 3. Normalized Treatment Effect (TE / pre-period SD)
  te_post <- taus_result[post_period == TRUE, .(
    avg_treatment_effect = mean(tau)
  ), by = .(unit_name, unit_type, outcome_model)]
  
  pre_sd <- taus_result[post_period == FALSE, .(
    pre_sd = sd(tau)
  ), by = .(unit_name, unit_type, outcome_model)]
  
  normalized_te_stats <- merge(te_post, pre_sd, by = c("unit_name", "unit_type", "outcome_model"))
  normalized_te_stats[, test_statistic_value := avg_treatment_effect / pre_sd]
  
  # Function to calculate p-values for a given test statistic
  calculate_pvalues <- function(test_stat_data, test_type = "treatment_effect") {
    treated_ratios <- test_stat_data[unit_type == "treated", .(outcome_model, treated_ratio = test_statistic_value)]
    
    if (nrow(treated_ratios) > 0) {
      p_values_result <- test_stat_data[unit_type == "control", {
        if (.N > 0) {
          treated_val <- treated_ratios[outcome_model == .BY$outcome_model, treated_ratio]
          if (length(treated_val) > 0) {
            control_vals <- test_statistic_value
            
            if (test_type == "rmse_ratio") {
              # For RMSE ratios: one-sided test (is treated unit an outlier in upper tail?)
              # Traditional Abadie approach - testing if treated unit is unusually volatile
              rank_val <- sum(control_vals >= treated_val) + 1
              total_val <- length(control_vals) + 1
              p_value_one_sided <- rank_val / total_val
              p_value_two_sided <- NA  # Keep as one-sided for RMSE ratios
              
            } else {
              # For treatment effects: standard two-sided test
              total_val <- length(control_vals) + 1
              
              # Count control units as extreme or more extreme in the same direction
              if (treated_val >= 0) {
                same_direction_extreme <- sum(control_vals >= treated_val)
              } else {
                same_direction_extreme <- sum(control_vals <= treated_val)  
              }
              
              # One-sided p-value
              p_value_one_sided <- (same_direction_extreme + 1) / total_val
              
              # Two-sided: double the one-sided p-value, capped at 1
              p_value_two_sided <- min(2 * p_value_one_sided, 1.0)
              
              rank_val <- same_direction_extreme + 1
            }
            
            list(
              treated_ratio = treated_val,
              rank = rank_val,
              total_units = total_val,
              p_value_one_sided = p_value_one_sided,
              p_value = p_value_two_sided
            )
          } else {
            list(treated_ratio = NA_real_, rank = NA_integer_, total_units = NA_integer_, 
                 p_value_one_sided = NA_real_, p_value = NA_real_)
          }
        } else {
          list(treated_ratio = NA_real_, rank = NA_integer_, total_units = NA_integer_, 
               p_value_one_sided = NA_real_, p_value = NA_real_)
        }
      }, by = outcome_model]
    } else {
      p_values_result <- data.table()
    }
    return(p_values_result)
  }
  
  # Calculate p-values for all test statistics with appropriate test types
  rmse_ratio_pvalues <- calculate_pvalues(rmse_ratio_stats, "rmse_ratio")
  te_pvalues <- calculate_pvalues(te_stats, "treatment_effect")
  normalized_te_pvalues <- calculate_pvalues(normalized_te_stats, "treatment_effect")
  
  # Calculate filtered results for RMSE ratio (only applicable for this test statistic)
  treated_pre_rmse <- rmse_ratio_stats[unit_type == "treated", .(outcome_model, treated_pre_rmse = pre_rmse)]
  
  if (nrow(treated_pre_rmse) > 0) {
    filtered_result <- rmse_ratio_stats[, {
        treated_rmse <- treated_pre_rmse[outcome_model == .BY$outcome_model, treated_pre_rmse]
        if (length(treated_rmse) > 0) {
          threshold_val <- rmspe_threshold * treated_rmse
          
          # Filter data: keep treated + controls with good pre-treatment fit
          filtered_data <- .SD[unit_type == "treated" | (unit_type == "control" & pre_rmse <= threshold_val)]
        
          # Calculate filtered p-value if we have both treated and control units
          if (nrow(filtered_data[unit_type == "treated"]) > 0 && nrow(filtered_data[unit_type == "control"]) > 0) {
            treated_ratio_val <- filtered_data[unit_type == "treated", test_statistic_value]
            control_ratios_val <- filtered_data[unit_type == "control", test_statistic_value]
            
            # For filtered RMSE ratios: use one-sided test (consistent with unfiltered)
            rank_filtered_val <- sum(control_ratios_val >= treated_ratio_val) + 1
            total_filtered_val <- length(control_ratios_val) + 1
            units_excluded_val <- nrow(.SD[unit_type == "control"]) - length(control_ratios_val)
            
            p_value_filtered_one_sided <- rank_filtered_val / total_filtered_val
            p_value_filtered_two_sided <- p_value_filtered_one_sided  # Keep as one-sided
            
            list(
              treated_ratio = treated_ratio_val,
              rank_filtered = rank_filtered_val,
              total_units_filtered = total_filtered_val,
              p_value_filtered_one_sided = p_value_filtered_one_sided,
              p_value_filtered = p_value_filtered_two_sided,
              rmspe_threshold = rmspe_threshold,
              units_excluded = units_excluded_val
            )
          } else {
            list(treated_ratio = NA_real_, rank_filtered = NA_integer_, 
                 total_units_filtered = NA_integer_, p_value_filtered_one_sided = NA_real_,
                 p_value_filtered = NA_real_, rmspe_threshold = rmspe_threshold, units_excluded = NA_integer_)
          }
        } else {
          list(treated_ratio = NA_real_, rank_filtered = NA_integer_, 
               total_units_filtered = NA_integer_, p_value_filtered_one_sided = NA_real_,
               p_value_filtered = NA_real_, rmspe_threshold = rmspe_threshold, units_excluded = NA_integer_)
        }
      }, by = outcome_model]
  } else {
    filtered_result <- data.table()
  }
  
  return(list(
    rmse_ratio = list(
      test_statistics = rmse_ratio_stats,
      p_values = rmse_ratio_pvalues,
      filtered_results = filtered_result
    ),
    treatment_effect = list(
      test_statistics = te_stats,
      p_values = te_pvalues
    ),
    normalized_te = list(
      test_statistics = normalized_te_stats,
      p_values = normalized_te_pvalues
    )
  ))
}


#' Placebo Inference for Synthetic Control Method
#'
#' @title Perform placebo inference with multiple test statistics
#' @description Performs placebo inference by estimating synthetic control for each control unit
#' and computing three test statistics: post/pre RMSE ratio, treatment effect, and normalized treatment effect.
#' Returns two-sided p-values for all test statistics.
#'
#' @param sc.pred List. Results from synthetic control estimation
#' @param dataset Data.frame. Original dataset used for estimation
#' @param cores Numeric. Number of cores for parallel processing
#' @param verbose Logical. Whether to print progress information
#'
#' @return List containing:
#' \itemize{
#'   \item taus: Treatment effects for all units
#'   \item rmse: RMSE values for all units
#'   \item abadie_significance: List with three elements:
#'     \itemize{
#'       \item rmse_ratio: Post/pre RMSE ratio test statistic and p-values
#'       \item treatment_effect: Treatment effect test statistic and p-values
#'       \item normalized_te: Normalized treatment effect test statistic and p-values
#'     }
#' }
#'
#' @export
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