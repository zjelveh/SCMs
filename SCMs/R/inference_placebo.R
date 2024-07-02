#' @title Perform Placebo Inference for Synthetic Control Method
#' @description This function performs placebo inference by applying the synthetic control method to each control unit as if it were treated.
#'
#' @param sc.pred List. Results from the synthetic control estimation for the actual treated unit.
#' @param dataset Data frame. The original dataset used for synthetic control estimation.
#' @param cores Numeric. Number of cores to use for parallel processing (not currently used in the function).
#' @param verbose Logical. Whether to print progress information (not currently used in the function).
#'
#' @return A list of data tables containing treatment effect estimates for each control unit and the actual treated unit.
#'
#' @import data.table
#'
#' @export
#'
#' @examples
#' # Example usage (replace with actual example when available)
#' # placebo_results <- inference_placebo(sc_results, my_data, cores = 1, verbose = TRUE)
inference_placebo <- function(
    sc.pred,
    dataset,
    cores,
    verbose
    ){

    # Get control units and remove treated unit from dataset
    control_units = sc.pred$data$specs$donors.units
    dataset = dataset[get(sc.pred$col_name_unit_name)!=sc.pred$name_treated_unit]
    
    results = list()
    taus = list()
    
    # Iterate over all control units
    for(cu in control_units){
        # Estimate synthetic control for each control unit
        est_sc <- estimate_sc(dataset=dataset,
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
        
        results[[cu]] = est_sc
        
        # Calculate treatment effects for each outcome model
        for(oc in names(est_sc$est.results$outcome_model)){
            sc_post = est_sc$est.results$outcome_model[[oc]]$Y.post.fit        
            sc_pre  = est_sc$est.results$Y.pre.fit
            
            # Calculate treatment effects
            tau = results[[cu]]$data$Y.post - sc_post 
            pre_tau = results[[cu]]$data$Y.pre - sc_pre
            all_tau = c(pre_tau, tau)
            
            # Store results for control unit
            taus[[length(taus) + 1]] = data.table(
                unit_name = rep(cu, length(all_tau)),
                period = sc.pred$min_period:sc.pred$end_period,
                tau=all_tau,
                post_period = (sc.pred$min_period:sc.pred$end_period) >= sc.pred$treated_period,
                unit_type = 'control',
                outcome_model = oc
            )
            
            # Calculate and store results for actual treated unit
            sc_post = sc.pred$est.results$outcome_model[[oc]]$Y.post.fit 
            sc_pre  = sc.pred$est.results$Y.pre.fit
            tau = sc.pred$data$Y.post - sc_post 
            pre_tau = sc.pred$data$Y.pre - sc_pre
            all_tau = c(pre_tau, tau)
            taus[[length(taus) + 1]] = data.table(
                unit_name = rep(sc.pred$name_treated_unit, length(all_tau)),
                period = sc.pred$min_period:sc.pred$end_period,
                tau=all_tau,
                post_period = (sc.pred$min_period:sc.pred$end_period) >= sc.pred$treated_period,
                unit_type = 'treated',
                outcome_model = oc
            )
        }
    }

    return(taus)
}