inference_placebo <- function(
    sc.pred,
    dataset,
    cores,
    verbose
    ){


    control_units = sc.pred$data$specs$donors.units
    dataset = dataset[get(sc.pred$col_name_unit_name)!=sc.pred$name_treated_unit]
    results = list()
    taus = list()
    for(cu in control_units){
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
        for(oc in names(est_sc$est.results$outcome_model)){
            sc_post = est_sc$est.results$outcome_model[[oc]]$Y.post.fit        
            # sc_post = results[[cu]]$data$P %*% results[[cu]]$est.results$w
            sc_pre  = est_sc$est.results$Y.pre.fit # fix for augsynth
            tau = results[[cu]]$data$Y.post - sc_post 
            pre_tau = results[[cu]]$data$Y.pre - sc_pre
            all_tau = c(pre_tau, tau)
            taus[[length(taus) + 1]] = data.table(
                unit_name = rep(cu, length(all_tau)),
                period = sc.pred$min_period:sc.pred$end_period,
                tau=all_tau,
                post_period = (sc.pred$min_period:sc.pred$end_period) >= sc.pred$treated_period,
                unit_type = 'control',
                outcome_model = oc
            )
            # Compute for actual treated
            sc_post = sc.pred$est.results$outcome_model[[oc]]$Y.post.fit 
            sc_pre  = sc.pred$est.results$Y.pre.fit # fix for augsynth
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