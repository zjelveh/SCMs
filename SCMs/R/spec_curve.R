most_similar <- function(
    dataset,
    outcome,
    col_name_unit_name,
    name_treated_unit,
    treated_period,
    min_period, 
    col_name_period){
  
  dataset = as.data.table(dataset)
  
  dataset_pre_treatment = dataset[between(get(col_name_period), min_period, (treated_period-1))]
  dataset_pre_treatment_control = dataset_pre_treatment[get(col_name_unit_name)!=name_treated_unit]
  dataset_pre_treatment_treated = dataset_pre_treatment[get(col_name_unit_name)==name_treated_unit] 
  
  ybar_treated = mean(dataset_pre_treatment_treated[[outcome]])
  ybar_controls = dataset_pre_treatment_control[, .(ybar_control=mean(get(outcome))), by=col_name_unit_name]
  ybar_controls[, diff_with_treated := abs(ybar_treated - ybar_control)]
  ybar_controls = ybar_controls[order(diff_with_treated)]
  control_units = ybar_controls[ 1 : round(nrow(ybar_controls) / 2) ][[col_name_unit_name]]
  dataset = dataset[get(col_name_unit_name)%in%c(name_treated_unit, control_units)]
  return(dataset)
}


spec_curve <- function(
    dataset,
    outcomes,
    col_name_unit_name,
    name_treated_unit,
    covagg,
    treated_period,
    min_period,
    end_period, 
    col_name_period,
    feature_weights = c('uniform', 'optimize'),
    num_pre_period_years = NA,
    outcome_models = c('none', 'augsynth', 'ridge', 'lasso', 'ols'),
    donor_sample = c('all', 'most_similar'),
    sim_function = most_similar,
    constraints=list(
        list(name='ols'),
        list(name='simplex'),
        list(name='lasso'),
        list(name='ridge'),
        list(name='L1-L2')
    )
){
    results = list()
    for(outc in outcomes){
        results[[outc]] = list()
            for(const in constraints){
                results[[outc]][[const$name]] = list()
                for(fw in feature_weights){
                    results[[outc]][[const$name]][[fw]] = list()
                    for(ca in covagg){ 
                        feature_names = paste0(
                            names(ca), '__',
                            paste0(ca, collapse='_')
                        )
                        results[[outc]][[const$name]][[fw]][[feature_names]] = list()
                        for(ds in donor_sample){
                            results[[outc]][[const$name]][[fw]][[feature_names]][[ds]] = list()
                                                        if(ds=='all'){
                                dataset2 = copy(dataset)
                            } else if(ds=='most_similar'){
                                dataset2 = most_similar(dataset,
                                                        outc,
                                                        col_name_unit_name,
                                                        name_treated_unit,
                                                        treated_period,
                                                        min_period, 
                                                        col_name_period)
                            }

                            for(ny in num_pre_period_years){
                                if(is.na(ny)){
                                    l_name = paste0('n_pp_years_', treated_period-min_period)
                                    results[[outc]][[const$name]][[fw]][[feature_names]][[ds]][[l_name]] = list()
                                } else{
                                    l_name = paste0('n_pp_years_', ny)
                                    min_period = treated_period - ny 
                                    results[[outc]][[const$name]][[fw]][[feature_names]][[ds]][[l_name]] = list()   
                                }
                            


                        cat(
                            outc,
                            const$name,
                            fw,
                            feature_names,
                            ds,
                            l_name,
                            '\n'
                        )
                        
                        sc.pred = estimate_sc(
                            dataset2,
                            outc,
                            ca,
                            col_name_unit_name,
                            name_treated_unit,
                            col_name_period,
                            treated_period,
                            min_period,
                            end_period,
                            outcome_models = outcome_models,
                            feature_weights = fw,
                            w.constr = const                    
                            )
                        results[[outc]][[const$name]][[fw]][[feature_names]][[ds]][[l_name]]$estimate = sc.pred 
                        
                        sc.infer = inference_sc(
                            sc.pred,
                            dataset2,
                            inference_type='placebo',
                                P = NULL,
                                u.missp      = TRUE,
                                u.sigma      = "HC1",
                                u.order      = 1,
                                u.lags       = 0,
                                u.design     = NULL,
                                u.alpha      = 0.05,
                                e.method     = "all",
                                e.order      = 1,
                                e.lags       = 0,
                                e.design     = NULL,
                                e.alpha      = 0.05,
                                sims         = 200,
                                rho          = NULL,
                                rho.max      = 0.2,
                                lgapp        = "generalized",
                                cores        = 1,
                                w.bounds     = NULL,
                                e.bounds     = NULL,
                                verbose      = TRUE
                                )
                        results[[outc]][[const$name]][[fw]][[feature_names]][[ds]][[l_name]]$infer = sc.infer 
                    }}
                    }
                }
            }
        }
    
    return(results)
}