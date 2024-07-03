#' @title Generate Specification Curve for Synthetic Control Method
#' @description This function generates a specification curve by running multiple variations of Synthetic Control Method (SCM) estimations across different specifications.
#'
#' @param dataset A data frame containing the panel data.
#' @param outcomes Character vector. Names of the outcome variables to analyze.
#' @param col_name_unit_name Character. Column name in dataset containing unit names.
#' @param name_treated_unit Character. Name of the treated unit.
#' @param covagg List of covariates used for matching.
#' @param treated_period Numeric. Time period when treatment starts for the treated unit.
#' @param min_period Numeric. Earliest time period in the dataset.
#' @param end_period Numeric. Latest time period in the dataset.
#' @param col_name_period Character. Column name in dataset containing time periods.
#' @param feature_weights Character vector. Methods for assigning weights to predictors. Default is c('uniform', 'optimize').
#' @param num_pre_period_years Numeric or NA. Number of pre-treatment periods to use. If NA, uses all available pre-treatment periods.
#' @param outcome_models Character vector. Outcome models to fit. Default is c('none', 'augsynth', 'ridge', 'lasso', 'ols').
#' @param donor_sample Character vector. Method for selecting donor units. Default is c('all', 'most_similar').
#' @param sim_function Function. Function to select most similar donor units. Default is most_similar.
#' @param constraints List of lists. Each inner list specifies a constraint type for the SCM estimation.
#'
#' @return A nested list containing results for each combination of specifications.
#'
#' @export
#'
#' @examples
#' # Example usage (replace with actual example when available)
#' # results <- spec_curve(dataset = my_data, outcomes = c("gdp", "unemployment"),
#' #                       col_name_unit_name = "state", name_treated_unit = "California",
#' #                       covagg = list(c("population", "education")), treated_period = 2000,
#' #                       min_period = 1990, end_period = 2010, col_name_period = "year")
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
  # Initialize results list
  results = list()
  
  # Iterate over all combinations of specifications
  for(outc in outcomes){
    results[[outc]] = list()
    for(const in constraints){
      results[[outc]][[const$name]] = list()
      for(fw in feature_weights){
        results[[outc]][[const$name]][[fw]] = list()
        for(ca in covagg){ 
          # Create feature names
          feature_names = paste0(
            names(ca), '__',
            paste0(ca, collapse='_')
          )
          results[[outc]][[const$name]][[fw]][[feature_names]] = list()
          for(ds in donor_sample){
            results[[outc]][[const$name]][[fw]][[feature_names]][[ds]] = list()
            # Select donor sample
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
              # Set pre-period years
              if(is.na(ny)){
                l_name = paste0('n_pp_years_', treated_period-min_period)
                results[[outc]][[const$name]][[fw]][[feature_names]][[ds]][[l_name]] = list()
              } else{
                l_name = paste0('n_pp_years_', ny)
                min_period = treated_period - ny 
                results[[outc]][[const$name]][[fw]][[feature_names]][[ds]][[l_name]] = list()   
              }
              
              # Print current specification
              cat(
                outc,
                const$name,
                fw,
                feature_names,
                ds,
                l_name,
                '\n'
              )
              
              # Estimate Synthetic Control
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
              
              # Perform inference
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
            }
          }
        }
      }
    }
  }
  
  return(results)
}