estimate_sc <- function(dataset,
                outcome,
                covagg,
                col_name_unit_name,
                name_treated_unit,
                col_name_period,
                treated_period,
                min_period,
                end_period,
                 outcome_models = c("None", "OLS", "Ridge", "Lasso", "AugSynth"),
                 feature_weights = c("uniform", "optimized"),
                 w.constr = NULL,
                 V = "separate",
                 V.mat = NULL,
                 solver = "ECOS",
                 P = NULL) {



  data = create_scm_dataset(
      dataset=dataset,
      outcome=outcome,
      covagg=covagg,
      col_name_unit_name=col_name_unit_name,
      name_treated_unit=name_treated_unit,
      col_name_period=col_name_period,
      treated_period=treated_period,
      min_period=min_period,
      end_period=end_period
  )


  if ((methods::is(data, "scdata") || methods::is(data, "scdataMulti")) == FALSE) {
    stop("data should be the object returned by running scdata or scdata_multi!")
  }

  scm_model <- scest(data = data,
    w.constr = w.constr,
    feature_weights=feature_weights,
    V = V,
    V.mat = V.mat,
    solver = solver)


  sc.pred <- run_outcome_models(
    scm_model = scm_model,
    scm_data = data,
    treated_unit = data$specs$treated.units,
    outcome_models = outcome_models,
    period_post = data$specs$period.post,
    Y = data$Y.donors.post,
    Z0 = data$Y.donors,
    Z1 = data$Y.pre
  )

  result <- list(
    data = data,
    est.results = sc.pred$est.results,
    outcome_models = outcome_models,
    feature_weights = feature_weights,
    covagg=covagg,
    w.constr = w.constr,
    V = V,
    V.mat = V.mat,
    solver = solver,
    outcome = outcome,
    col_name_unit_name = col_name_unit_name,
    name_treated_unit = name_treated_unit,
    col_name_period = col_name_period,
    treated_period = treated_period,
    min_period = min_period,
    end_period = end_period
  )

  if (methods::is(data, "scdata") == TRUE) {
    class.type <- "scpi_data"
  } else if (methods::is(data, "scdataMulti") == TRUE) {
    class.type <- "scpi_data_multi"
  }

  class(result) <- "scpi"
  if (class.type == "scpi_data") {
    result$data$specs$class.type <- "scpi_scpi"
  } else if (class.type == "scpi_data_multi") {
    result$data$specs$class.type <- "scpi_scpi_multi"
  }

  return(result)
}
