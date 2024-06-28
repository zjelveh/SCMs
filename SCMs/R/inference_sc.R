inference_sc <- function(
    sc.pred,
    dataset,
    inference_type='scpi',
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
){
      
  if (inference_type == "scpi") {
    inference.results <- inference_scpi(
      sc.pred = sc.pred,
      u.missp      = u.missp,
      u.sigma      = u.sigma,
      u.order      = u.order,
      u.lags       = u.lags,
      u.design     = u.design,
      u.alpha      = u.alpha,
      e.method = e.method,
      e.order = e.order,
      e.lags = e.lags,
      e.design = e.design,
      e.alpha = e.alpha,
      rho = rho,
      rho.max = rho.max,
      lgapp = lgapp, 
      w.bounds = w.bounds,
      e.bounds = e.bounds,
      sims = sims,
      cores = cores,
      verbose = verbose
    )
    result <- list(
    data = sc.pred$data,
    est.results = sc.pred$est.results,
    inference.results = inference.results
  )
  
  class(result) <- "scpi"
  if (class.type == "scpi_data") {
    result$data$specs$class.type <- "scpi_scpi"
  } else if (class.type == "scpi_data_multi") {
    result$data$specs$class.type <- "scpi_scpi_multi"
  }
  } else{
    inference.results <- inference_placebo(
      sc.pred = sc.pred,
      dataset = dataset, 
      cores = cores,
      verbose = verbose
    )
    result = inference.results
  }
  
  
  
  return(result)
  
}

