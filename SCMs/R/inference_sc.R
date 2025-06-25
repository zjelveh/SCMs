#' @title Perform Inference for Synthetic Control Method
#' @description This function performs inference for synthetic control estimates using either SCPI (Synthetic Control Prediction Intervals) or placebo methods.
#'
#' @param sc.pred List. Results from the synthetic control estimation.
#' @param dataset Data frame. The original dataset used for synthetic control estimation.
#' @param inference_type Character. Type of inference to perform: "scpi" or "placebo". Default is "scpi".
#' @param P Matrix. Prediction matrix (not used in the current implementation).
#' @param u.missp Logical. Whether to account for model misspecification. Default is TRUE.
#' @param u.sigma Character. Type of heteroskedasticity-consistent standard errors. Default is "HC1".
#' @param u.order Numeric. Order of the misspecification. Default is 1.
#' @param u.lags Numeric. Number of lags for misspecification. Default is 0.
#' @param u.design Matrix. Design matrix for misspecification (optional).
#' @param u.alpha Numeric. Significance level for uncertainty intervals. Default is 0.05.
#' @param e.method Character. Method for estimating errors. Default is "all".
#' @param e.order Numeric. Order of the error process. Default is 1.
#' @param e.lags Numeric. Number of lags for error process. Default is 0.
#' @param e.design Matrix. Design matrix for error process (optional).
#' @param e.alpha Numeric. Significance level for error intervals. Default is 0.05.
#' @param sims Numeric. Number of simulations for inference. Default is 200.
#' @param rho Numeric. Correlation parameter (optional).
#' @param rho.max Numeric. Maximum correlation. Default is 0.2.
#' @param lgapp Character. Type of lag polynomial approximation. Default is "generalized".
#' @param cores Numeric. Number of cores to use for parallel processing. Default is 1.
#' @param w.bounds List. Bounds for weights (optional).
#' @param e.bounds List. Bounds for errors (optional).
#' @param verbose Logical. Whether to print progress information. Default is TRUE.
#'
#' @return A list containing the inference results, either from SCPI or placebo methods.
#'
#' @export
#'
#' @examples
#' # Example usage (replace with actual example when available)
#' # inference_results <- inference_sc(sc_results, my_data, inference_type = "scpi")
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
  
  # Check if sc.pred structure is valid
  if (is.null(sc.pred)) {
    return(NULL)
  }
  
  if (is.null(sc.pred$data) || is.null(sc.pred$data$specs)) {
    return(NULL)
  }
  
  # Check column name
  col_name_unit <- sc.pred$data$specs$col.name.unit
  
  if (is.null(col_name_unit)) {
    # Try to infer the unit column name from common patterns
    possible_names <- c("country", "state", "unit", "id", "unit_id", "unit_name", "ori9", "stateid")
    col_name_unit <- NULL
    for (name in possible_names) {
      if (name %in% names(dataset)) {
        col_name_unit <- name
        break
      }
    }
    
    if (is.null(col_name_unit)) {
      return(NULL)
    }
  }
  
  if (!col_name_unit %in% names(dataset)) {
    return(NULL)
  }
  
  # Start by defining class.type based on sc.pred
  if (!is.null(sc.pred$data$specs) && !is.null(sc.pred$data$specs$class.type)) {
    class.type <- sc.pred$data$specs$class.type
  } else {
    # Default class type if not found in data
    class.type <- "scpi_data"
  }
  
  if (inference_type == "scpi") {
    # Perform SCPI inference
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
    
    # Prepare result structure for SCPI
    result <- list(
      data = sc.pred$data,
      est.results = sc.pred$est.results,
      inference.results = inference.results
    )
    
    class(result) <- "scpi"
    
    # Set class type based on data structure
    if (class.type == "scpi_data") {
      result$data$specs$class.type <- "scpi_scpi"
    } else if (class.type == "scpi_data_multi") {
      result$data$specs$class.type <- "scpi_scpi_multi"
    }
  } else {
    # Perform placebo inference
    inference.results <- tryCatch({
      inference_placebo(
        sc.pred = sc.pred,
        dataset = dataset, 
        cores = cores,
        verbose = verbose
      )
    }, error = function(e) {
      return(NULL)
    })
    
    # Now inference.results contains both taus and rmse
    result <- list(
      data = sc.pred$data,
      est.results = sc.pred$est.results,
      inference.results = inference.results$taus,
      rmse = inference.results$rmse
    )
  }
  
  return(result)
}