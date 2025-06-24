#' SCMS: Synthetic Control Methods with Specification Curve Analysis
#'
#' This package provides tools for implementing Synthetic Control Methods (SCM),
#' performing inference, conducting specification curve analysis, and visualizing results.
#'
#' @keywords internal
"_PACKAGE"


#' @section Important Functions:
#' \itemize{
#'   \item \code{scdata}: Prepare data for synthetic control analysis
#'   \item \code{scest}: Estimate the Synthetic Control Model
#'   \item \code{inference_sc}: Conduct inference on SCM results (placebo-based)
#'   \item \code{inference_placebo}: Alternative placebo inference method
#'   \item \code{spec_curve}: Perform specification curve analysis
#'   \item \code{scplot}: Plot the results of the Synthetic Control Model
#'   \item \code{plot_spec_curve}: Visualize the specification curve
#'   \item \code{create_scm_dataset}: Create datasets for SCM analysis
#'   \item \code{estimate_sc}: High-level estimation wrapper
#' }
#'
#' @section Basic Workflow:
#' The typical workflow involves:
#' 1. Prepare your data with \code{scdata()}
#' 2. Estimate the synthetic control with \code{scest()}
#' 3. Conduct inference with \code{inference_sc()} or \code{inference_placebo()}
#' 4. Visualize results with \code{scplot()}
#' 
#' Example:
#' \preformatted{
#' # Prepare data
#' data_prep <- scdata(df, id.var = "state", time.var = "year", 
#'                     outcome.var = "outcome", treatment.identifier = "treated_state")
#'                     
#' # Estimate synthetic control
#' scm_result <- scest(data_prep, w.constr = list(name = "simplex", Q = 1))
#' 
#' # Conduct inference
#' inference_result <- inference_sc(data_prep, scm_result)
#' 
#' # Plot results
#' scplot(scm_result)
#' }
#'
#' @section Specification Curve Analysis:
#' For robustness analysis across multiple specifications:
#' \preformatted{
#' spec_results <- spec_curve(data_prep, 
#'                           w.constr.options = list(simplex = list(name = "simplex", Q = 1),
#'                                                  lasso = list(name = "lasso", Q = 0.1)),
#'                           cores = 4)
#' plot_spec_curve(spec_results)
#' }
#'
#' @section Advanced Features:
#' \itemize{
#'   \item Multiple constraint types: simplex, lasso, ridge, OLS
#'   \item Various outcome models: OLS, Ridge, Lasso, Augmented Synthetic Control
#'   \item Parallel processing for specification curve analysis
#'   \item SCPI-based inference for prediction intervals
#'   \item Placebo-based inference for statistical significance
#' }
#'
NULL
