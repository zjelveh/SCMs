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
#'                     outcome.var = "outcome",
#'                     period.pre = 1990:1999, period.post = 2000:2005,
#'                     unit.tr = "treated_state",
#'                     unit.co = c("control1", "control2"))
#'                     
#' # Estimate synthetic control
#' scm_result <- scest(data_prep, w.constr = list(name = "simplex", Q = 1))
#' 
#' # Conduct inference
#' inference_result <- inference_sc(scm_result, df)
#' 
#' # Plot results
#' scplot(scm_result)
#' }
#'
#' @section Specification Curve Analysis:
#' For robustness analysis across multiple specifications:
#' \preformatted{
#' spec_results <- spec_curve(
#'   dataset = df,
#'   outcomes = "outcome",
#'   col_name_unit_name = "state",
#'   name_treated_unit = "treated_state",
#'   covagg = list(
#'     "Outcome Only" = list(
#'       label = "Outcome Only",
#'       operations = list(
#'         list(var = "outcome_var", partition_periods = list(type = "by_period"))
#'       )
#'     )
#'   ),
#'   treated_period = 2000,
#'   min_period = 1990,
#'   end_period = 2005,
#'   col_name_period = "year",
#'   constraints = list(
#'     list(name = "simplex"),
#'     list(name = "lasso", Q = 0.1)
#'   ),
#'   cores = 4
#' )
#' plot_spec_curve(spec_results)
#' }
#'
#' @section Advanced Features:
#' \itemize{
#'   \item Multiple constraint types: simplex, lasso, ridge, OLS
#'   \item Various outcome models: OLS, Ridge, Lasso, Augmented Synthetic Control
#'   \item Parallel processing for specification curve analysis
#'   \item Placebo-based inference for statistical significance
#'   \item Bootstrap null inference for additional uncertainty assessment
#' }
#'
NULL
