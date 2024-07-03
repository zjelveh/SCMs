#' SCMS: Synthetic Control Methods with Specification Curve Analysis
#'
#' This package provides tools for implementing Synthetic Control Methods (SCM),
#' performing inference, conducting specification curve analysis, and visualizing results.
#'
#' @keywords internal
"_PACKAGE"


#' @importFrom stats predict fitted coef
#' @importFrom kernlab alpha
#' @importFrom data.table data.table as.data.table

#' @importFrom CVXR Variable Minimize Problem solve
#' @rawNamespace import(stats, except = c(lag, filter))
#' @rawNamespace import(dplyr, except = c(lag, filter, first))
#' @rawNamespace import(CVXR, except = c(id, size))
#' @rawNamespace import(kernlab, except = c(alpha, predict, fitted, coef))
#'
#' @section Important Functions:
#' \itemize{
#'   \item \code{estimate_scm}: Estimates the Synthetic Control Model
#'   \item \code{perform_inference}: Conducts inference on SCM results
#'   \item \code{run_spec_curve}: Performs specification curve analysis
#'   \item \code{plot_scm}: Plots the results of the Synthetic Control Model
#'   \item \code{plot_spec_curve}: Visualizes the specification curve
#' }
#'
#' @section Estimating SCM:
#' Use \code{estimate_scm()} to fit a Synthetic Control Model to your data.
#' Example:
#' \preformatted{
#' scm_result <- estimate_scm(data, treatment_unit, control_units, time_variable, outcome)
#' }
#'
#' @section Inference:
#' Use \code{perform_inference()} to conduct statistical inference on your SCM results.
#' Example:
#' \preformatted{
#' inference_result <- perform_inference(scm_result, method = "placebo")
#' }
#'
#' @section Specification Curve:
#' Use \code{run_spec_curve()} to perform specification curve analysis.
#' Example:
#' \preformatted{
#' spec_curve <- run_spec_curve(data, model_variations)
#' }
#'
#' @section Plotting:
#' Use \code{plot_scm()} and \code{plot_spec_curve()} to visualize your results.
#' Examples:
#' \preformatted{
#' plot_scm(scm_result)
#' plot_spec_curve(spec_curve)
#' }
#'
NULL
