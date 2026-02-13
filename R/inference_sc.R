#' @title Perform Placebo Inference for Synthetic Control Results
#' @description Wrapper around \code{\link{inference_placebo}} for rank-based
#' placebo inference.
#'
#' @param sc.pred List. Synthetic-control estimation result, typically from
#'   \code{\link{estimate_sc}}.
#' @param dataset Data frame. Original panel dataset.
#' @param cores Integer. Number of cores for placebo runs. Default is 1.
#' @param verbose Logical. Whether to print progress information. Default is
#'   TRUE.
#' @param expected_direction Character. Expected direction of treatment effect:
#'   \code{"negative"}, \code{"positive"}, or \code{"two_sided"}. Default is
#'   \code{"negative"}.
#'
#' @return A list with \code{data}, \code{est.results}, \code{inference.results},
#' \code{rmse}, and \code{abadie_significance}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' toy <- data.frame(
#'   unit = rep(c("treated", "donor1", "donor2"), each = 6),
#'   year = rep(1:6, 3),
#'   y = c(10, 11, 12, 14, 15, 16,
#'         9, 10, 11, 12, 13, 14,
#'         11, 12, 13, 13, 14, 15)
#' )
#' sc_fit <- estimate_sc(
#'   dataset = toy,
#'   outcome = "y",
#'   covagg = list(
#'     list(var = "outcome_var", partition_periods = list(type = "by_period"))
#'   ),
#'   col_name_unit_name = "unit",
#'   name_treated_unit = "treated",
#'   col_name_period = "year",
#'   treated_period = 5,
#'   min_period = 1,
#'   end_period = 6,
#'   outcome_models = "none"
#' )
#' inf <- inference_sc(sc_fit, toy, cores = 1, verbose = FALSE)
#' names(inf)
#' }
inference_sc <- function(
    sc.pred,
    dataset,
    cores        = 1,
    verbose      = TRUE,
    expected_direction = "negative"
){
  if (is.null(sc.pred) || is.null(sc.pred$data) || is.null(sc.pred$est.results)) {
    stop("sc.pred must contain 'data' and 'est.results'")
  }

  if (is.null(sc.pred$data$specs) || is.null(sc.pred$data$specs$col.name.unit)) {
    stop("sc.pred$data$specs must include col.name.unit metadata")
  }

  if (!sc.pred$data$specs$col.name.unit %in% names(dataset)) {
    stop("Unit-id column '", sc.pred$data$specs$col.name.unit, "' not found in dataset")
  }

  required_fields <- c(
    "outcome", "covagg", "col_name_unit_name", "name_treated_unit",
    "col_name_period", "treated_period", "min_period", "end_period",
    "feature_weights", "outcome_models"
  )
  missing_fields <- required_fields[vapply(required_fields, function(x) is.null(sc.pred[[x]]), logical(1))]
  if (length(missing_fields) > 0) {
    stop(
      "sc.pred is missing required metadata for placebo inference: ",
      paste(missing_fields, collapse = ", "),
      ". Pass the object returned by estimate_sc()."
    )
  }

  if (is.null(sc.pred$est.results$outcome_model)) {
    stop("sc.pred$est.results$outcome_model is missing. Pass the object returned by estimate_sc().")
  }

  placebo_results <- inference_placebo(
    sc.pred = sc.pred,
    dataset = dataset,
    cores = cores,
    verbose = verbose,
    expected_direction = expected_direction
  )

  result <- list(
    data = sc.pred$data,
    est.results = sc.pred$est.results,
    inference.results = placebo_results$taus,
    rmse = placebo_results$rmse,
    abadie_significance = placebo_results$abadie_significance
  )

  class(result) <- "scm_inference"
  return(result)
}
