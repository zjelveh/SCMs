#' @title Estimate Synthetic Control Model
#' @description This function processes and estimates a Synthetic Control Model (SCM) based on the provided dataset and parameters.
#'
#' @param dataset A data frame containing the panel data.
#' @param outcome Character. Name of the outcome variable (column in dataset).
#' @param covagg List of covariates used for matching (columns in dataset).
#' @param col_name_unit_name Character. Column name in dataset containing unit names (e.g., state names).
#' @param name_treated_unit Character. Name of the treated unit (e.g., the state that received treatment).
#' @param col_name_period Character. Column name in dataset containing time periods.
#' @param treated_period Numeric. Time period when treatment starts for the treated unit.
#' @param min_period Numeric. Earliest time period in the dataset.
#' @param end_period Numeric. Latest time period in the dataset.
#' @param outcome_models Character vector. Outcome models to fit. Default is c("None", "OLS", "Ridge", "Lasso", "AugSynth").
#' @param feature_weights Character vector. Method for assigning weights to predictors. Default is c("uniform", "optimized").
#' @param w.constr Optional. Constraints on the weights.
#' @param V Character. Covariance matrix estimation method. Default is "separate".
#' @param V.mat Optional. Pre-computed covariance matrix.
#' @param solver Character. Solver to use for the optimization problem. Default is "ECOS".
#' @param P Optional. Number of factors to use (for AugSynth).
#'
#' @return A list of class "scpi" containing the estimated SCM results and input parameters.
#'
#' @export
#'
#' @examples
#' # Example usage (replace with actual example when available)
#' # result <- estimate_sc(dataset = my_data, outcome = "gdp", covagg = c("population", "unemployment"),
#' #                       col_name_unit_name = "state", name_treated_unit = "California",
#' #                       col_name_period = "year", treated_period = 2000, min_period = 1990, end_period = 2010)
estimate_sc <- function(dataset,           # Input dataset (panel data format)
                        outcome,           # Name of the outcome variable (column in dataset)
                        covagg,            # List of covariates used for matching (columns in dataset)
                        col_name_unit_name, # Column name in dataset containing unit names (e.g., state names)
                        name_treated_unit, # Name of the treated unit (e.g., the state that received treatment)
                        col_name_period,    # Column name in dataset containing time periods
                        treated_period,    # Time period when treatment starts for the treated unit
                        min_period,        # Earliest time period in the dataset
                        end_period,        # Latest time period in the dataset
                        outcome_models = c("None", "OLS", "Ridge", "Lasso", "AugSynth"), # Outcome models to fit
                        feature_weights = c("uniform", "optimized"),  # Method for assigning weights to predictors
                        w.constr = NULL,    # Constraints on the weights (optional)
                        V = "separate",     # Covariance matrix estimation method (usually "separate")
                        V.mat = NULL,       # Pre-computed covariance matrix (optional)
                        solver = "ECOS"){#,    # Solver to use for the optimization problem
                        #P = NULL) {         # Number of factors to use (optional, for AugSynth)


  data <- create_scm_dataset(       # Create a formatted dataset for SCM
    dataset = dataset, 
    outcome = outcome, 
    covagg = covagg,
    col_name_unit_name = col_name_unit_name,
    name_treated_unit = name_treated_unit,
    col_name_period = col_name_period,
    treated_period = treated_period,
    min_period = min_period,
    end_period = end_period
  )

  # Error handling: Check if data preparation was successful
  if (!methods::is(data, "scdata") & !methods::is(data, "scdataMulti")) {
    stop("data should be the object returned by running scdata or scdata_multi!")
  }

  # 2. SCM Estimation
  scm_model <- scest(data = data,     # Estimate the Synthetic Control
    w.constr = w.constr,
    feature_weights = feature_weights,
    V = V,
    V.mat = V.mat,
    solver = solver
  )

  # 3. Outcome Model Estimation
  sc.pred <- run_outcome_models( # Run various outcome models on the SCM results
    scm_model = scm_model,
    scm_data = data,
    treated_unit = data$specs$treated.units,
    outcome_models = outcome_models,
    period_post = data$specs$period.post,
    Y = data$Y.donors.post,
    Z0 = data$Y.donors,
    Z1 = data$Y.pre
  )

  # 4. Prepare and return results
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

  # Set the appropriate class type
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