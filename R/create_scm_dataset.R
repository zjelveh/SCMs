#' @title Create Dataset for Synthetic Control Method
#' @description This function prepares a dataset for use with the Synthetic Control Method (SCM) by formatting and standardizing the input data.
#'
#' @param dataset Data frame or data.table containing the panel data.
#' @param outcome Character. Name of the outcome variable column.
#' @param covagg List of covariates used for matching.
#' @param col_name_unit_name Character. Name of the column containing unit names.
#' @param name_treated_unit Character. Name of the treated unit.
#' @param col_name_period Character. Name of the column containing time periods.
#' @param treated_period Numeric. Time period when treatment starts.
#' @param min_period Numeric or NULL. Earliest time period to include. If NULL, uses the minimum in the dataset.
#' @param end_period Numeric or NULL. Latest time period to include. If NULL, uses the maximum in the dataset.
#' @param constant Logical. Whether to include a constant term. Default is FALSE.
#'
#' @return A list containing the prepared data for SCM analysis.
#'
#'
#' @export
#'
#' @examples
#' # Example usage (replace with actual example when available)
#' # scm_data <- create_scm_dataset(dataset = my_data, outcome = "gdp", 
#' #                                covagg = list(
#' #                                  population = list(var = "population"),
#' #                                  education = list(var = "education")
#' #                                ),
#' #                                col_name_unit_name = "state", name_treated_unit = "California",
#' #                                col_name_period = "year", treated_period = 2000,
#' #                                min_period = 1990, end_period = 2010)
create_scm_dataset <- function(dataset,
                               outcome,
                               covagg,
                               col_name_unit_name,
                               name_treated_unit,
                               col_name_period,
                               treated_period,
                               min_period,
                               end_period,
                               constant = FALSE) {

  # Convert input to data.table and create a proper copy to avoid shallow copy warnings
  dataset <- data.table::copy(data.table::as.data.table(dataset))

  # Store original unit names before any transformation
  dataset[, original_unit_name := get(col_name_unit_name)]
  
  # Apply consistent name transformation (same as scdata function)
  # This ensures all operations work with consistent naming
  name_treated_unit <- gsub(' +|\\.+|-+', '_', name_treated_unit)
  
  # Transform unit names in the dataset using proper data.table syntax
  if(is.character(dataset[[col_name_unit_name]])){
    dataset[, (col_name_unit_name) := gsub(' +|\\.+|-+', '_', get(col_name_unit_name))]
  }

  # Create unit name and number columns (now with transformed names)
  dataset[, unit_name := get(col_name_unit_name)]
  dataset[, unit_numbers := as.numeric(as.factor(unit_name))]
  
  # Create treatment indicator (now works because names are consistent)
  dataset[, trt := ifelse(unit_name == name_treated_unit, 1, 0)]
  
  # Set maximum untreated period
  max_untreated_period <- treated_period - 1
  
  # Get unique control and treated units
  unit.co <- unique(dataset[trt == 0][[col_name_unit_name]])
  unit.tr <- unique(dataset[trt == 1][[col_name_unit_name]])
  
  # Set min and max periods if not provided
  if (is.null(min_period)) {
    min_period <- min(dataset[[col_name_period]])
  }
  if (is.null(end_period)) {
    end_period <- max(dataset[[col_name_period]])
  }
  if (is.null(max_untreated_period)) {
    max_untreated_period <- treated_period - 1
  }

  # Create SCM dataset using scpi package

  # Check which variables actually exist in the dataset
  dataset_vars <- names(dataset)
  
  # Check for variables with periods that might be missing
  period_vars <- dataset_vars[grepl("\\.", dataset_vars)]

  scpi_data <- scdata(
    df = as.data.frame(dataset),
    id.var = col_name_unit_name,
    outcome.var = outcome,
    time.var = col_name_period,
    period.pre = (min_period:max_untreated_period),
    period.post = (treated_period:end_period),
    unit.tr = unit.tr,
    unit.co = unit.co,
    constant = constant,
    covagg = covagg
  )

  # Calculate standard deviations and standardize data
  X0_sds <- apply(scpi_data$B, 1, sd)
  scpi_data$A_original <- copy(scpi_data$A)
  scpi_data$B_original <- copy(scpi_data$B)
  scpi_data$X0_sds <- X0_sds
  scpi_data$A <- scpi_data$A / X0_sds
  scpi_data$B <- scpi_data$B / X0_sds

  return(scpi_data)
}
