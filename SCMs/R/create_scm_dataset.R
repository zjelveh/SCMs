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
  require(data.table)
  dataset <- data.table::as.data.table(dataset)
  dataset[, unit_name := get(col_name_unit_name)]
  
  dataset[, unit_numbers := as.numeric(as.factor(unit_name))]
  dataset[, trt := ifelse(unit_name == name_treated_unit, 1, 0)]
  max_untreated_period <- treated_period - 1
  unit.co <- unique(dataset[trt == 0][[col_name_unit_name]])
  unit.tr <- unique(dataset[trt == 1][[col_name_unit_name]])


  if (is.null(min_period)) {
    min_period <- min(dataset[[period.var]])
  }
  if (is.null(end_period)) {
    end_period <- max(dataset[[period.var]])
  }
  if (is.null(max_untreated_period)) {
    max_untreated_period <- treat_period - 1
  }

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

  X0_sds <- apply(scpi_data$B, 1, sd)
  scpi_data$A_original <- copy(scpi_data$A)
  scpi_data$B_original <- copy(scpi_data$B)
  scpi_data$X0_sds <- X0_sds
  scpi_data$A <- scpi_data$A / X0_sds
  scpi_data$B <- scpi_data$B / X0_sds

  return(scpi_data)
}
