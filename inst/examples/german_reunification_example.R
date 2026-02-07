# German Reunification Example
#
# Dataset-backed SCM workflow on included `scpi_germany.csv`.

library(SCMs)
library(data.table)

dataset <- fread(system.file("extdata/scpi_germany.csv", package = "SCMs"))

# Keep rows with observed outcome and use sanitized unit ids for robust matching.
dataset <- dataset[!is.na(gdp)]
dataset[, unit_id := gsub("[^A-Za-z0-9_]", "_", country)]

treated_unit <- "West_Germany"
control_units <- setdiff(unique(dataset$unit_id), treated_unit)

scm_data <- scdata(
  df = dataset,
  id.var = "unit_id",
  time.var = "year",
  outcome.var = "gdp",
  period.pre = 1960:1990,
  period.post = 1991:2003,
  unit.tr = treated_unit,
  unit.co = control_units,
  covagg = list(
    list(var = "outcome_var", partition_periods = list(type = "by_period"))
  )
)

fit <- scest(
  data = scm_data,
  w.constr = list(name = "simplex"),
  feature_weights = "optimize"
)

cat("Estimated donor weights:\n")
print(fit$est.results$w)
scplot(fit)

sc_result <- estimate_sc(
  dataset = dataset,
  outcome = "gdp",
  covagg = list(
    list(var = "outcome_var", partition_periods = list(type = "by_period"))
  ),
  col_name_unit_name = "unit_id",
  name_treated_unit = treated_unit,
  col_name_period = "year",
  treated_period = 1991,
  min_period = 1960,
  end_period = 2003,
  outcome_models = "none",
  feature_weights = "optimize",
  w.constr = list(name = "simplex")
)

placebo <- inference_sc(
  sc_result,
  dataset = dataset,
  verbose = FALSE
)

print(placebo$abadie_significance)
