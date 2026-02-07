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

placebo <- inference_sc(
  fit,
  dataset = dataset,
  inference_type = "placebo",
  verbose = FALSE
)

print(placebo$abadie_significance)
