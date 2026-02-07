# Quick Start Example
#
# Minimal end-to-end example:
#   1) build scm data with scdata()
#   2) estimate with scest()
#   3) run placebo inference with estimate_sc() + inference_sc()

library(SCMs)

set.seed(123)
example_data <- expand.grid(
  unit = c("treated", "c1", "c2", "c3"),
  year = 2001:2010,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

example_data$outcome <-
  50 + 0.8 * (example_data$year - 2001) +
  ifelse(example_data$unit == "treated" & example_data$year >= 2007, 3, 0) +
  rnorm(nrow(example_data), 0, 0.5)

example_data$cov1 <-
  20 + 0.2 * (example_data$year - 2001) +
  rnorm(nrow(example_data), 0, 0.3)

scm_data <- scdata(
  df = example_data,
  id.var = "unit",
  time.var = "year",
  outcome.var = "outcome",
  period.pre = 2001:2006,
  period.post = 2007:2010,
  unit.tr = "treated",
  unit.co = c("c1", "c2", "c3"),
  covagg = list(
    list(var = "outcome_var", partition_periods = list(type = "by_period")),
    list(var = "cov1", compute = "mean")
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
  dataset = example_data,
  outcome = "outcome",
  covagg = list(
    list(var = "outcome_var", partition_periods = list(type = "by_period")),
    list(var = "cov1", compute = "mean")
  ),
  col_name_unit_name = "unit",
  name_treated_unit = "treated",
  col_name_period = "year",
  treated_period = 2007,
  min_period = 2001,
  end_period = 2010,
  outcome_models = "none",
  feature_weights = "optimize",
  w.constr = list(name = "simplex")
)

placebo <- inference_sc(
  sc_result,
  dataset = example_data,
  verbose = FALSE
)

print(placebo$abadie_significance)
