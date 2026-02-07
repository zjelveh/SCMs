# Generate the README spec-curve figure used in README.md
# Run from package root:
#   Rscript inst/examples/generate_readme_spec_curve_image.R

if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(SCMs)
}
library(ggplot2)

set.seed(42)
d <- expand.grid(
  unit = c("treated", "c1", "c2", "c3", "c4"),
  year = 2001:2010,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

d$unit_id <- gsub("[^A-Za-z0-9_]", "_", d$unit)

d$outcome <-
  100 + 0.5 * (d$year - 2001) +
  ifelse(d$unit == "treated" & d$year >= 2007, 2, 0) +
  rnorm(nrow(d), 0, 0.3)

d$population <-
  50 + 0.1 * (d$year - 2001) +
  rnorm(nrow(d), 0, 0.2)

spec_results <- spec_curve(
  dataset = d,
  outcomes = "outcome",
  col_name_unit_name = "unit_id",
  name_treated_unit = "treated",
  covagg = list(
    "Outcome Per Period" = list(
      label = "Outcome Per Period",
      operations = list(
        list(var = "outcome_var", partition_periods = list(type = "by_period"))
      )
    ),
    "Outcome Per Period + Population Mean" = list(
      label = "Outcome Per Period + Population Mean",
      operations = list(
        list(var = "outcome_var", partition_periods = list(type = "by_period")),
        list(var = "population", compute = "mean")
      )
    )
  ),
  treated_period = 2007,
  min_period = 2001,
  end_period = 2010,
  col_name_period = "year",
  feature_weights = c("uniform", "optimize"),
  outcome_models = c("none", "ridge"),
  donor_sample = "all",
  constraints = list(list(name = "simplex")),
  constants = FALSE,
  cores = 1,
  verbose = FALSE,
  inference_type = "placebo"
)

plot_obj <- plot_spec_curve(
  long_data = spec_results,
  name_treated_unit = "treated",
  outcomes = "outcome",
  show_shap = FALSE,
  show_pvalues = TRUE,
  test_statistic = "treatment_effect"
)

ggsave(
  filename = "inst/examples/readme_spec_curve.png",
  plot = plot_obj$final_plot,
  width = 9,
  height = 10,
  units = "in",
  dpi = 140
)

cat("Saved: inst/examples/readme_spec_curve.png\n")
