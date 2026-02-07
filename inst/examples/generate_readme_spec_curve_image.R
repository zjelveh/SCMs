# Generate the README spec-curve figure used in README.md
# Run from package root:
#   Rscript inst/examples/generate_readme_spec_curve_image.R

if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(SCMs)
}
library(ggplot2)

basque <- data.table::fread("inst/extdata/basque.csv")
clean_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_{2,}", "_", x)
  x <- gsub("_+$", "", x)
  x
}
basque$region_id <- clean_id(basque$regionname)

treated_unit <- clean_id("Basque Country (Pais Vasco)")

spec_results <- spec_curve(
  dataset = basque,
  outcomes = "gdpcap",
  col_name_unit_name = "region_id",
  name_treated_unit = treated_unit,
  covagg = list(
    "Outcome Per Period" = list(
      label = "Outcome Per Period",
      operations = list(
        list(var = "outcome_var", partition_periods = list(type = "by_period"))
      )
    ),
    "Outcome Mean" = list(
      label = "Outcome Mean",
      operations = list(
        list(
          var = "outcome_var",
          partition_periods = list(type = "all"),
          compute = "mean"
        )
      )
    )
  ),
  treated_period = 1970,
  min_period = 1955,
  end_period = 1997,
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
  name_treated_unit = treated_unit,
  outcomes = "gdpcap",
  show_shap = TRUE,
  show_predictions = TRUE,
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
