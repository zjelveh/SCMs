# Generate README spec-curve figures from the Basque paper specification grid.
# Recomputes the full grid (320 specifications with placebo inference) on all but one core.
# Run from package root:
#   Rscript inst/examples/generate_readme_spec_curve_image.R
#
# Outputs:
#   - inst/examples/readme_spec_curve.png (README inline image; white background)
#   - inst/examples/readme_spec_curve.svg (README inline vector image)
#   - inst/examples/readme_spec_curve.pdf (vector companion)

if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(SCMs)
}
library(ggplot2)
library(data.table)

# Print warnings as they occur rather than a deferred "There were 50 or more warnings".
options(warn = 1)
n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

treated_unit <- "Basque_Country_Pais_Vasco"
outcome <- "gdpcap"
spec_width <- 10
spec_height <- 9
png_out <- "inst/examples/readme_spec_curve.png"
svg_out <- "inst/examples/readme_spec_curve.svg"
pdf_out <- "inst/examples/readme_spec_curve.pdf"

build_basque_covagg <- function() {
  all_vars <- c(
    "gdpcap", "school.illit", "school.prim", "school.med", "school.high",
    "school.post.high", "invest", "sec.agriculture", "sec.energy",
    "sec.industry", "sec.construction", "sec.services.venta",
    "sec.services.nonventa", "popdens"
  )

  make_per_period_ops <- function(vars) {
    lapply(vars, function(v) {
      list(var = v, partition_periods = list(type = "by_period"))
    })
  }

  make_mean_ops <- function(vars) {
    lapply(vars, function(v) {
      list(var = v, compute = "mean")
    })
  }

  list(
    "Outcome Per-Period" = list(
      label = "Outcome Per-Period",
      operations = list(
        list(var = "outcome_var", partition_periods = list(type = "by_period"))
      )
    ),
    "All Meaned" = list(
      label = "All Meaned",
      operations = make_mean_ops(all_vars)
    ),
    "Outcome Per-Period + All Meaned" = list(
      label = "Outcome Per-Period + All Meaned",
      operations = c(
        list(list(var = "outcome_var", partition_periods = list(type = "by_period"))),
        make_mean_ops(all_vars)
      )
    ),
    "All Per-Period" = list(
      label = "All Per-Period",
      operations = make_per_period_ops(all_vars)
    )
  )
}

basque <- fread("inst/extdata/basque.csv")
basque[, regionname := gsub("\\(|\\)", "", regionname)]
basque[, regionname := gsub(" ", "_", regionname)]
basque <- basque[regionno > 1]

params <- list(
  outcomes = outcome,
  col_name_unit_name = "regionname",
  name_treated_unit = treated_unit,
  treated_period = 1970,
  min_period = 1955,
  end_period = 1997,
  col_name_period = "year",
  inference_type = "placebo",
  feature_weights = c("uniform", "optimize"),
  donor_sample = c("all"),
  constant = c(FALSE, TRUE),
  outcome_models = c("none", "augsynth", "lasso", "ridge", "ols"),
  constraints = list(
    list(name = "pensynth"),
    list(name = "simplex"),
    list(name = "lasso"),
    list(name = "ridge")
  ),
  covagg = build_basque_covagg(),
  inference_config = list(bootstrap_n_replications = 100, verbose = FALSE)
)

results <- run_spec_curve_analysis(
  dataset = basque,
  params = params,
  cores = n_cores,
  output_dir = NULL
)

xgb_config <- create_xgboost_config(
  dataset_name = "basque_readme",
  treated_unit_name = treated_unit,
  outcome_filter = outcome,
  spec_features = c("feat", "outcome_model", "const", "fw", "constant"),
  treated_unit_only = TRUE
)
shap_result <- run_xgboost_shap_analysis(results$results, xgb_config, compute_loo = TRUE)

plot_obj <- plot_spec_curve(
  long_data = results,
  name_treated_unit = treated_unit,
  outcomes = outcome,
  show_shap = TRUE,
  shap_label_type = "signed",
  shap_values = shap_result$shapley,
  show_predictions = TRUE,
  predictions = shap_result$predictions,
  show_pvalues = TRUE,
  null_distribution = "placebo",
  test_statistic = "treatment_effect"
)

final_plot <- plot_obj$final_plot

ggsave(
  filename = png_out,
  plot = final_plot,
  width = spec_width,
  height = spec_height,
  units = "in",
  dpi = 150,
  bg = "white"
)

ggsave(
  filename = svg_out,
  plot = final_plot,
  width = spec_width,
  height = spec_height,
  units = "in",
  # Cairo SVG embeds glyph outlines, so the rich-text labels render exactly;
  # svglite leaves text to the browser's fonts and misplaces them.
  device = grDevices::svg,
  bg = "white"
)

ggsave(
  filename = pdf_out,
  plot = final_plot,
  width = spec_width,
  height = spec_height,
  units = "in",
  device = "pdf",
  bg = "white"
)

cat("Saved:", png_out, "\n")
cat("Saved:", svg_out, "\n")
cat("Saved:", pdf_out, "\n")
