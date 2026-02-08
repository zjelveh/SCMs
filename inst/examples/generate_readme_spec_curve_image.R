# Generate README spec-curve figures from the Basque paper specification grid.
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

treated_unit <- "Basque_Country_Pais_Vasco"
outcome <- "gdpcap"
spec_width <- 10
spec_height <- 12
png_out <- "inst/examples/readme_spec_curve.png"
svg_out <- "inst/examples/readme_spec_curve.svg"
pdf_out <- "inst/examples/readme_spec_curve.pdf"

paper_results_path <- normalizePath("../../data/basque_results/basque_full_results.rdata", mustWork = FALSE)
paper_shap_path <- normalizePath("../../data/basque_results/basque_shap_result.rdata", mustWork = FALSE)

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

load_paper_cache <- function() {
  if (!file.exists(paper_results_path) || !file.exists(paper_shap_path)) {
    return(NULL)
  }

  load_env <- new.env(parent = emptyenv())
  load(paper_results_path, envir = load_env)
  load(paper_shap_path, envir = load_env)

  if (!exists("results", envir = load_env) || !exists("shap_result", envir = load_env)) {
    stop("Paper cache files found but expected objects ('results', 'shap_result') were not present.")
  }

  results_obj <- get("results", envir = load_env)
  shap_obj <- get("shap_result", envir = load_env)

  if (!inherits(results_obj, "spec_curve")) {
    stop("Cached 'results' is not a spec_curve object.")
  }
  if (!is.list(shap_obj) || is.null(shap_obj$shapley) || is.null(shap_obj$predictions)) {
    stop("Cached 'shap_result' is missing shapley and/or predictions.")
  }

  list(results = results_obj, shap = shap_obj)
}

cache <- load_paper_cache()

if (is.null(cache)) {
  message("Paper cache not found; recomputing Basque full grid (paper spec space).")
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
    expected_direction = "negative",
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
    cores = max(1L, min(4L, parallel::detectCores(logical = FALSE))),
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
} else {
  message("Using cached paper Basque results/shap objects.")
  results <- cache$results
  shap_result <- cache$shap
}

plot_obj <- plot_spec_curve(
  long_data = results,
  name_treated_unit = treated_unit,
  outcomes = outcome,
  show_shap = TRUE,
  shap_label_type = "signed",
  richtext_feature_labels = FALSE,
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
  device = "svg",
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
