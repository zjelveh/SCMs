# SCMs: Synthetic Control Methods

R package for synthetic control workflows, including specification-curve analysis and SHAP-based diagnostics.

## Installation

```bash
# From the package root
R CMD INSTALL --no-multiarch --no-demo --no-docs .
```

```r
# Development workflow
library(devtools)
load_all(".")
```

## Quick Start

```r
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
    list(
      var = "outcome_var",
      partition_periods = list(type = "by_period")
    ),
    list(
      var = "cov1",
      compute = "mean"
    )
  )
)

fit <- scest(
  data = scm_data,
  w.constr = list(name = "simplex"),
  feature_weights = "optimize"
)

summary(fit)
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
```

## Covariate Aggregation (`covagg`)

`covagg` tells SCMs which pre-period moments become matching features.

Each `covagg` entry is one operation with up to four fields:
- `var`: variable name (or `"outcome_var"`)
- `select_periods`: which pre-periods are used
- `partition_periods`: how selected periods are grouped into features
- `compute`: aggregation function (`"mean"` by default)

Operation format:

```r
covagg = list(
  list(
    var = "outcome_var",
    partition_periods = list(type = "by_period")
  ),
  list(
    var = "population",
    compute = "mean"
  ),
  list(
    var = "income",
    select_periods = list(type = "explicit", periods = 2001:2004),
    compute = "mean"
  )
)
```

For `spec_curve()`, `covagg` is a named list of specifications, each with `operations`:

```r
covagg = list(
  "Outcome Only" = list(
    label = "Outcome Only",
    operations = list(
      list(var = "outcome_var", partition_periods = list(type = "by_period"))
    )
  ),
  "Outcome + Population Mean" = list(
    label = "Outcome + Population Mean",
    operations = list(
      list(var = "outcome_var", partition_periods = list(type = "by_period")),
      list(var = "population", compute = "mean")
    )
  )
)
```

Rules:
- `select_periods = NULL` means all pre-periods.
- `partition_periods = list(type = "by_period")` creates one feature per period.
- `partition_periods = list(type = "all")` creates one aggregated feature.
- Use `"outcome_var"` to reference the current outcome.

## Specification Curve + SHAP

```r
library(SCMs)

set.seed(42)
d <- expand.grid(
  unit = c("treated", "c1", "c2", "c3", "c4"),
  year = 2001:2010,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Keep IDs simple/alphanumeric for robust matching across internal reshaping paths.
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
    "Outcome Per Period + Pop Mean" = list(
      label = "Outcome Per Period + Pop Mean",
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

xgb_cfg <- create_xgboost_config(
  dataset_name = "toy",
  treated_unit_name = "treated",
  outcome_filter = "outcome",
  spec_features = c("feat", "outcome_model", "const", "fw", "data_sample")
)

shap_results <- run_xgboost_shap_analysis(
  spec_results$results,
  xgb_cfg,
  compute_loo = FALSE
)

plot_obj <- plot_spec_curve(
  long_data = spec_results,
  name_treated_unit = "treated",
  outcomes = "outcome",
  shap_values = shap_results$shapley,
  show_pvalues = TRUE,
  test_statistic = "treatment_effect"
)
```

## Where `covagg` Is Documented

- `COVAGG.md`
- `?scdata`
- `?spec_curve`

## Data Requirements

Input data must be long format with:
- unit identifier
- time identifier
- numeric outcome
- optional numeric covariates

## Contact

Zubin Jelveh - zjelveh@umd.edu

## License

GPL-3
