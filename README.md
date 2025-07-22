# SCMs: Synthetic Control Methods with Advanced Analytics

A comprehensive R package for Synthetic Control Methods featuring specification curve analysis with SHAP interpretability. The package's highlight is specification curves where each point is colored by its Shapley value, showing how much each specification feature contributes to the predicted treatment effect.

## Installation

```r
# Install from source (ensure you're in the package directory)
devtools::install()
```

## Quick Start

```r
library(SCMs)
library(data.table)

# Load German reunification dataset
dataset <- fread(system.file("extdata/scpi_germany.csv", package = "SCMs"))

# Create SCM data structure
scm_data <- scdata(
  df = dataset,
  id.var = "country", time.var = "year", outcome.var = "gdp",
  period.pre = 1960:1990, period.post = 1991:2003,
  unit.tr = "West Germany",
  unit.co = setdiff(unique(dataset$country), "West Germany")
)

# Estimate synthetic control
results <- scest(data = scm_data, w.constr = list(name = "simplex"))

# View results and run inference
summary(results)
plot(results)
inference_results <- inference_sc(results, dataset, inference_type = "placebo")
```

## Covariate Aggregation Guide

The package's flexible covariate system allows sophisticated feature engineering for synthetic control matching. Each specification is a named list with variable name and aggregation method:

### Basic Aggregation Types

```r
covagg = list(
  # Average over all pre-treatment periods (default)
  gdp_mean = list(var = "gdp"),
  
  # Separate feature for each period
  gdp_per_period = list(var = "gdp", each = TRUE),
  
  # Use specific periods only
  gdp_late_pre = list(var = "gdp", periods = c(1985, 1988, 1990)),
  
  # Average over full pre-treatment period (explicit)
  trade_average = list(var = "trade", average = "full_pre")
)
```

### Advanced Aggregation Methods

```r
covagg = list(
  # Rolling N-period averages
  gdp_rolling_3yr = list(var = "gdp", rolling = 3),
  
  # First N periods of pre-treatment
  early_investment = list(var = "investment", first = 5),
  
  # Last N periods before treatment
  late_trade = list(var = "trade", last = 3),
  
  # Every N periods (useful for long time series)
  gdp_every_5yr = list(var = "gdp", every_n = 5),
  
  # Period-over-period growth rates
  gdp_growth = list(var = "gdp", growth = "period_over_period"),
  
  # Volatility measures (standard deviation)
  gdp_volatility = list(var = "gdp", volatility = "sd")
)
```

### Complex Aggregation Examples

```r
# Multiple aggregations for same variable
covagg = list(
  # GDP in different forms
  gdp_early = list(var = "gdp", periods = c(1970, 1975, 1980)),
  gdp_trend = list(var = "gdp", rolling = 5),
  gdp_volatility = list(var = "gdp", volatility = "sd"),
  
  # Trade patterns
  trade_baseline = list(var = "trade", periods = 1970),
  trade_growth = list(var = "trade", growth = "period_over_period"),
  
  # Investment timing
  invest_early = list(var = "investment", first = 10),
  invest_late = list(var = "investment", last = 5),
  
  # Population every 5 years
  pop_milestones = list(var = "population", every_n = 5)
)
```

### Realistic Example: Economic Development Analysis

```r
# Comprehensive covariate specification for development analysis
economic_covariates = list(
  # GDP measures
  gdp_1970 = list(var = "gdp", periods = 1970),           # Baseline level
  gdp_1980s = list(var = "gdp", periods = 1980:1989),     # Decade average
  gdp_growth = list(var = "gdp", growth = "period_over_period"), # Growth trend
  gdp_volatility = list(var = "gdp", volatility = "sd"),  # Economic stability
  
  # Trade patterns
  trade_openness = list(var = "trade", rolling = 3),      # Smoothed trend
  trade_late_pre = list(var = "trade", last = 5),         # Recent pattern
  
  # Investment and development
  investment_early = list(var = "investment", first = 10), # Early development
  investment_trend = list(var = "investment", rolling = 5), # Investment momentum
  
  # Demographics and education
  pop_milestones = list(var = "population", every_n = 5),  # Population checkpoints
  education_recent = list(var = "education", last = 3),    # Recent education levels
  
  # Infrastructure
  infrastructure_each = list(var = "infrastructure", each = TRUE) # Period-specific
)
```

## Specification Curve Analysis with SHAP

```r
# Run comprehensive specification curve analysis
spec_results <- spec_curve(
  dataset = dataset,
  outcomes = "gdp",
  col_name_unit_name = "country",
  name_treated_unit = "West Germany",
  covagg = economic_covariates,  # Use complex covariate specification
  treated_period = 1991,
  min_period = 1960,
  end_period = 2003,
  col_name_period = "year",
  feature_weights = c("uniform", "optimize"),
  outcome_models = c("none", "augsynth", "lasso"),
  constraints = list(
    list(name = "simplex"),
    list(name = "lasso", Q = 0.1),
    list(name = "ridge", Q = 0.1)
  ),
  inference_type = "placebo",
  cores = 4
)

# Configure CatBoost + SHAP analysis
catboost_config <- create_catboost_config(
  dataset_name = "german_reunification",
  treated_unit_name = "West Germany",
  outcome_filter = "gdp",
  spec_features = c("feat", "outcome_model", "const", "data_sample", "fw")
)

# Generate SHAP values
shap_results <- run_catboost_shap_analysis(spec_results$results, catboost_config)

# Create specification curve with SHAP coloring
plot_spec_curve(
  spec_results,
  name_treated_unit = "West Germany",
  outcomes = "gdp",
  shap_values = shap_results$shapley,  # Color by SHAP importance
  test_statistic = "treatment_effect"
)
```

## Performance and Advanced Features

- **High-Performance Computing**: Built-in clarabel solver for faster optimization
- **Parallel Processing**: Handle 1000+ specifications efficiently with multi-core support
- **Comprehensive Inference**: Multiple test statistics (RMSE ratio, treatment effect, normalized)
- **Machine Learning Integration**: CatBoost + SHAP for specification interpretability
- **Flexible Data Processing**: Sophisticated covariate aggregation system

Enable high-performance optimization:
```r
# Use clarabel solver for faster computation
options(SCMs.prefer_clarabel = TRUE)
```

## Data Requirements

Dataset should be in long format with:
- **Unit identifier**: Column for units (countries, states, etc.)
- **Time identifier**: Numeric time periods
- **Outcome variable**: Numeric variable of interest
- **Covariates**: Additional numeric matching variables

```
  country  year   gdp  trade  investment
  Germany  1960  100    50      20
  Germany  1961  105    55      22
  France   1960   95    45      18
  France   1961  100    48      19
```

## Examples and Documentation

- **Quick Start**: `inst/examples/quick_start_example.R`
- **Complete Workflow**: `inst/examples/german_reunification_example.R`

Both examples demonstrate the full pipeline from data preparation through specification curve analysis with SHAP interpretability.

## Citation

```
Jelveh, Z. (2024). SCMs: Synthetic Control Methods with Advanced Analytics.
R package version 0.1.0.
```

## Contact

Zubin Jelveh - zjelveh@umd.edu

## License

GPL-3