# SCMs: Synthetic Control Methods with Advanced Analytics

A comprehensive R package for Synthetic Control Methods (SCM) featuring specification curve analysis, machine learning interpretability, and robust statistical inference.

## Features

### Core Synthetic Control Methods
- **Multiple Estimation Methods**: OLS, Ridge, Lasso, and Augmented Synthetic Control
- **Flexible Data Processing**: Covariate processing with period-specific and aggregated features
- **Robust Inference**: Placebo-based inference and bootstrap methods
- **Parallel Processing**: High-performance computing support for large-scale analyses

### Advanced Analytics
- **Specification Curve Analysis**: Systematic exploration of researcher degrees of freedom
- **CatBoost + SHAP Analysis**: Machine learning approach to understand what drives variation in treatment effects
- **Advanced Visualizations**: Publication-ready plots including specification curves and SHAP distributions

### Data Processing Capabilities
- **Period-Specific Features**: Create separate variables for each time period
- **Custom Period Selection**: Use specific periods or ranges for feature creation
- **Flexible Aggregation**: Mean aggregation across specified periods
- **Robust Data Handling**: Comprehensive error checking and data validation

## Installation

```r
# Install from source (ensure you're in the package directory)
devtools::install()

# Or install with dependencies
devtools::install(dependencies = TRUE)
```

## Quick Start

### Basic Synthetic Control Analysis

```r
library(SCMs)
library(data.table)

# Load the included German reunification dataset
dataset <- fread(system.file("extdata/scpi_germany.csv", package = "SCMs"))

# Create SCM data structure
scm_data <- scdata(
  df = dataset,
  id.var = "country",
  time.var = "year", 
  outcome.var = "gdp",
  period.pre = 1960:1990,
  period.post = 1991:2003,
  unit.tr = "West Germany",
  unit.co = setdiff(unique(dataset$country), "West Germany"),
  covagg = list(
    gdp_each = list(var = "gdp", each = TRUE)
  )
)

# Basic synthetic control estimation
results <- scest(
  data = scm_data,
  w.constr = list(name = "simplex"),
  feature_weights = "uniform"
)

# View results
summary(results)
plot(results)

# Run placebo inference
inference_results <- inference_sc(results, dataset, inference_type = "placebo")

# View inference results
print(inference_results$abadie_significance)
```

### Specification Curve Analysis

```r
# Define specification parameters for German reunification analysis
spec_results <- spec_curve(
  dataset = dataset,
  outcomes = "gdp",
  col_name_unit_name = "country",
  name_treated_unit = "West Germany",
  covagg = list(
    gdp_periods = list(
      var = "gdp", 
      each = TRUE
    ),
    gdp_aggregate = list(
      var = "gdp",
      periods = c(1980, 1985, 1990)
    ),
    invest_final = list(
      var = "invest80",
      periods = 1980
    )
  ),
  treated_period = 1991,
  min_period = 1960,
  end_period = 2003,
  col_name_period = "year",
  feature_weights = c("uniform", "optimize"),
  outcome_models = c("none", "augsynth", "lasso"),
  constraints = list(
    list(name = "simplex"),
    list(name = "lasso")
  ),
  inference_type = "placebo",
  cores = 4
)

# Visualize results with different test statistics
plot_spec_curve(
  spec_results,
  name_treated_unit = "West Germany",
  outcomes = "gdp",
  test_statistic = "treatment_effect"  # Options: 'treatment_effect', 'rmse_ratio', 'normalized_te'
)
```

### CatBoost + SHAP Analysis

```r
# Configure CatBoost analysis for interpretability
catboost_config <- create_catboost_config(
  dataset_name = "german_reunification_analysis",
  treated_unit_name = "West Germany", 
  outcome_filter = "gdp",
  spec_features = c("feat", "outcome_model", "const", "data_sample", "fw"),
  treated_unit_only = TRUE
)

# Run CatBoost + SHAP analysis
shap_results <- run_catboost_shap_analysis(spec_results$results, catboost_config)

# Create specification curve with SHAP coloring
plot_spec_curve(
  spec_results,
  name_treated_unit = "West Germany",
  outcomes = "gdp",
  shap_values = shap_results$shapley,  # Color points by SHAP values
  test_statistic = "treatment_effect"
)

# Create advanced visualizations
plot_shapley_distributions(shap_results)
plot_interaction_heatmap(shap_results)
```

## Covariate Specification Format

The package uses a flexible covariate specification format. Each specification is a named list containing:

```r
covagg = list(
  specification_name = list(
    var = "variable_name",           # Required: variable name in dataset
    each = TRUE,                     # Optional: create separate feature for each period
    periods = c(1980, 1985, 1990),   # Optional: specific periods to use (default: all pre-treatment)
    aggfunc = "mean"                 # Optional: aggregation function (default: "mean")
  )
)
```

### Examples:

```r
# Create separate features for GDP in each year
gdp_each = list(var = "gdp", each = TRUE)

# Use GDP from specific years only 
gdp_specific = list(var = "gdp", periods = c(1980, 1985, 1990))

# Aggregate investment over all pre-treatment periods
investment_mean = list(var = "investment")

# Use trade data from 1985-1990 only
trade_late = list(var = "trade", periods = 1985:1990)
```

## Package Structure

### Core Functions
- `scdata()`: Data preparation and validation for synthetic control analysis
- `scest()`: Main synthetic control estimation function  
- `inference_sc()`: Statistical inference for SC estimates (placebo and bootstrap methods)

### Specification Curve Analysis
- `spec_curve()`: Core specification curve computation engine
- `run_spec_curve_analysis()`: Complete specification curve pipeline with configuration management
- `plot_spec_curve()`: Specification curve visualization with multiple test statistics

### Machine Learning Analytics  
- `run_catboost_shap_analysis()`: CatBoost modeling with SHAP values
- `plot_shapley_distributions()`: SHAP value visualizations
- `calculate_specification_interactions()`: Feature interaction analysis

### Data Processing
- `create_scm_dataset()`: High-level data formatting wrapper
- `process_covariates()`: Flexible covariate aggregation dispatcher

## Key Dependencies

- **Core Analysis**: `CVXR`, `Matrix`, `optimx`
- **Data Processing**: `data.table`, `dplyr`, `tidyr`, `stringr`
- **Machine Learning**: `catboost`
- **Visualization**: `ggplot2`, `cowplot`, `ggrepel`
- **Parallel Computing**: `foreach`, `doParallel`

## Advanced Features

### Specification Curve Analysis
- Systematic exploration of all reasonable model specifications
- Parallel processing for large specification spaces (1000+ specifications)
- Comprehensive visualization with multiple test statistics
- Statistical summaries and robustness checks

### CatBoost + SHAP Integration
- Leave-one-out cross-validation for robust model evaluation
- SHAP values for interpretable machine learning
- Feature importance analysis across specifications
- Advanced interaction effect visualization

### Flexible Inference Methods
- **Placebo Tests**: Comprehensive placebo-based inference with multiple test statistics
  - RMSE ratio tests
  - Treatment effect tests  
  - Normalized treatment effect tests
- **Bootstrap Methods**: Robust uncertainty quantification

### Multiple Test Statistics
The package supports three types of placebo test statistics:

1. **RMSE Ratio**: Post/pre-treatment RMSE ratio (Abadie et al. 2010)
2. **Treatment Effect**: Average post-treatment effect size
3. **Normalized Treatment Effect**: Treatment effect normalized by pre-treatment RMSE

## Examples

The package includes complete working examples:

- **Quick Start**: `inst/examples/quick_start_example.R` - Basic usage with German data
- **Complete Pipeline**: `inst/examples/german_reunification_example.R` - Full workflow demonstration

Both examples use the included `scpi_germany.csv` dataset and demonstrate:
- Basic synthetic control estimation with different constraints
- Multiple test statistics and inference methods
- Specification curve analysis with SHAP interpretability
- Advanced visualizations and robustness checks

## Data Requirements

Your dataset should be in long format with:
- **Unit identifier**: Column identifying each unit (e.g., countries, states)
- **Time identifier**: Column identifying time periods (numeric)
- **Outcome variable**: Numeric outcome of interest
- **Covariates**: Additional numeric variables for matching

Example data structure:
```
  country  year   gdp  trade  investment
1 Germany  1960  100    50      20
2 Germany  1961  105    55      22
3 France   1960   95    45      18
4 France   1961  100    48      19
```

## Performance Notes

- **Specification Curves**: Can handle 1000+ specifications efficiently with parallel processing
- **Memory Usage**: Scales well with dataset size due to efficient data.table operations
- **Recommended**: Use multiple cores (`cores = 4` or higher) for large specification spaces

## Citation

If you use this package in your research, please cite:

```
Jelveh, Z. (2024). SCMs: Synthetic Control Methods with Advanced Analytics. 
R package version 0.1.
```

## License

GPL-3

## Contact

Zubin Jelveh - zjelveh@umd.edu