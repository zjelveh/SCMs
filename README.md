# SCMs: Synthetic Control Methods with Advanced Analytics

A comprehensive R package for Synthetic Control Methods (SCM) featuring specification curve analysis, machine learning interpretability, and robust statistical inference.

## Features

### Core Synthetic Control Methods
- **Multiple Estimation Methods**: OLS, Ridge, Lasso, and Augmented Synthetic Control
- **Flexible Data Processing**: Comprehensive covariate processing with multiple aggregation methods
- **Robust Inference**: SCPI (Synthetic Control Prediction Intervals) and placebo-based inference
- **Parallel Processing**: High-performance computing support for large-scale analyses

### Advanced Analytics
- **Specification Curve Analysis**: Systematic exploration of researcher degrees of freedom
- **CatBoost + SHAP Analysis**: Machine learning approach to understand what drives variation in treatment effects
- **Advanced Visualizations**: Publication-ready plots including specification curves, SHAP distributions, and interaction heatmaps

### Data Processing Capabilities
- **Multiple Aggregation Methods**: Averages, growth rates, volatility measures, rolling windows, and more
- **Flexible Time Processing**: Handle specific periods, ranges, first/last N periods
- **Robust Data Handling**: Comprehensive error checking and data validation

## Installation

```r
# Install from source
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

# Basic synthetic control estimation
results <- scest(
  dataset, 
  treated_unit = "West Germany",
  treated_period = 1991,
  outcome = "gdp",
  unit_col = "country",
  time_col = "year"
)

# View results
summary(results)
plot(results)

# Run placebo inference with multiple test statistics
inference_results <- inference_sc(results, dataset, inference_type = "placebo")

# View p-values for different test statistics
print(inference_results$abadie_significance)
```

### Specification Curve Analysis

```r
# Define specification parameters for German reunification analysis
params <- list(
  cores = 4,
  outcomes = "gdp",
  col_name_unit_name = "country",
  name_treated_unit = "West Germany",
  covagg = list(
    parsimonious = list(
      label = 'GDP + Trade Average',
      gdp_avg = list(var = "gdp", average = "full_pre"),
      trade_avg = list(var = 'trade', average = 'full_pre')
    ),
    gdp_only = list(
      label = 'GDP Only',
      gdp_each = list(var = "gdp", each = TRUE)
    )
  ),
  treated_period = 1991,  # German reunification
  min_period = 1960,
  end_period = 2003,
  col_name_period = "year",
  feature_weights = c("uniform", "optimize"),
  outcome_models = c("none", "augsynth", "lasso"),
  constraints = list(
    list(name = "simplex"),
    list(name = "lasso")
  ),
  inference_type = 'placebo'
)

# Run comprehensive specification curve analysis
spec_results <- run_spec_curve_analysis(dataset = dataset, params = params)

# Visualize results with different test statistics
plot_spec_curve(
  long_data = spec_results,
  name_treated_unit = 'West Germany',
  outcomes = 'gdp',
  test_statistic = 'treatment_effect'  # or 'rmse_ratio', 'normalized_te'
)
```

### CatBoost + SHAP Analysis

```r
# Configure CatBoost analysis for interpretability
catboost_config <- create_catboost_config(
  dataset_name = "german_reunification_analysis",
  treated_unit_name = "West Germany",
  outcome_filter = "gdp",
  spec_features = c('feat', 'outcome_model', 'const', 'data_sample', 'fw'),
  treated_unit_only = TRUE
)

# Run CatBoost + SHAP analysis
shap_results <- run_catboost_shap_analysis(spec_results$results, catboost_config)

# Create specification curve with SHAP coloring
plot_spec_curve(
  long_data = spec_results,
  name_treated_unit = 'West Germany',
  outcomes = 'gdp',
  shap_values = shap_results$shapley,  # Color points by SHAP values
  test_statistic = 'treatment_effect'
)

# Create advanced visualizations
plot_shapley_distributions(shap_results)
plot_interaction_heatmap(shap_results)
```

## Package Structure

### Core Functions
- `scest()`: Main synthetic control estimation function
- `create_scm_dataset()`: Data preparation and validation
- `inference_sc()`: Statistical inference for SC estimates

### Specification Curve Analysis
- `run_spec_curve_analysis()`: Complete specification curve pipeline
- `create_analysis_config()`: Configuration management
- `plot_spec_curve()`: Specification curve visualization

### Machine Learning Analytics
- `run_catboost_shap_analysis()`: CatBoost modeling with SHAP values
- `plot_shapley_distributions()`: SHAP value visualizations
- `calculate_specification_interactions()`: Feature interaction analysis

### Data Processing
- `process_covariates()`: Multiple covariate aggregation methods
- `process_averages_dt()`, `process_growth_dt()`, `process_volatility_dt()`: Specific processing functions
- `process_rolling_dt()`, `process_ranges_dt()`: Time-based processing

## Key Dependencies

- **Core Analysis**: `CVXR`, `Matrix`, `optimx`
- **Data Processing**: `data.table`, `dplyr`, `tidyr`
- **Machine Learning**: `catboost`
- **Visualization**: `ggplot2`, `cowplot`
- **Parallel Computing**: `foreach`, `doParallel`, `doSNOW`

## Advanced Features

### Specification Curve Analysis
- Systematic exploration of all reasonable model specifications
- Parallel processing for large specification spaces
- Comprehensive visualization of results
- Statistical summaries and robustness checks

### CatBoost + SHAP Integration
- Leave-one-out cross-validation for robust model evaluation
- SHAP values for interpretable machine learning
- Feature importance analysis across specifications
- Advanced interaction effect visualization

### Flexible Inference Methods
- **SCPI**: Prediction intervals with finite-sample guarantees
- **Placebo Tests**: Comprehensive placebo-based inference
- **Bootstrap Methods**: Robust uncertainty quantification

## Examples

The package includes complete working examples:

- **Quick Start**: `inst/examples/quick_start_example.R` - Basic usage with German data
- **Complete Pipeline**: `inst/examples/german_reunification_example.R` - Full workflow demonstration combining all package features

Both examples use the included `scpi_germany.csv` dataset and demonstrate:
- Basic synthetic control estimation
- Multiple test statistics (RMSE ratio, treatment effect, normalized treatment effect)  
- Specification curve analysis with SHAP interpretability
- Advanced visualizations

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