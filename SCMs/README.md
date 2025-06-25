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
- **XGBoost + SHAP Analysis**: Machine learning approach to understand what drives variation in treatment effects
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

# Load and prepare data
data <- create_scm_dataset(your_data, 
                          unit_col = "unit_id", 
                          time_col = "year",
                          outcome_col = "outcome")

# Run synthetic control estimation  
results <- scest(data, treated_unit = "TREATED_ID", treated_period = 2015)

# View results
summary(results)
plot(results)
```

### Specification Curve Analysis

```r
# Create analysis configuration
config <- create_analysis_config("study_name", custom_params = list(
  name_treated_unit = "TREATED_ID",
  treated_period = 2015,
  min_period = 2010,
  end_period = 2020,
  cores = 4
))

# Run comprehensive specification curve analysis
spec_results <- run_spec_curve_analysis(config)

# Visualize results
plot_spec_curve(spec_results)
```

### XGBoost + SHAP Analysis

```r
# Configure XGBoost analysis
xgb_config <- create_xgboost_config(
  dataset_name = "my_study",
  file_path = "spec_curve_results.csv",
  treated_unit_name = "TREATED_ID",
  outcome_filter = "outcome_variable"
)

# Run XGBoost + SHAP analysis
xgb_results <- run_xgboost_shap_analysis(xgb_config)

# Create advanced visualizations
plot_shapley_distributions(xgb_results)
plot_interaction_heatmap(xgb_results)
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
- `run_xgboost_shap_analysis()`: XGBoost modeling with SHAP values
- `plot_shapley_distributions()`: SHAP value visualizations
- `calculate_specification_interactions()`: Feature interaction analysis

### Data Processing
- `process_covariates()`: Multiple covariate aggregation methods
- `process_averages_dt()`, `process_growth_dt()`, `process_volatility_dt()`: Specific processing functions
- `process_rolling_dt()`, `process_ranges_dt()`: Time-based processing

## Key Dependencies

- **Core Analysis**: `CVXR`, `Matrix`, `optimx`
- **Data Processing**: `data.table`, `dplyr`, `tidyr`
- **Machine Learning**: `xgboost`, `fastshap`
- **Visualization**: `ggplot2`, `cowplot`
- **Parallel Computing**: `foreach`, `doParallel`, `doSNOW`

## Advanced Features

### Specification Curve Analysis
- Systematic exploration of all reasonable model specifications
- Parallel processing for large specification spaces
- Comprehensive visualization of results
- Statistical summaries and robustness checks

### XGBoost + SHAP Integration
- Leave-one-out cross-validation for robust model evaluation
- SHAP values for interpretable machine learning
- Feature importance analysis across specifications
- Advanced interaction effect visualization

### Flexible Inference Methods
- **SCPI**: Prediction intervals with finite-sample guarantees
- **Placebo Tests**: Comprehensive placebo-based inference
- **Bootstrap Methods**: Robust uncertainty quantification

## Examples

See `examples/complete_pipeline_example.R` for a comprehensive workflow demonstration combining all package features.

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