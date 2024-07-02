# SCMs
The `SCMs` (Synthetic Control Methods) package provides an R implementation for exploring multiple specifications of synthetic control methods. This package is a fork of the original `scpi` package, with additional functionality incorporated from the `Synth` and `AugSynth` packages. The name "SCMs" reflects its capability to handle and compare various synthetic control method specifications.


## Features
- Implementation of the original Synthetic Control Method
- Augmented Synthetic Control Methods (from `AugSynth`)
- Additional outcome models (ols, ridge, lasso)
- Various weight constraint options (ols, ridge, lasso)
- Optimized feature weighting (from `Synth`)
- Placebo inference
- Specification curve visualiztion for SCM results

## Installation
To install/update in R from GitHub:
```R
# install.packages("devtools")
devtools::install_github("zjelveh/SCMs")
```


## Usage
Here are basic examples of how to use the package:

Estimating a single synthetic control:

```R
library(SCMs)

# Estimate the synthetic control
sc_result <- estimate_sc(
    dataset = scm_data,
    outcome = "gdp",
    covagg = list(c("population", "education")),
    col_name_unit_name = "state",
    name_treated_unit = "California",
    col_name_period = "year",
    treated_period = 2000,
    min_period = 1990,
    end_period = 2010,
    outcome_models = c("none", "ridge", "lasso"),
    feature_weights = "optimize"
)

# Perform inference
inference_result <- inference_sc(
    sc_result, 
    dataset = my_data,
    inferece_type='placeb'
)
```

Exploring multiple SCM specifications:

```R
library(SCMs)
# Define multiple specifications
spec_curve_results <- spec_curve(
    dataset = my_data,
    outcomes = c("gdp", "unemployment"),
    col_name_unit_name = "state",
    name_treated_unit = "California",
    covagg = list(
        every_other_period = c("population", "education"),
        every_period = c("income", "inflation")
        ),
    treated_period = 2000,
    min_period = 1990,
    end_period = 2010,
    col_name_period = "year",
    feature_weights = c('uniform', 'optimize'),
    num_pre_period_years = c(NA, 5, 10),
    outcome_models = c('none', 'augsynth', 'ridge', 'lasso', 'ols'),
    donor_sample = c('all', 'most_similar'),
    constraints = list(
        list(name='simplex'),
        list(name='lasso'),
        list(name='ridge')
    )
)

# Plot the specification curve
plot_spec_curve(
    spec_curve_results,
    outcomes = c("gdp", "unemployment"),
    name_treated_unit = "California",
    normalize_outcomes = TRUE,
    rmse_threshold = 0.1
)
```


## Documentation
For detailed documentation on each function, please refer to the help files in R:
```R
?estimate_sc
?inference_sc
?plot_spec_curve
```

## Contributing
Contributions to the SCMs package are welcome. Please feel free to submit a Pull Request.

## References
### Original SCPI Package

Cattaneo, Feng, Palomba and Titiunik (2024): scpi: Uncertainty Quantification for Synthetic Control Methods.
Journal of Statistical Software, forthcoming.

### Synth Package

Abadie, A., Diamond, A., Hainmueller, J. (2011). Synth: An R Package for Synthetic Control Methods in Comparative Case Studies. Journal of Statistical Software, 42(13), 1–17.

### AugSynth Package

Ben-Michael, E., Feller, A., & Rothstein, J. (2021). The Augmented Synthetic Control Method. Journal of the American Statistical Association, 116(536), 1789-1803.

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.
You can copy this content, paste it into a text editor, and then use find-and-replace to convert the HTML entities back to their original characters:

