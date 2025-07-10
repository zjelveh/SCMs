#' Quick Start Example - German Reunification
#' 
#' A minimal example showing basic SCMs package usage with the German
#' reunification dataset included in the package.

library(SCMs)
library(data.table)

# Load the included German dataset
dataset <- fread(system.file("extdata/scpi_germany.csv", package = "SCMs"))

# Basic synthetic control estimation using the main user interface
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
inference_results <- inference_sc(
  results, 
  dataset, 
  inference_type = "placebo"
)

# The inference results now include three test statistics:
# - RMSE ratio (traditional Abadie approach) 
# - Treatment effect (direct comparison)
# - Normalized treatment effect (effect / pre-period SD)

print("P-values for different test statistics:")
print(inference_results$abadie_significance)