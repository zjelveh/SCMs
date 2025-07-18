library(devtools)
load_all('.')

# Create simple synthetic data for testing
set.seed(123)
n_units <- 5
n_periods <- 10
n_post <- 5

# Create synthetic panel data - ensure complete data for all units/periods
test_data <- expand.grid(
  unit = paste0("unit_", 1:n_units),
  time = 1:(n_periods + n_post)
)
test_data$outcome <- rnorm(nrow(test_data), 
                          mean = as.numeric(gsub("unit_", "", test_data$unit)) + test_data$time * 0.1, 
                          sd = 0.1)

print("Test data created")
print(head(test_data))

# Create scdata object
test_scdata <- scdata(
  df = test_data,
  id.var = 'unit',
  outcome.var = 'outcome', 
  time.var = 'time',
  period.pre = 1:n_periods,
  period.post = (n_periods + 1):(n_periods + n_post),
  unit.tr = 'unit_1',
  unit.co = paste0("unit_", 2:n_units)
)

print('Testing pensynth with fixed lambda...')
result <- scest(data = test_scdata, w.constr = list(name = 'pensynth', lambda = 0.1))
print('Success!')
print(paste("Weights:", paste(round(result$est.results$w, 3), collapse = ", ")))