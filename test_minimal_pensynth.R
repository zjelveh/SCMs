library(devtools)
load_all('.')

# Create minimal test data - just 3 units, 5 pre periods, 2 post periods
# All data present, no missing values
test_data <- data.frame(
  unit = rep(c("treated", "control1", "control2"), each = 7),
  time = rep(1:7, 3),
  outcome = c(
    # Treated unit (gets treatment at time 6)
    1.0, 1.1, 1.2, 1.3, 1.4, 2.0, 2.5,  # Jump at treatment
    # Control 1 - steady growth
    0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
    # Control 2 - different trend  
    1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8
  )
)

print("Test data:")
print(test_data)

# Create scdata
test_scdata <- scdata(
  df = test_data,
  id.var = 'unit',
  outcome.var = 'outcome', 
  time.var = 'time',
  period.pre = 1:5,
  period.post = 6:7,
  unit.tr = 'treated',
  unit.co = c('control1', 'control2')
)

print("scdata created successfully")
print(paste("Y.donors dim:", paste(dim(test_scdata$Y.donors), collapse="x")))
print(paste("Y.pre dim:", paste(dim(test_scdata$Y.pre), collapse="x")))

# Test pensynth with fixed lambda
print("Testing pensynth with fixed lambda...")
result <- scest(data = test_scdata, w.constr = list(name = 'pensynth', lambda = 0.1))
print("SUCCESS!")
print(paste("Weights:", paste(round(result$est.results$w, 3), collapse = ", ")))