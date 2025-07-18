library(devtools)
load_all('.')

cat("=== PENSYNTH DEBUGGING SCRIPT ===\n\n")

# Test 1: Create minimal synthetic data
cat("TEST 1: Creating minimal synthetic data...\n")
set.seed(123)
n_units <- 4
n_periods <- 6
n_post <- 3

data_list <- list()
for (unit in 1:n_units) {
  for (period in 1:(n_periods + n_post)) {
    data_list <- append(data_list, list(data.frame(
      unit = paste0("unit_", unit),
      time = period,
      outcome = rnorm(1, mean = unit + period * 0.1, sd = 0.1)
    )))
  }
}
test_data <- do.call(rbind, data_list)

cat("Synthetic data structure:\n")
print(str(test_data))
cat("Sample data:\n")
print(head(test_data, 12))

# Test 2: Create scdata object and inspect it
cat("\nTEST 2: Creating scdata object...\n")
tryCatch({
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
  cat("scdata creation: SUCCESS\n")
  
  # Inspect scdata structure
  cat("scdata specs:\n")
  print(test_scdata$specs)
  cat("Y.donors dimensions:", dim(test_scdata$Y.donors), "\n")
  cat("Y.pre dimensions:", dim(test_scdata$Y.pre), "\n")
  cat("Y.donors sample:\n")
  print(test_scdata$Y.donors[1:3, 1:min(3, ncol(test_scdata$Y.donors))])
  cat("Y.pre sample:\n")
  print(test_scdata$Y.pre[1:3, ])
  
}, error = function(e) {
  cat("scdata creation FAILED:", e$message, "\n")
  return()
})

# Test 3: Test scest with simplex constraint (baseline)
cat("\nTEST 3: Testing baseline simplex constraint...\n")
tryCatch({
  baseline_result <- scest(data = test_scdata, w.constr = list(name = 'simplex'))
  cat("Baseline simplex: SUCCESS\n")
  cat("Baseline weights:", paste(round(baseline_result$est.results$w, 3), collapse = ", "), "\n")
}, error = function(e) {
  cat("Baseline simplex FAILED:", e$message, "\n")
  return()
})

# Test 4: Test individual pensynth components
cat("\nTEST 4: Testing pensynth components...\n")

# Test lambda CV function directly
cat("Testing lambda CV function...\n")
tryCatch({
  # Extract data from scdata object
  Y.donors <- test_scdata$Y.donors
  Y.pre <- test_scdata$Y.pre
  
  cat("Data for CV:\n")
  cat("Y.donors class:", class(Y.donors), "dim:", dim(Y.donors), "\n")
  cat("Y.pre class:", class(Y.pre), "dim:", dim(Y.pre), "\n")
  
  # Create simple V matrix for testing
  V_test <- diag(nrow(test_scdata$B))
  cat("V matrix dim:", dim(V_test), "\n")
  
  # Test CV function
  optimal_lambda <- estimate_optimal_lambda_cv(Y.donors, Y.pre, V_test, nlambda = 10)
  cat("Lambda CV: SUCCESS, optimal lambda =", optimal_lambda, "\n")
  
}, error = function(e) {
  cat("Lambda CV FAILED:", e$message, "\n")
  cat("Traceback:\n")
  traceback()
})

# Test 5: Test fit_pensynth_single_lambda function
cat("\nTEST 5: Testing single lambda fit...\n")
tryCatch({
  A_test <- test_scdata$A
  Z_test <- test_scdata$B
  V_test <- diag(nrow(Z_test))
  lambda_test <- 0.1
  
  cat("Test data for single lambda:\n")
  cat("A dimensions:", length(A_test), "\n")
  cat("Z dimensions:", dim(Z_test), "\n")
  cat("V dimensions:", dim(V_test), "\n")
  cat("Lambda:", lambda_test, "\n")
  
  weights_test <- fit_pensynth_single_lambda(A_test, Z_test, V_test, lambda_test)
  cat("Single lambda fit: SUCCESS\n")
  cat("Test weights:", paste(round(weights_test, 3), collapse = ", "), "\n")
  
}, error = function(e) {
  cat("Single lambda fit FAILED:", e$message, "\n")
  cat("Traceback:\n")
  traceback()
})

# Test 6: Test full pensynth with fixed lambda
cat("\nTEST 6: Testing pensynth with FIXED lambda...\n")
tryCatch({
  pensynth_result <- scest(data = test_scdata, w.constr = list(name = 'pensynth', lambda = 0.1))
  cat("Pensynth fixed lambda: SUCCESS\n")
  cat("Pensynth weights:", paste(round(pensynth_result$est.results$w, 3), collapse = ", "), "\n")
  
}, error = function(e) {
  cat("Pensynth fixed lambda FAILED:", e$message, "\n")
  cat("Traceback:\n")
  traceback()
})

# Test 7: Test full pensynth with CV lambda
cat("\nTEST 7: Testing pensynth with CV lambda...\n")
tryCatch({
  pensynth_cv_result <- scest(data = test_scdata, w.constr = list(name = 'pensynth'))
  cat("Pensynth CV lambda: SUCCESS\n")
  cat("Pensynth CV weights:", paste(round(pensynth_cv_result$est.results$w, 3), collapse = ", "), "\n")
  cat("Optimal lambda was:", pensynth_cv_result$est.results$w.constr$lambda, "\n")
  
}, error = function(e) {
  cat("Pensynth CV lambda FAILED:", e$message, "\n")
  cat("Traceback:\n")
  traceback()
})

# Test 8: Test with Basque data (if available)
cat("\nTEST 8: Testing with Basque data...\n")
tryCatch({
  # Load basque data
  basque_data <- read.csv('inst/extdata/basque.csv')
  
  # Check data structure
  cat("Basque data loaded, regions:", length(unique(basque_data$regionname)), "\n")
  cat("Years range:", range(basque_data$year), "\n")
  
  # Check for Basque Country
  basque_units <- unique(basque_data$regionname)
  basque_unit_name <- basque_units[grepl("Basque", basque_units)]
  cat("Basque unit name:", basque_unit_name, "\n")
  
  # Try with a simple specification
  basque_scdata <- scdata(
    df = basque_data,
    id.var = 'regionname',
    outcome.var = 'gdpcap', 
    time.var = 'year',
    period.pre = 1960:1969,  # Shorter pre period
    period.post = 1970:1972,  # Very short post period
    unit.tr = basque_unit_name,
    unit.co = basque_units[1:3]  # Just first 3 control units
  )
  
  cat("Basque scdata creation: SUCCESS\n")
  
  # Test with fixed lambda
  basque_result <- scest(data = basque_scdata, w.constr = list(name = 'pensynth', lambda = 0.1))
  cat("Basque pensynth: SUCCESS\n")
  
}, error = function(e) {
  cat("Basque data test FAILED:", e$message, "\n")
})

cat("\n=== DEBUGGING COMPLETE ===\n")
cat("Please run this script and report which tests fail and their error messages.\n")