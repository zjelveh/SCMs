library(devtools)
load_all('.')

cat("=== DEBUGGING PENSYNTH CV FUNCTION ===\n\n")

# Create test data similar to what the spec curve would use
set.seed(123)

# Simulate what create_scm_dataset would create
test_data <- data.frame(
  regionname = rep(c("treated", "control1", "control2", "control3"), each = 20),
  year = rep(1:20, 4),
  gdpcap = rnorm(80, mean = rep(c(10, 8, 9, 7), each = 20), sd = 1)
)

cat("Test data created:\n")
print(head(test_data, 8))

# Create scdata object like the spec curve does
tryCatch({
  scm_data <- create_scm_dataset(
    dataset = test_data,
    outcome = "gdpcap",
    covagg = list(per_period = "gdpcap"),
    col_name_unit_name = "regionname",
    name_treated_unit = "treated",
    col_name_period = "year", 
    treated_period = 11,
    min_period = 1,
    end_period = 20
  )
  
  cat("create_scm_dataset: SUCCESS\n")
  cat("Data structure created:\n")
  cat("- Y.donors dim:", dim(scm_data$Y.donors), "\n")
  cat("- Y.pre dim:", dim(scm_data$Y.pre), "\n")
  cat("- A dim:", length(scm_data$A), "\n")
  cat("- B dim:", dim(scm_data$B), "\n")
  
  cat("\nY.donors content:\n")
  print(scm_data$Y.donors)
  cat("\nY.pre content:\n") 
  print(scm_data$Y.pre)
  
}, error = function(e) {
  cat("create_scm_dataset FAILED:", e$message, "\n")
  stop("Cannot proceed without scdata")
})

# Test the CV function step by step
cat("\n=== TESTING CV FUNCTION STEP BY STEP ===\n")

Y.donors <- scm_data$Y.donors
Y.pre <- scm_data$Y.pre

# Create V matrix
V <- diag(nrow(scm_data$B))
cat("V matrix created, dim:", dim(V), "\n")

# Test data preparation in CV function
cat("\nTesting data preparation...\n")
if (!is.matrix(Y.donors)) Y.donors <- as.matrix(Y.donors)
if (!is.matrix(Y.pre)) Y.pre <- as.matrix(Y.pre)

T_pre <- nrow(Y.donors)
cat("T_pre =", T_pre, "\n")

if (T_pre < 3) {
  cat("ERROR: Not enough periods for CV\n")
  stop()
}

# Test lambda sequence creation
lambda_seq <- exp(seq(log(1e-6), log(10), length.out = 10))
cat("Lambda sequence:", paste(round(lambda_seq, 6), collapse = ", "), "\n")

# Test single iteration of CV loop
cat("\nTesting single CV iteration (holdout period 1)...\n")
t_holdout <- 1
train_periods <- setdiff(1:T_pre, t_holdout)
cat("Train periods:", paste(train_periods, collapse = ", "), "\n")

Y_train_donors <- Y.donors[train_periods, , drop = FALSE]
Y_train_treated <- Y.pre[train_periods, , drop = FALSE]
Y_val_donors <- Y.donors[t_holdout, , drop = FALSE]
Y_val_treated <- Y.pre[t_holdout, , drop = FALSE]

cat("Training data dims:\n")
cat("- Y_train_donors:", dim(Y_train_donors), "\n")
cat("- Y_train_treated:", dim(Y_train_treated), "\n")
cat("Validation data dims:\n")
cat("- Y_val_donors:", dim(Y_val_donors), "\n") 
cat("- Y_val_treated:", dim(Y_val_treated), "\n")

# Test A, Z conversion
A_train <- colMeans(Y_train_treated)
Z_train <- t(Y_train_donors)

cat("Converted training data:\n")
cat("- A_train length:", length(A_train), "value:", A_train, "\n")
cat("- Z_train dim:", dim(Z_train), "\n")
cat("- Z_train content:\n")
print(Z_train)

# Test single lambda fit
cat("\nTesting single lambda fit with lambda =", lambda_seq[1], "...\n")
tryCatch({
  weights <- fit_pensynth_single_lambda(A_train, Z_train, V, lambda_seq[1])
  cat("Single lambda fit: SUCCESS\n")
  cat("Weights:", paste(round(weights, 3), collapse = ", "), "\n")
  
  # Test prediction
  pred <- Y_val_donors %*% weights
  error <- as.numeric((Y_val_treated - pred)^2)
  cat("Prediction:", pred, "Actual:", Y_val_treated, "Error:", error, "\n")
  
}, error = function(e) {
  cat("Single lambda fit FAILED:", e$message, "\n")
  cat("Detailed error investigation:\n")
  
  # Check if it's a V matrix issue
  cat("V matrix properties:\n")
  cat("- V dim:", dim(V), "\n")
  cat("- V diagonal:", paste(round(diag(V)[1:5], 3), collapse = ", "), "...\n")
  cat("- Any NA in V:", any(is.na(V)), "\n")
  cat("- Any Inf in V:", any(is.infinite(V)), "\n")
  
  # Check if it's a penalty calculation issue
  cat("\nPenalty calculation debug:\n")
  tryCatch({
    J <- ncol(Z_train)
    cat("J (num donors):", J, "\n")
    
    v_weights <- diag(V)
    cat("v_weights length:", length(v_weights), "sample:", paste(round(v_weights[1:3], 3), collapse = ", "), "\n")
    
    Z_weighted <- Z_train * sqrt(v_weights)
    A_weighted <- A_train * sqrt(v_weights)
    cat("Weighted data created\n")
    
    Delta <- apply(Z_weighted - c(A_weighted), 2, crossprod)
    cat("Delta computed, length:", length(Delta), "values:", paste(round(Delta, 3), collapse = ", "), "\n")
    
    if (any(is.na(Delta)) || any(is.infinite(Delta))) {
      cat("ERROR: Delta contains NA or Inf values\n")
    }
    
  }, error = function(e2) {
    cat("Penalty calculation failed:", e2$message, "\n")
  })
  
  # Check if it's a CVXR issue
  cat("\nTesting basic CVXR setup:\n")
  tryCatch({
    x <- CVXR::Variable(J)
    base_objective <- CVXR::quad_form(A_train - Z_train %*% x, diag(length(A_train)))
    cat("Basic CVXR objective created\n")
    
    constraints <- list(CVXR::sum_entries(x) == 1, x >= 0)
    prob <- CVXR::Problem(CVXR::Minimize(base_objective), constraints)
    cat("Basic CVXR problem created\n")
    
    sol <- CVXR::solve(prob, solver = "ECOS", verbose = TRUE)
    cat("Basic CVXR solve status:", sol$status, "\n")
    
  }, error = function(e3) {
    cat("Basic CVXR failed:", e3$message, "\n")
  })
})

cat("\n=== DEBUGGING COMPLETE ===\n")
cat("Please run this and report the detailed error information.\n")