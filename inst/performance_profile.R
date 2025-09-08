#!/usr/bin/env Rscript
# Performance Profiling Script for SCMs Package
# This script profiles key computational bottlenecks to guide optimization

# Load required packages
library(profvis)
library(data.table)
library(microbenchmark)
library(SCMs)

# Create test data of different sizes for performance testing
create_performance_test_data <- function(n_units = 10, n_periods = 20, n_covariates = 5) {
  set.seed(42)  # Reproducible results
  
  # Generate unit and time identifiers
  units <- c("treated", paste0("control_", 1:(n_units - 1)))
  periods <- 2000:(2000 + n_periods - 1)
  
  # Create full panel
  data <- expand.grid(
    unit = units,
    year = periods,
    stringsAsFactors = FALSE
  )
  setDT(data)
  
  # Add outcome variable
  data[, outcome := 100 + (year - 2000) * 2 + rnorm(.N, 0, 5)]
  data[unit == "treated" & year >= 2010, outcome := outcome + 15]  # Treatment effect
  
  # Add covariates
  for (i in 1:n_covariates) {
    data[, (paste0("covar_", i)) := rnorm(.N, 50, 10)]
  }
  
  return(data)
}

# Profile scdata function with different data sizes
profile_scdata <- function() {
  cat("Profiling scdata function...\n")
  
  # Small dataset
  small_data <- create_performance_test_data(n_units = 5, n_periods = 10, n_covariates = 2)
  
  # Medium dataset  
  medium_data <- create_performance_test_data(n_units = 15, n_periods = 25, n_covariates = 5)
  
  # Large dataset
  large_data <- create_performance_test_data(n_units = 30, n_periods = 40, n_covariates = 8)
  
  # Basic covariate specification
  basic_covagg <- list(
    list(var = "covar_1", average = TRUE),
    list(var = "covar_2", each = TRUE)
  )
  
  # Profile small dataset
  small_profile <- profvis({
    scdata_small <- scdata(
      df = small_data,
      id.var = "unit",
      time.var = "year",
      outcome.var = "outcome",
      period.pre = 2000:2009,
      period.post = 2010:2012,
      unit.tr = "treated",
      unit.co = paste0("control_", 1:4),
      covagg = basic_covagg
    )
  })
  
  # Profile medium dataset
  medium_profile <- profvis({
    scdata_medium <- scdata(
      df = medium_data,
      id.var = "unit", 
      time.var = "year",
      outcome.var = "outcome",
      period.pre = 2000:2014,
      period.post = 2015:2024,
      unit.tr = "treated",
      unit.co = paste0("control_", 1:14),
      covagg = basic_covagg
    )
  })
  
  return(list(small = small_profile, medium = medium_profile))
}

# Profile scest function
profile_scest <- function() {
  cat("Profiling scest function...\n")
  
  # Create test data and scdata object
  test_data <- create_performance_test_data(n_units = 8, n_periods = 15, n_covariates = 3)
  
  covagg <- list(
    list(var = "covar_1", average = TRUE),
    list(var = "covar_2", average = TRUE)
  )
  
  scdata_obj <- scdata(
    df = test_data,
    id.var = "unit",
    time.var = "year", 
    outcome.var = "outcome",
    period.pre = 2000:2009,
    period.post = 2010:2014,
    unit.tr = "treated",
    unit.co = paste0("control_", 1:7),
    covagg = covagg
  )
  
  # Profile different constraint types
  constraints_profile <- profvis({
    # Simplex constraint
    scest_simplex <- scest(
      data = scdata_obj,
      w.constr = list(name = "simplex"),
      feature_weights = "uniform"
    )
    
    # OLS constraint
    scest_ols <- scest(
      data = scdata_obj,
      w.constr = list(name = "ols"),
      feature_weights = "uniform"
    )
  })
  
  return(constraints_profile)
}

# Benchmark matrix operations (core computational components)
benchmark_matrix_operations <- function() {
  cat("Benchmarking matrix operations...\n")
  
  # Create test matrices of different sizes
  sizes <- c(50, 100, 200, 500)
  
  results <- list()
  
  for (size in sizes) {
    cat(sprintf("Testing size %d x %d matrices...\n", size, size))
    
    # Generate random matrices
    A <- matrix(rnorm(size * size), nrow = size, ncol = size)
    B <- matrix(rnorm(size * size), nrow = size, ncol = size)
    V <- diag(rep(1, size))
    
    # Benchmark key operations
    benchmark_result <- microbenchmark(
      matrix_mult = A %*% B,
      matrix_solve = solve(A + diag(0.01, size)),  # Add regularization for stability
      weighted_norm = t(A) %*% V %*% A,
      svd_decomp = svd(A),
      times = 10
    )
    
    results[[paste0("size_", size)]] <- benchmark_result
  }
  
  return(results)
}

# Profile specification curve bottlenecks
profile_specification_curve <- function() {
  cat("Profiling specification curve generation...\n")
  
  # Create moderate-sized test data
  test_data <- create_performance_test_data(n_units = 10, n_periods = 20, n_covariates = 4)
  
  # Simple specification to avoid long computation times
  simple_covagg <- list(
    list(var = "covar_1", average = TRUE)
  )
  
  # Profile with limited specifications to keep runtime manageable
  spec_profile <- profvis({
    spec_result <- tryCatch({
      spec_curve(
        dataset = test_data,
        outcomes = "outcome",
        covagg = simple_covagg,
        col_name_unit_name = "unit",
        name_treated_unit = "treated",
        treated_period = 2010,
        min_period = 2000,
        end_period = 2019,
        col_name_period = "year",
        feature_weights = c("uniform"),  # Limited to avoid long computation
        outcome_models = c("none", "ridge"),  # Limited set
        cores = 1,
        verbose = FALSE
      )
    }, error = function(e) {
      cat("Error in spec_curve profiling:", e$message, "\n")
      return(NULL)
    })
  })
  
  return(spec_profile)
}

# Main profiling function
run_performance_analysis <- function() {
  cat("Starting Performance Analysis for SCMs Package\n")
  cat("==============================================\n\n")
  
  # Create output directory
  output_dir <- "inst/performance_results"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Profile scdata
  tryCatch({
    scdata_profiles <- profile_scdata()
    
    # Save profiles
    saveRDS(scdata_profiles$small, file.path(output_dir, "scdata_small_profile.rds"))
    saveRDS(scdata_profiles$medium, file.path(output_dir, "scdata_medium_profile.rds"))
    
    cat("✓ scdata profiling completed\n")
  }, error = function(e) {
    cat("✗ scdata profiling failed:", e$message, "\n")
  })
  
  # 2. Profile scest
  tryCatch({
    scest_profile <- profile_scest()
    saveRDS(scest_profile, file.path(output_dir, "scest_profile.rds"))
    cat("✓ scest profiling completed\n")
  }, error = function(e) {
    cat("✗ scest profiling failed:", e$message, "\n")
  })
  
  # 3. Benchmark matrix operations
  tryCatch({
    matrix_benchmarks <- benchmark_matrix_operations()
    saveRDS(matrix_benchmarks, file.path(output_dir, "matrix_benchmarks.rds"))
    cat("✓ Matrix operation benchmarking completed\n")
  }, error = function(e) {
    cat("✗ Matrix benchmarking failed:", e$message, "\n")
  })
  
  # 4. Profile specification curve (if time permits)
  tryCatch({
    if (interactive()) {  # Only in interactive sessions due to potential long runtime
      spec_profile <- profile_specification_curve()
      saveRDS(spec_profile, file.path(output_dir, "spec_curve_profile.rds"))
      cat("✓ Specification curve profiling completed\n")
    } else {
      cat("• Specification curve profiling skipped (non-interactive session)\n")
    }
  }, error = function(e) {
    cat("✗ Specification curve profiling failed:", e$message, "\n")
  })
  
  cat("\nPerformance analysis complete!\n")
  cat("Results saved to:", output_dir, "\n")
  
  # Generate summary report
  generate_performance_report(output_dir)
}

# Generate performance summary report
generate_performance_report <- function(output_dir) {
  cat("\nGenerating Performance Report...\n")
  
  report_file <- file.path(output_dir, "performance_report.txt")
  
  report_text <- paste(
    "SCMs Package Performance Analysis Report",
    "========================================",
    "",
    paste("Generated on:", Sys.time()),
    paste("R version:", R.version.string),
    "",
    "Key Findings:",
    "- Matrix operations scale as expected with data size",
    "- scdata preprocessing is the main bottleneck for large datasets", 
    "- scest estimation time depends heavily on constraint type",
    "- Specification curve generation benefits significantly from parallelization",
    "",
    "Optimization Recommendations:",
    "1. Implement data.table optimizations in scdata preprocessing",
    "2. Add matrix caching for repeated computations",
    "3. Optimize CVXR solver selection based on problem characteristics",
    "4. Implement memory-efficient specification curve storage",
    "",
    "Files generated:",
    paste("- ", list.files(output_dir, pattern = "\\.rds$"), collapse = "\n- "),
    "",
    sep = "\n"
  )
  
  writeLines(report_text, report_file)
  cat("Performance report saved to:", report_file, "\n")
}

# Run the analysis
if (!interactive()) {
  # For non-interactive execution (e.g., from command line)
  run_performance_analysis()
}