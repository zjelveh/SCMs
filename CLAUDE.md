# SCMs R Package - Project Context

## 🚨 CRITICAL DEVELOPMENT PRINCIPLE - HIGHEST PRIORITY 🚨

**FAIL HARD, NOT SOFT**
- **NO FALLBACKS**: Never implement fallback mechanisms or graceful degradation without explicit user approval
- **NO APPROXIMATIONS**: Never substitute approximations for exact functionality (e.g., correlation-based "interactions" when user requests true SHAP interactions)
- **HARD FAILURES ONLY**: When functionality cannot be implemented as requested, fail with clear error messages
- **SOLVE THE ROOT CAUSE**: Fix underlying issues rather than working around them
- **ASK PERMISSION**: If exact functionality is impossible, explicitly ask user for permission to use alternatives
- **PRIORITY**: Maximizing long-term robustness over short-term functionality
- **NO BACKWARD COMPATIBILITY**: Do not write backward compatible code. This repo has not been released so we don't need to worry about breaking changes. Remove any existing backward compatibility code.

## Project Overview
This is an R package for Synthetic Control Methods (SCM) with advanced analytics including specification curve analysis, CatBoost+SHAP integration, and comprehensive inference methods.
The highlight of the package is the specification curve where the feature for each specification are colored by their shapley value. 
The shapley value represents how much that feature is relevant to a predicted treatment effect.

## Code Review Focus Areas
- **Code Organization**: Ensure functions are in logical files, remove redundancy
- **Efficiency**: Identify verbose or inefficient code patterns
- **Refactoring**: Explore whether it makes sense to break supporting_functions.R into two or more separate files?
- **Documentation**: Verify docstrings match function behavior, update man pages
- **Package Integrity**: Ensure DESCRIPTION dependencies are accurate and minimal
- **File Management**: Remove debugging/temporary files and unused functions
- **Backward Compatibility Removal**: Identify and remove all backward compatibility code, fallback mechanisms, and legacy support since this package has not been released

## Key Package Structure
- **Core SC Methods**: `scest.R`, `estimate_sc.R`, `inference_*.R`
- **Data Processing**: `create_scm_dataset.R`, `covariate_processing.R`, `data_processing.R`
- **Specification Curve**: `spec_curve.R`, `spec_curve_analysis.R`
- **Machine Learning**: `catboost_shap_analysis.R` (note: moved from xgboost)
- **Visualization**: `plot_spec_curve.R`, `visualization_advanced.R`
- **Supporting**: `supporting_functions.R`, `outcome_models.R`

## Documentation Standards
- All exported functions must have complete roxygen2 documentation
- Examples should be runnable and demonstrate key functionality
- Parameter descriptions should be precise and include expected data types

## Dependencies to Verify
- Remove unused imports from DESCRIPTION
- Ensure all used packages are listed
- Check version requirements are appropriate
- Note: Package uses catboost, not xgboost (recent change)

## Clarabel Solver Integration
The package supports the high-performance clarabel solver as an alternative to CVXR and kernlab::ipop for optimization tasks.

### Usage
```r
# Enable clarabel for all optimization (W weights and V weights)
options(SCMs.prefer_clarabel = TRUE)

# Use with any SCM operation
result <- scest(data, w.constr = list(name = "simplex"))
result <- scest(data, w.constr = list(name = "ols"))
result <- scest(data, feature_weights = "optimize")  # Uses clarabel for V optimization
```

### Supported Operations
- **W Optimization**: OLS, Simplex constraints (Ridge, Lasso, L1-L2 planned)
- **V Optimization**: Feature weight optimization when `feature_weights = "optimize"`
- **Pensynth**: Penalized synthetic control with cross-validation

### Performance Benefits
- Faster convergence for large problems
- Better numerical stability
- Unified solver architecture across all constraint types
- Sparse matrix support for efficiency

### Solver Tracking and Monitoring
The package tracks which solver was used for each optimization, especially useful for parallel execution:

```r
# Enable clarabel globally
options(SCMs.prefer_clarabel = TRUE)

# Run SCM estimation
result <- scest(data, w.constr = list(name = "simplex"))

# Check which solvers were used
solver_info <- get_solver_info(result)
print_solver_info(solver_info)

# For specification curve analysis
spec_results <- run_spec_curve_analysis(...)
solver_summary <- collect_solver_info(spec_results)
print(solver_summary)
```

### Solver Information Tracking
- **W Optimization**: Tracks solver used, constraint type, and fallback reasons
- **V Optimization**: Tracks solver used for feature weight optimization
- **Parallel Execution**: Fallback warnings and solver usage visible in parallel runs
- **Specification Curves**: Aggregate solver statistics across all specifications

## Build Commands
- Only run these when necessary.
- `devtools::document()` - Update documentation
- `devtools::check()` - Package checks
- `devtools::install()` - Install package

IMPORTANT: 
- Prioritize 100% success versus graceful failure. Avoid graceful failure if it conflicts the goal of 100% success. 
- Do not circumvent issues just to make the code run. Solve the issue.
- Do not make local fixes if a global fix would lead to a more robust codebase
- **CRITICAL**: Do not implement approximations or workarounds without explicit user approval. When the user requests specific functionality (e.g., "true SHAP interactions"), only implement that exact functionality or explicitly ask for permission to use alternatives. Do not silently substitute correlation-based approximations for true interaction effects.


