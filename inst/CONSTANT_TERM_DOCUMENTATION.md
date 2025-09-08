# Constant Term (Intercept) Support in SCMs Package

## Overview

The SCMs package now supports global constant/intercept terms in synthetic control estimation. This feature allows the synthetic control to have a level shift relative to the donor pool, following the methodology used in the SCPI package.

## Key Features

- **Global constant term**: Estimated alongside donor weights
- **All constraint types supported**: Simplex, OLS, Lasso, Ridge, L1-L2
- **Flexible integration**: Works with existing optimization framework
- **Multiple solver support**: CVXR, clarabel, ECOS automatically supported

## Usage

### 1. Data Preparation with Constant

```r
# Create SCM data with constant term enabled
scm_data <- scdata(
  df = your_data,
  id.var = "country",
  time.var = "year", 
  outcome.var = "gdp",
  period.pre = 1960:1990,
  period.post = 1991:2003,
  unit.tr = "treated_unit",
  unit.co = c("donor1", "donor2", "donor3"),
  constant = TRUE,  # Enable constant term
  covagg = list(
    gdp_avg = list(var = "gdp", average = "full_pre"),
    investment_avg = list(var = "investment", average = "full_pre"),
    trade_periods = list(var = "trade", periods = c(1985, 1990))
  )
)
```

**Requirements for constant term:**
- Must specify multiple features (more than one covariate specification)
- Cannot use constant with only a single feature (mathematical non-identifiability)

### 2. Estimation with Different Constraints

#### Simplex Constraint (Standard SC)
```r
result_simplex <- scest(
  data = scm_data, 
  w.constr = list(name = "simplex")
)

# Access results
donor_weights <- result_simplex$est.results$w            # Sum to 1
constant_term <- result_simplex$est.results$constant_term # Unconstrained
```

#### OLS (No Constraints) 
```r
result_ols <- scest(
  data = scm_data,
  w.constr = list(name = "ols")
)

# All coefficients (donors + constant) are unconstrained
constant_term <- result_ols$est.results$constant_term
```

#### Lasso/Ridge with Constant
```r
result_lasso <- scest(
  data = scm_data,
  w.constr = list(name = "lasso", Q = 0.1)
)

# Penalty applies to both donor weights and constant term
constant_term <- result_lasso$est.results$constant_term
```

## Mathematical Framework

### Without Constant
- Optimization: `min ||A - B*w||_V^2 + penalty(w)`
- Variables: `w = [w1, w2, ..., wJ]` (J donor weights)

### With Constant  
- Optimization: `min ||A_exp - B_exp*[w;c]||_V^2 + penalty([w;c])`
- Variables: `[w1, w2, ..., wJ, c]` (J donor weights + constant)
- Matrices expanded:
  - `A_exp = [A_original; 1]` 
  - `B_exp = [B_original; 1,1,...,1]`

### Constraint Behavior

| Constraint Type | Donor Weights | Constant Term |
|----------------|---------------|---------------|
| **Simplex** | `sum(w) = 1, w >= 0` | Unconstrained |
| **OLS** | Unconstrained | Unconstrained |
| **Lasso** | `||w||_1 <= Q` | Included in L1 penalty |
| **Ridge** | `||w||_2^2 <= Q` | Included in L2 penalty |

## Result Structure

After estimation with constant term:

```r
result$est.results$b              # Full coefficient vector [w1,...,wJ,c]
result$est.results$w              # Donor weights [w1,...,wJ] 
result$est.results$constant_term  # Constant coefficient c
result$est.results$r              # Additional variables (includes constant)
```

## Interpretation

- **Positive constant**: Synthetic control has higher baseline level than donors
- **Negative constant**: Synthetic control has lower baseline level than donors  
- **Zero constant**: No level shift between synthetic control and donors

## Examples

### Basic Usage
```r
# 1. Prepare data with multiple features and constant
data_const <- scdata(df, ..., constant = TRUE, covagg = list(...))

# 2. Estimate with preferred constraint
result <- scest(data_const, w.constr = list(name = "simplex"))

# 3. Extract results
print(paste("Donor weights:", paste(result$est.results$w, collapse = ", ")))
print(paste("Constant term:", result$est.results$constant_term))
print(paste("Sum of weights:", sum(result$est.results$w)))
```

### Comparison Without/With Constant
```r
# Without constant
data_no_const <- scdata(df, ..., constant = FALSE)
result_no_const <- scest(data_no_const, w.constr = list(name = "simplex"))

# With constant  
data_const <- scdata(df, ..., constant = TRUE, covagg = list(...))
result_const <- scest(data_const, w.constr = list(name = "simplex"))

# Compare fit quality
cat("RMSE without constant:", sqrt(mean((result_no_const$est.results$res)^2)))
cat("RMSE with constant:", sqrt(mean((result_const$est.results$res)^2)))
```

## Technical Implementation Details

### Files Modified
- **`SCMs/R/scdata.R`**: Matrix construction (lines 847-853, 873-879)
- **`SCMs/R/scest.R`**: Result extraction (lines 258-263, 286)

### Key Changes
1. **Matrix Construction**: When `constant=TRUE`, adds constant rows to matrices A and B
2. **Result Processing**: Extracts constant term from optimization vector `b[J+1]`
3. **Infrastructure**: Leverages existing `J + KM` variable framework

### Solver Compatibility
- All existing solvers automatically support constant terms
- No changes needed to optimization problem formulation
- CVXR, clarabel, ECOS all work seamlessly

## Notes and Limitations

1. **Feature Requirement**: Must have multiple features to use constant term
2. **Identifiability**: Single feature + constant creates over-parameterized system  
3. **Constraint Logic**: Only donor weights constrained in simplex; constant always free
4. **Performance**: Minimal computational overhead from additional variable

## Migration from SCPI

If migrating from SCPI package:
- Same `constant = TRUE/FALSE` parameter in `scdata()`
- Similar mathematical formulation and results
- Compatible constraint specifications
- Results accessible via `$est.results$constant_term`