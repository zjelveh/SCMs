# Constant Term Support in SCMs

## Overview

`SCMs` supports a global constant (intercept) term in synthetic-control estimation.
This allows a level shift between treated and synthetic units while preserving the
main donor-weight structure.

## Current Implementation

When `constant = TRUE` is set in `scdata()`:

- `A` and `B` remain feature matrices.
- A separate constant design matrix `C` is created.
- Estimation uses the form:

`min ||A - Z w - C gamma||_V^2`

where:

- `w` are donor weights
- `gamma` is the constant coefficient

The constant is handled explicitly in the optimizer output and stored in
`est.results$constant_term`.

## Minimal Example

```r
library(SCMs)

scm_data <- scdata(
  df = panel_df,
  id.var = "unit",
  time.var = "year",
  outcome.var = "outcome",
  period.pre = 2000:2009,
  period.post = 2010:2015,
  unit.tr = "treated",
  unit.co = c("c1", "c2", "c3"),
  constant = TRUE,
  covagg = list(
    list(var = "outcome_var", partition_periods = list(type = "by_period")),
    list(var = "population", compute = "mean")
  )
)

fit <- scest(
  data = scm_data,
  w.constr = list(name = "simplex")
)

fit$est.results$w
fit$est.results$constant_term
```

## Constraint Behavior

- `simplex`: donor weights are simplex-constrained; constant is estimated separately.
- `ols`, `lasso`, `ridge`, `pensynth`: supported under the same estimation interface.

## Notes

- Use the operation-style `covagg` schema documented in `?scdata` and `COVAGG.md`.
- If `covagg` is empty, `scdata()` defaults to outcome-per-period matching.
- The authoritative implementation details are in `R/scdata.R` and `R/scest.R`.
