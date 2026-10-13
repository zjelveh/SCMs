# SCMs 0.1.1

Resubmission of 0.1.0, which CRAN returned for changes in February 2026.

## Bug fixes

* Covariates observed in only some pre-treatment years are aggregated over the observed years.
  Previously `mean()` received the missing values, the feature became `NA`, and it was dropped
  without a message, so a specification meant to match on many covariates could keep only one.
* Placebo inference keeps control units whose names contain spaces, hyphens or periods.
  Previously these units were silently left out of the placebo distribution.
* Two-sided placebo p-values (`expected_direction = "two_sided"`) rank units by absolute effect.
  Previously the unit with the largest effect received p = 0.
* Outcome models no longer share one mutable synthetic-control fit (the augmented SCM model
  altered the fit used by the others), and ridge/lasso cross-validation folds are deterministic.
* `run_spec_curve_analysis()` passes `cores` to `spec_curve()`, which previously ran
  single-threaded unless `params` contained its own `cores` entry.
* Parallel workers run the same SCMs code as the calling session (fork clusters on Unix); on
  Windows, a mismatch with the installed package is an error.

## Silent failures are now errors

* A failed placebo run stops with the unit's name instead of being skipped.
* `scdata()` names donors that have no pre-treatment data instead of failing with an unrelated
  error.
* Ridge penalty estimation stops instead of falling back to default values.
* `plot_spec_curve()` stops when SHAP values do not cover the specification grid.

## Changes

* `expected_direction` defaults to `"two_sided"` in `inference_sc()`, `inference_placebo()`,
  `spec_curve()` and `run_spec_curve_analysis()`.
* New curve-level statistic `curve_stat = "median_tau"` (signed median). One-sided curve-level
  tests (`two_sided = FALSE`) require it, because the other statistics are non-negative.
* `scest()` exposes the inner interior-point solver tolerances (`margin.ipop`, `sigf.ipop`,
  `bound.ipop`).
* SHAP attribution encodes every level of each categorical specification choice.
* `plot_spec_curve()` always draws Panel B labels as rich text (`richtext_feature_labels` is
  removed) and no longer prints SHAP merge diagnostics.
* Requires 'CVXR' 1.9.0 or later, which replaced `solve()` with `psolve()`. 'ECOSolveR' is now
  imported directly, because 'CVXR' no longer provides the default ECOS solver.
