# Covariate Aggregation (`covagg`) Guide

`covagg` controls how pre-treatment matching features are generated.

## Core Operation Signature

For `scdata()`, `covagg` is a list of operations:

```r
covagg = list(
  list(
    var = "outcome_var",                              # or explicit column name
    select_periods = list(type = "all_pre"),         # optional
    partition_periods = list(type = "by_period"),    # optional
    compute = "mean"                                  # optional
  )
)
```

Each operation has four fields:
- `var`: variable name (`"outcome_var"` resolves to `outcome.var`).
- `select_periods`: chooses periods from the pre-window.
- `partition_periods`: groups selected periods into one or more features.
- `compute`: function name or function used within each group.

If `covagg` is empty, `scdata()` defaults to outcome per-period matching.

## Selectors

Supported `select_periods` forms:
- `NULL`: all pre-periods.
- `list(type = "all_pre")`
- `list(type = "explicit", periods = c(...))`
- `list(type = "every_n", n = 2, offset = 1)`
- `list(type = "predicate", fn = function(periods) ...)`
- `function(periods) ...` (returns logical mask or numeric periods)

## Partitions

Supported `partition_periods` forms:
- `NULL` or `list(type = "all")`: one group.
- `list(type = "by_period")`: one group per period.
- `list(type = "fixed_width", width = 3, step = 3)`
- `list(type = "custom", groups = list(c(...), c(...)))`
- `function(selected_periods) ...` (returns list of numeric vectors)

## Common Recipes

Outcome only (per period):

```r
covagg = list(
  list(var = "outcome_var", partition_periods = list(type = "by_period"))
)
```

All variables meaned:

```r
covagg = list(
  list(var = "outcome_var", compute = "mean"),
  list(var = "population", compute = "mean"),
  list(var = "income", compute = "mean")
)
```

Outcome per period + covariates meaned:

```r
covagg = list(
  list(var = "outcome_var", partition_periods = list(type = "by_period")),
  list(var = "population", compute = "mean"),
  list(var = "income", compute = "mean")
)
```

## `spec_curve()` Shape

For `spec_curve()`, `covagg` is a named list of specifications, each with `operations`:

```r
covagg = list(
  "Outcome Only" = list(
    label = "Outcome Only",
    operations = list(
      list(var = "outcome_var", partition_periods = list(type = "by_period"))
    )
  ),
  "Outcome + Means" = list(
    label = "Outcome + Means",
    operations = list(
      list(var = "outcome_var", partition_periods = list(type = "by_period")),
      list(var = "population", compute = "mean")
    )
  )
)
```

## Troubleshooting

`unsupported fields`:
- You passed removed keys (`per_period`, `pre_period_mean`, `each`, `first`, `last`, `average`, etc.).

`must define a non-empty operations list`:
- In `spec_curve()`, each top-level specification must include `operations = list(...)`.

`selected periods are outside pre-period window`:
- Your selector returned periods not in `period.pre`.
