# Covariate Aggregation (`covagg`) Guide

This guide explains how `covagg` turns pre‑treatment data into matching features.
It is aligned with the current code, and is written for both human readers and
LLMs.

## 1) What `covagg` Controls (Mental Model)

`covagg` defines **how you construct matching variables** from the pre‑period.
Each `covagg` entry produces either:
- **One feature** (e.g., mean of GDP over pre‑periods), or
- **Many features** (one per pre‑period, when `each = TRUE`).

## 2) Where You Use It

- `spec_curve(...)`: you define a **set** of matching‑variable constructions.
- `scdata(...)`: you define a **single** matching‑variable construction.

## 3) Choose a Format

There are two supported formats, and **both work in `spec_curve(...)` and `scdata(...)`**.

**Use Format A (Simplified)** for clarity and grid design.  
**Use Format B (Traditional)** when you need fine‑grained control.

Do **not** mix the two formats inside a single `covagg` entry.

## Format A: Simplified Specs (Recommended for `spec_curve`)

Each specification is a list with optional `label`, plus:
- `per_period`: variables included **one per pre‑period**
- `pre_period_mean`: variables included as **one mean feature** over all pre‑periods

Use the sentinel `"outcome_var"` to refer to the outcome variable.

```r
covagg = list(
  "Outcome Only" = list(
    label = "Outcome Only",
    per_period = "outcome_var"
  ),
  "All Mean" = list(
    label = "All Mean",
    pre_period_mean = c("population", "income", "outcome_var")
  ),
  "Outcome Per Period + All Mean" = list(
    label = "Outcome Per Period + All Mean",
    per_period = "outcome_var",
    pre_period_mean = c("population", "income")
  )
)
```

## Format B: Traditional Specs (Detailed Control)

Each entry is **one** covariate aggregation rule and must include:
- `var` (required): a single variable name

Optional fields:
- `periods`: explicit period values to use
- `first`: first N pre‑periods
- `last`: last N pre‑periods
- `each`: if `TRUE`, create one feature per period
- `agg_fun`: aggregation function (default: `mean`)
- `label`: optional display label

```r
covagg = list(
  gdp_mean = list(var = "gdp"),                    # mean over all pre‑periods
  trade_last5 = list(var = "trade", last = 5),     # mean over last 5 pre‑periods
  gdp_each = list(var = "gdp", each = TRUE),       # one feature per pre‑period
  inv_med_early = list(var = "investment", first = 4, agg_fun = median),
  gdp_selected = list(var = "gdp", periods = c(1980, 1985), each = TRUE)
)
```

## 4) `scdata(...)` Notes (Important)

`scdata(...)` accepts **both** formats. When given simplified specs, it internally
converts them to traditional entries before aggregation.

Traditional fields supported by `scdata(...)` are:
- `var` (required)
- `periods`, `first`, or `last` (optional)
- `each` (optional)
- `agg_fun` or `aggfunc` (optional; default is `mean`)

`agg_fun`/`aggfunc` can be a function (recommended) or a function name string.
If the function supports `na.rm`, it will be used.

## 5) Rules (Keep These Straight)

- Use **only one** of `periods`, `first`, or `last` per spec.
- If none are provided, SCMs uses **all pre‑periods**.
- `each = TRUE` creates **one feature per period** (no aggregation).
- `agg_fun` must return a **single numeric** value.
- Provide `label` if you want clean display names in plots.

## 6) Not Supported (Built‑In)

Rolling windows, growth rates, and volatility are **not** built‑in.  
If you need them, **precompute** them and pass the derived variables via `var`.

## 7) Quick Troubleshooting

- **“Invalid covagg format”**: you mixed simplified and traditional fields in one spec.
- **“cannot have multiple period specifications”**: remove overlap between `periods`,
  `first`, and `last`.
- **Missing variable warning**: the variable is not present in your data or pre‑period.
