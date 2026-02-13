# CRAN Resubmission Punch List (`SCMs` 0.1.0)

Source: CRAN feedback email from Benjamin Altmann (Feb 12, 2026) on submission `SCMs_0.1.0.tar.gz`.

## Scope Summary
- `DESCRIPTION` policy formatting: **small**
- References added to `DESCRIPTION`: **small**
- Example/vignette execution policy (`\dontrun{}` / commented examples): **medium**
- Console output cleanup (`cat()`/`print()` in non-interactive core paths): **medium-large**
- Global assignment (`<<-`) removal: **critical, localized**

## Action Punch List

| ID | CRAN request | Scope | Target files (initial) | Action | Status |
|---|---|---|---|---|---|
| 1 | Quote software/package/API names in single quotes in `Title`/`Description` | Small | `DESCRIPTION`, `R/package.R` (roxygen source), regenerated `man/SCMs-package.Rd` | Use quoted names, correct case (`'xgboost'`, `'SCMs'`, etc.), keep wording CRAN-style | DONE |
| 2 | Add method references in `Description` field using CRAN format (`authors (year) <doi:...>`) | Small | `DESCRIPTION` | Add 2-3 key references with DOI/URL in approved syntax | DONE |
| 3 | Vignettes should execute meaningful code | Medium | `vignettes/*.Rmd` | Ensure each vignette runs executable chunks (not just narrative) under CRAN-friendly timing | DONE |
| 4 | Replace unnecessary `\dontrun{}` with `\donttest{}`; unwrap very fast examples | Medium | `R/*.R` roxygen blocks; regenerated `man/*.Rd` | Convert wrappers based on runtime; keep `\dontrun{}` only when truly required | DONE |
| 5 | Remove commented-out example code (`inference_sc.Rd`) | Small | `R/inference_sc.R` roxygen source, `man/inference_sc.Rd` | Replace commented examples with executable toy example(s) + `\donttest{}` if needed | DONE |
| 6 | Avoid unsuppressible console output (`cat`/`print`) in package functions | Medium-large | `R/memory_management.R`, `R/scplot.R`, plus other `R/*.R` with `cat(` | Replace with `message()`/`warning()` or `if (verbose)`-guarded output; keep print methods as print methods | DONE |
| 7 | Do not write to global env (`<<-`) | Critical | `R/memory_management.R` | Remove/replace `<<-` with local assignment and explicit returns/state updates | DONE |
| 8 | Rebuild docs and run checks before resubmission | Medium | full package | `devtools::document()`, `R CMD build`, `R CMD check --as-cran`, optionally win-builder (`R-devel`, `R-release`) | DONE |

## Validation Checklist Before Resubmit
- [x] `DESCRIPTION` updated (quoted software names + references).
- [x] `Roxygen` regenerated (`man/` and `NAMESPACE` refreshed).
- [x] `rg "\\\\dontrun\\{"` reviewed and justified.
- [x] `rg "<<-" R` returns no policy-violating global writes.
- [x] `rg -n "^[[:space:]]*cat\\(|^[[:space:]]*print\\(" R` reviewed for non-interactive core functions.
- [x] `R CMD check --as-cran` now at `1 WARNING, 4 NOTE` (no ERRORs).
- [x] Win-builder checks reviewed (`R-devel` pretest NOTE only).

Current residual check items are environment/tooling only:
- WARNING: `qpdf` not installed (PDF size reduction check).
- NOTE: CRAN incoming `New submission`.
- NOTE: local time verification (`unable to verify current time`) on this machine.
- NOTE: local `libxml` build mismatch warning (`compiled against libxml 210 using older 209`).
- NOTE: HTML manual tools unavailable (`tidy`, `V8`).

## Draft Resubmission Reply (Template)

Thank you for the review. We have addressed all points:

1. Updated `DESCRIPTION` formatting to use single-quoted software/API names and corrected case-sensitive package naming.
2. Added method references in CRAN-preferred format (`authors (year) <doi:...>` / `<https:...>`).
3. Revised vignette/example structure to include executable code paths suitable for user testing and CRAN checks.
4. Replaced unnecessary `\\dontrun{}` blocks with `\\donttest{}` (and unwrapped very short examples).
5. Removed commented-out example code (including `inference_sc` examples) and replaced with executable examples.
6. Refactored console output in package functions to avoid unsuppressible `cat()/print()` in non-interactive code paths.
7. Removed global assignment (`<<-`) patterns.

Checks run before resubmission:
- `R CMD check --as-cran` locally (Linux)
- win-builder checks for `R-devel` and `R-release`
