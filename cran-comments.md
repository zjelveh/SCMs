## Resubmission

This is a resubmission of SCMs 0.1.0. Thank you for the review of 2026-02-12. Each request is
quoted below with our response.

This upload also fixes the arXiv reference flagged by the incoming pre-tests on 2026-10-13: it now
uses the arXiv DOI <doi:10.48550/arXiv.1705.07874>.

> Please always write package names, software names and API (application programming interface)
> names in single quotes in title and description. [...] Please note that package names are case
> sensitive.

Software names in the Description are now in single quotes with their exact case ('xgboost').

> If there are references describing the methods in your package, please add these in the
> description field of your DESCRIPTION file in the form authors (year) <doi:...>

Added references with DOIs for the synthetic control, augmented synthetic control,
specification-curve, 'xgboost' and SHAP methods.

> Some of your vignettes don't execute any code.

All three vignettes now execute code; a few long-running chunks are shown but not run.

> Please replace \dontrun with \donttest. Please unwrap the examples if they are executable in
> < 5 sec, or replace dontrun{} with \donttest{}.

\dontrun{} is no longer used. Examples that run in under 5 seconds are unwrapped; longer ones use
\donttest{}.

> Some code lines in examples are commented out. Please never do that. [...] Examples in comments
> in: inference_sc.Rd

inference_sc.Rd now has an executable toy example, wrapped in \donttest{}.

> You write information messages to the console that cannot be easily suppressed. [...]
> -> R/memory_management.R; R/scplot.R; etc...

Outside print and summary methods, functions no longer call cat() or print(); progress output
uses message(), which users can suppress. This includes R/memory_management.R and R/scplot.R.

> Please do not modify the global environment (e.g. by using <<-) in your functions. [...]
> -> R/memory_management.R

Removed: no function assigns with <<-.

The version is now 0.1.1. It also includes bug fixes made since the first submission, listed in
NEWS.md.

## Test environments

* Local: Ubuntu 22.04.4 LTS, R 4.4.1, with all dependencies at their current CRAN versions
* win-builder: R-devel (2026-10-10 r90659) and R-release (4.6.1)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.
* Possibly misspelled words in DESCRIPTION: these are author names in the method references and
  the acronym SHAP (SHapley Additive exPlanations).
