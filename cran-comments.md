## Resubmission

This is a resubmission of SCMs 0.1.0, reviewed on 2026-02-12. Thank you for the review. In
response:

1. `Title` and `Description` put software, package and API names in single quotes, with
   correct case.
2. `Description` cites the methods in the form authors (year) <doi:...>.
3. The vignettes execute their code.
4. `\dontrun{}` is no longer used. Examples that take more than a few seconds are wrapped in
   `\donttest{}`; fast examples are unwrapped.
5. The commented-out example code (`inference_sc`) is replaced by an executable example.
6. Functions other than print and summary methods no longer write to the console with
   `cat()` or `print()`; progress output uses `message()`, which users can suppress.
7. The package no longer assigns to the global environment (`<<-` removed).

The version is now 0.1.1. It also includes bug fixes made since the first submission, listed in
NEWS.md.

## Test environments

* Local: Ubuntu 22.04.4 LTS, R 4.4.1
* win-builder: R-devel and R-release (PENDING)
* macOS builder: R-release (PENDING)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

Locally, two further NOTEs come from the build machine only ("unable to verify current time",
and a libxml version message).
