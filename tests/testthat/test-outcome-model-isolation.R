# Regression tests for the run_outcome_models in-place mutation bug.
#
# Until 2026-08-28, the augsynth branch of run_outcome_models() overwrote
# scm_model$est.results$Y.post.fit and $Y.pre.fit in place. Three silent
# consequences followed, none of which any existing test caught:
#   1. the 'none' entry was assigned from those fields at the end of the
#      function, so "no outcome model" actually returned the Ridge ASCM fit;
#   2. the two-step ols/ridge/lasso corrections were applied on top of the
#      Ridge ASCM fit rather than the uncorrected fit;
#   3. every model inherited augsynth's pre-period fit, making pre_rmse
#      identical across all five models and feeding the pre-fit-weighted
#      placebo inference from the wrong baseline.
#
# These tests fail loudly if any of that returns.

make_panel <- function(n_units = 12, n_per = 10, treat_at = 7, seed = 42) {
  set.seed(seed)
  units <- paste0("u", seq_len(n_units))
  dt <- data.table::CJ(unit = units, period = seq_len(n_per))
  dt[, y := 10 + as.integer(sub("u", "", unit)) + period + stats::rnorm(.N, sd = 0.4)]
  dt[, x1 := 5 + 0.5 * as.integer(sub("u", "", unit)) + stats::rnorm(.N, sd = 0.3)]
  dt[, treat_period := treat_at]
  dt[]
}

fit_all <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dt <- make_panel()
  estimate_sc(
    dataset = dt, outcome = "y",
    covagg = list(y_mean = list(var = "y"), x1_mean = list(var = "x1")),
    col_name_unit_name = "unit", name_treated_unit = "u1",
    col_name_period = "period", treated_period = 7,
    min_period = 1, end_period = 10,
    outcome_models = c("none", "ols", "ridge", "lasso", "augsynth"),
    feature_weights = "uniform", w.constr = list(name = "simplex"),
    constant = FALSE
  )
}

test_that("'none' is the uncorrected fit, not the augsynth fit", {
  om <- fit_all()$est.results$outcome_model
  expect_false(isTRUE(all.equal(
    as.numeric(om$none$Y.post.fit), as.numeric(om$augsynth$Y.post.fit)
  )), info = "'none' must not equal 'augsynth' -- the in-place mutation is back")
})

test_that("every outcome model carries its own Y.pre.fit", {
  om <- fit_all()$est.results$outcome_model
  for (m in c("none", "ols", "ridge", "lasso", "augsynth")) {
    expect_false(is.null(om[[m]]$Y.pre.fit),
                 info = paste("missing per-model Y.pre.fit for", m))
    expect_gt(length(as.numeric(om[[m]]$Y.pre.fit)), 0)
  }
})

test_that("pre_rmse is not shared across all five outcome models", {
  r  <- fit_all()
  om <- r$est.results$outcome_model
  ypre <- as.numeric(r$data$Y.pre)
  rmse <- vapply(c("none", "ols", "ridge", "lasso", "augsynth"), function(m) {
    sqrt(mean((ypre - as.numeric(om[[m]]$Y.pre.fit))^2))
  }, numeric(1))
  expect_gt(length(unique(round(rmse, 10))), 1,
            label = "distinct pre_rmse values across outcome models")
})

test_that("the augsynth branch does not mutate the shared SCM fit", {
  with_aug <- fit_all()$est.results
  dt <- make_panel()
  without_aug <- estimate_sc(
    dataset = dt, outcome = "y",
    covagg = list(y_mean = list(var = "y"), x1_mean = list(var = "x1")),
    col_name_unit_name = "unit", name_treated_unit = "u1",
    col_name_period = "period", treated_period = 7,
    min_period = 1, end_period = 10,
    outcome_models = "none", feature_weights = "uniform",
    w.constr = list(name = "simplex"), constant = FALSE
  )$est.results
  # Requesting augsynth must not change the base fit, nor the 'none' entry.
  expect_equal(as.numeric(with_aug$Y.post.fit), as.numeric(without_aug$Y.post.fit))
  expect_equal(as.numeric(with_aug$outcome_model$none$Y.post.fit),
               as.numeric(without_aug$outcome_model$none$Y.post.fit))
})

test_that("ridge and lasso outcome models are reproducible across seeds", {
  a <- fit_all(seed = 1)$est.results$outcome_model
  b <- fit_all(seed = 999)$est.results$outcome_model
  for (m in c("ridge", "lasso")) {
    expect_equal(as.numeric(a[[m]]$Y.post.fit), as.numeric(b[[m]]$Y.post.fit),
                 info = paste(m, "must not depend on the RNG seed"))
  }
})

test_that("deterministic_foldid is stable and restores the RNG stream", {
  f1 <- SCMs:::deterministic_foldid(40)
  f2 <- SCMs:::deterministic_foldid(40)
  expect_identical(f1, f2)
  expect_equal(sort(unique(f1)), 1:10)
  set.seed(123); before <- stats::runif(1)
  set.seed(123); invisible(SCMs:::deterministic_foldid(40)); after <- stats::runif(1)
  expect_equal(before, after, info = "must not disturb the caller's RNG stream")
})
