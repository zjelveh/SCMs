#' Deterministic cross-validation fold assignment
#'
#' @description Builds a reproducible `foldid` vector for [glmnet::cv.glmnet()].
#'   `cv.glmnet` otherwise draws folds from the session RNG, which made the ridge
#'   and lasso outcome models irreproducible run to run: the specification-curve
#'   pipeline sets no seed, so every run selected a different penalty and produced
#'   different bias-corrected estimates (differences up to ~50 in the Philadelphia
#'   grid). Supplying `foldid` removes the dependency entirely rather than relying
#'   on a caller to seed correctly.
#'
#'   The draw is derived from a fixed seed and the caller's RNG stream is restored
#'   afterwards, so this function is deterministic AND leaves global random state
#'   untouched (important because it is called inside parallel workers).
#'
#' @param n Integer. Number of observations.
#' @param nfolds Integer. Requested number of folds; capped at `n`.
#' @param seed Integer. Fixed seed for the fold draw.
#' @return Integer vector of length `n` with fold assignments.
#' @keywords internal
deterministic_foldid <- function(n, nfolds = 10L, seed = 20260206L) {
  nf <- max(2L, min(as.integer(nfolds), as.integer(n)))
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  sample(rep_len(seq_len(nf), n))
}

run_outcome_models = function(scm_model,
                              scm_data,
                              treated_unit,
                              outcome_models,
                              period_post,
                              Y,
                              Z0,
                              Z1){
  
  
  scm_model$est.results$outcome_model = list()

  # Capture the UNCORRECTED SCM fits before any outcome model runs.
  # Previously the augsynth branch below overwrote scm_model$est.results$Y.post.fit
  # and $Y.pre.fit in place, so (a) the 'none' entry, assigned from those fields at the
  # end of this function, silently returned the Ridge ASCM fit rather than the
  # uncorrected fit; (b) the two-step ols/ridge/lasso corrections were applied on top
  # of the Ridge ASCM fit instead of the base fit; and (c) every model inherited
  # augsynth's pre-period fit, making pre_rmse identical across all five models and
  # feeding the pre-fit-weighted placebo inference from the wrong baseline.
  base_Y_post_fit <- scm_model$est.results$Y.post.fit
  base_Y_pre_fit  <- scm_model$est.results$Y.pre.fit

  if('augsynth' %in% outcome_models){
    ridge_aug_model = get_ridgeAugSCM_weights(scm_model,
                                              scm_data,
                                              lambda=NULL,
                                              lambda_min_ratio = 1e-8,
                                              n_lambda = 20,
                                              lambda_max = NULL,
                                              holdout_length = 1,
                                              min_1se = T)
    temp1 = Y %*% ridge_aug_model$ridge_scm_weights

    scm_model$est.results$augsynth = ridge_aug_model

    temp2 = Y %*% ridge_aug_model$ridge_scm_weights
    scm_model$est.results$outcome_model$augsynth = list()
    scm_model$est.results$outcome_model$augsynth$Y.post.fit = temp2
    # augsynth uses its own donor weights, so it has its own pre-period fit.
    # Stored per-model; the shared scm_model$est.results$Y.pre.fit is left untouched.
    scm_model$est.results$outcome_model$augsynth$Y.pre.fit = Z0 %*% ridge_aug_model$ridge_scm_weights
  }
  
  outcome_models = tolower(outcome_models)
  # for outcome models (except augsynth)
  other_outcome_models = outcome_models[!outcome_models %in% c('none', 'augsynth')]
  
  Y_control = t(Y)
  X_control = t(Z0)
  X_treated = t(Z1)
  
  yhat_list = list()
  for(outcome_model in other_outcome_models){
    yhat_list[[outcome_model]] = list()
    yhat_list[[outcome_model]]$control = Y_control - Y_control
    yhat_list[[outcome_model]]$treated = (Y_control - Y_control)[1, ,drop=FALSE]
  }
  
  # train and predict
  for(period in period_post){
    period = paste0(treated_unit, '.' , period)
    
    y_control = Y_control[, period, drop=F]
    
    for(outcome_model in other_outcome_models){
      if(outcome_model%in%c('ols', 'lm')){
        ols_model <- glmnet::glmnet(X_control, y_control, alpha = 0, lambda = 0)
        coefs = coef(ols_model)
      }
      
      if(outcome_model%in%c('ridge')){
        ridge_model <- glmnet::cv.glmnet(X_control, y_control, alpha = 0, grouped=FALSE,
                                         foldid = deterministic_foldid(nrow(X_control)))
        coefs = coef(ridge_model)
      }
      
      if(outcome_model%in%c('lasso')){
        lasso_model <- glmnet::cv.glmnet(X_control, y_control, alpha = 1, grouped=FALSE,
                                         foldid = deterministic_foldid(nrow(X_control)))
        coefs = coef(lasso_model)
      }
      yhat_list[[outcome_model]]$control[, period] = (cbind(1, X_control) %*% coefs)[,1]
      yhat_list[[outcome_model]]$treated[, period] = (cbind(1, X_treated) %*% coefs)[1]  

    }

  }

  for(outcome_model in other_outcome_models){
    temp1 = t(yhat_list[[outcome_model]]$treated)
    temp2 = t(yhat_list[[outcome_model]]$control) %*% scm_model$est.results$w 
    scm_model$est.results$outcome_model[[outcome_model]] = list()
    scm_model$est.results$outcome_model[[outcome_model]]$Y.post.fit = base_Y_post_fit + (temp1 - temp2)
    # The two-step correction adjusts post-period predictions only and does not change
    # the donor weights, so the pre-period fit is the uncorrected SCM's.
    scm_model$est.results$outcome_model[[outcome_model]]$Y.pre.fit = base_Y_pre_fit
  }
    

  if(length(other_outcome_models)>0){
    scm_model$est.results$yhats = yhat_list
  }
  
  # create entry for no outcome model -- the uncorrected SCM fit
  scm_model$est.results$outcome_model$none = list()
  scm_model$est.results$outcome_model$none$Y.post.fit = base_Y_post_fit
  scm_model$est.results$outcome_model$none$Y.pre.fit  = base_Y_pre_fit

  return(scm_model)
}
