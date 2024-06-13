run_outcome_models = function(scm_model,
                              scm_data,
                              treated_unit,
                              outcome_models,
                              period_post,
                              Y,
                              Z0,
                              Z1){
  
  library(glmnet)
  
  if('augsynth' %in% outcome_models){
    ridge_aug_model = get_ridgeAugSCM_weights(scm_model,
                                              scm_data,
                                              lambda=NULL,
                                              lambda_min_ratio = 1e-8,
                                              n_lambda = 20,
                                              lambda_max = NULL,
                                              holdout_length = 1,
                                              min_1se = T)
    temp1 = scm_data$Y.donors.post %*% ridge_aug_model$ridge_scm_weights
    
    scm_model$est.results$Y.post.fit[,1] = t(temp1)
    scm_model$est.results$augsynth = ridge_aug_model
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
  
  for(period in period_post){
    period = paste0(treated_unit, '.' , period)
    
    y_control = Y_control[, period, drop=F]
    
    for(outcome_model in other_outcome_models){
      if(outcome_model%in%c('ols', 'lm')){
        ols_model <- glmnet::glmnet(X_control, y_control, alpha = 0, lambda = 0)
        coefs = coef(ols_model)
        yhat_list[[outcome_model]]$control[, period] = (cbind(1, X_control) %*% coefs)[,1]
        yhat_list[[outcome_model]]$treated[, period] = (cbind(1, X_treated) %*% coefs)[1]  
      }
      
      if(outcome_model%in%c('ridge')){
        ridge_model <- glmnet::cv.glmnet(X_control, y_control, alpha = 0, grouped=FALSE)
        coefs = coef(ridge_model)
        yhat_list[[outcome_model]]$control[, period] = (cbind(1, X_control) %*% coefs)[,1]
        yhat_list[[outcome_model]]$treated[, period] = (cbind(1, X_treated) %*% coefs)[1]        
      }
      
      if(outcome_model%in%c('lasso')){
        lasso_model <- glmnet::cv.glmnet(X_control, y_control, alpha = 1, grouped=FALSE)
        coefs = coef(lasso_model)
        yhat_list[[outcome_model]]$control[, period] = (cbind(1, X_control) %*% coefs)[,1]
        yhat_list[[outcome_model]]$treated[, period] = (cbind(1, X_treated) %*% coefs)[1]
      }
    }

  }
  
for(outcome_model in other_outcome_models){
      if(outcome_model%in%c('ols', 'lm')){
        temp1 = t(yhat_list[[outcome_model]]$treated)
        temp2 = t(yhat_list[[outcome_model]]$control) %*% scm_model$est.results$w 
        scm_model$est.results$Y.post.fit[,1] = scm_model$est.results$Y.post.fit[,1] + (temp1 - temp2) 
      }
      
      if(outcome_model%in%c('ridge')){
        temp1 = t(yhat_list[[outcome_model]]$treated)
        temp2 = t(yhat_list[[outcome_model]]$control) %*% scm_model$est.results$w 
        scm_model$est.results$Y.post.fit[,1] = scm_model$est.results$Y.post.fit[,1] + (temp1 - temp2)
      }
      
      if(outcome_model%in%c('lasso')){
        temp1 = t(yhat_list[[outcome_model]]$treated)
        temp2 = t(yhat_list[[outcome_model]]$control) %*% scm_model$est.results$w 
        scm_model$est.results$Y.post.fit[,1] = scm_model$est.results$Y.post.fit[,1] + (temp1 - temp2)
      }
    }
  

  if(length(other_outcome_models)>0){
    scm_model$est.results$yhats = yhat_list
  }
  
  return(scm_model)
}
