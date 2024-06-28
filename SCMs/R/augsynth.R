get_ridgeAugSCM_weights <- function(scm_model,
                                    synth_obj,
                                    lambda=NULL,
                                    lambda_min_ratio = 1e-8,
                                    n_lambda = 20,
                                    lambda_max = NULL,
                                    holdout_length = 1,
                                    min_1se = T,
                                    ...) {
  
  extra_params = list(...)
  if (length(extra_params) > 0) {
    warning("Unused parameters in using ridge augmented weights: ", paste(names(extra_params), collapse = ", "))
  }
  
  # will convert so that rows are units, columns are pre-period characteristics
  X_c = t(synth_obj$B) # pre-period characteristics, donors
  X_1 = t(synth_obj$A) # pre-period characteristics, treated
  Z_c = t(synth_obj$Y.donors) # pre-period outcomes, donors
  Z_1 = t(synth_obj$Y.pre) # pre-period outcomes, treated
  scm_weights = scm_model$est.results$w
  
  
  lambda_errors <- NULL
  lambda_errors_se <- NULL
  lambdas <- NULL
  
  
  cv_out <- cv_lambda(Z_c,
                      Z_1,
                      scm_weights,
                      holdout_length,
                      lambda_max,
                      lambda_min_ratio,
                      n_lambda,
                      min_1se)
  
  lambda <- cv_out$lambda
  lambda_errors <- cv_out$lambda_errors
  lambda_errors_se <- cv_out$lambda_errors_se
  lambdas <- cv_out$lambdas
  
  # get ridge weights
  ridge_w <- t(t(X_1) - t(X_c) %*% scm_weights) %*%
    solve(t(X_c) %*% X_c  + lambda * diag(ncol(X_c))) %*% t(X_c)
  
  
  ## combine weights
  weights <- scm_weights + t(ridge_w)
  
  
  # ridge_weights <- weights
  # lambda <- out$lambda
  # lambdas <- out$lambdas
  # lambda_errors <- out$lambda_errors
  # lambda_errors_se <- out$lambda_errors_se
  
  
  # l2_imbalance <- sqrt(sum((synth_data$X0 %*% weights - synth_data$X1)^2))
  #
  # ## primal objective value scaled by least squares difference for mean
  # uni_w <- matrix(1/ncol(synth_data$X0), nrow=ncol(synth_data$X0), ncol=1)
  # unif_l2_imbalance <- sqrt(sum((synth_data$X0 %*% uni_w - synth_data$X1)^2))
  # scaled_l2_imabalance <- l2_imbalance / unif_l2_imbalance
  
  
  output <- list(ridge_weights = ridge_w,
                 ridge_scm_weights = weights,
                 lambda = lambda
  )
  
  
  return(output)
}



#' Choose best lambda with CV
#' @param X_c Matrix of control lagged outcomes
#' @param X_1 Vector of treated leagged outcomes
#' @param scm_weights Synth SCM weights
#' @param holdout_length Length of conseuctive holdout period for when tuning lambdas
#' @param lambda_max Initial (largest) lambda, if NULL sets it to be (1+norm(X_1-X_c))^2
#' @param lambda_min_ratio Ratio of the smallest to largest lambda when tuning lambda values
#' @param n_lambda Number of lambdas to consider between the smallest and largest lambda value
#' @param min_1se If TRUE, chooses the maximum lambda within 1 standard error of the lambda
#' @noRd
#' @return \itemize{
#'          \item{"lambda"}{Value of the ridge hyperparameter}
#'          \item{"lambdas"}{List of lambda values evaluated to tune ridge regression}
#'          \item{"lambda_errors"}{"The MSE associated with each lambda term in lambdas."}
#'          \item{"lambda_errors_se"}{"The SE of the MSE associated with each lambda term}
#' }

cv_lambda <- function(X_c,
                      X_1,
                      scm_weights,
                      holdout_length,
                      lambda_max,
                      lambda_min_ratio,
                      n_lambda,
                      min_1se) {
  if(is.null(lambda_max)) {
    lambda_max <- get_lambda_max(X_c)
  }
  
  lambdas <- create_lambda_list(lambda_max, lambda_min_ratio, n_lambda)
  
  lambda_out <- get_lambda_errors(lambdas,
                                  X_c,
                                  X_1,
                                  scm_weights,
                                  holdout_length)
  
  lambda_errors <- lambda_out$lambda_errors
  lambda_errors_se <- lambda_out$lambda_errors_se
  
  lambda <- choose_lambda(lambdas, lambda_errors, lambda_errors_se, min_1se)
  
  return(list(lambda = lambda, lambda_errors = lambda_errors,
              lambda_errors_se = lambda_errors_se, lambdas = lambdas))
}


#' Choose max lambda as largest eigenvalue of control X
#' @param X_c matrix of control lagged outcomes
#' @noRd
#' @return max lambda
get_lambda_max <- function(X_c) {
  svd(X_c)$d[1] ^ 2
}

#' Create list of lambdas
#' @param lambda_min_ratio Ratio of the smallest to largest lambda when tuning lambda values
#' @param n_lambda Number of lambdas to consider between the smallest and largest lambda value
#' @param lambda_max Initial (largest) lambda, if NULL sets it to be (1+norm(X_1-X_c))^2
#' @noRd
#' @return List of lambdas
create_lambda_list <- function(lambda_max,
                               lambda_min_ratio,
                               n_lambda) {
  scaler <- (lambda_min_ratio) ^ (1/n_lambda)
  lambdas <- lambda_max * (scaler ^ (seq(0:n_lambda) - 1))
  return(lambdas)
}

#' Choose either the lambda that minimizes CV MSE or largest lambda within 1 se of min
#' @param lambdas list of lambdas
#' @param lambda_errors The MSE associated with each lambda term in lambdas.
#' @param lambda_errors_se The SE of the MSE associated with each lambda
#' @param min_1se If TRUE, chooses the maximum lambda within 1 standard error of the lambda that minimizes the CV error, if FALSE chooses the optimal lambda; default TRUE
#' @noRd
#' @return optimal lambda
choose_lambda <- function(lambdas,
                          lambda_errors,
                          lambda_errors_se,
                          min_1se) {
  # lambda with smallest error
  min_idx <- which.min(lambda_errors)
  min_error <- lambda_errors[min_idx]
  min_se <- lambda_errors_se[min_idx]
  lambda_min <- lambdas[min_idx]
  # max lambda with error within one se of min
  lambda_1se <- max(lambdas[lambda_errors <= min_error + min_se])
  return(if(min_1se) lambda_1se else lambda_min)
}




################################################################################
## Function to calculate error on different lambda values if using Ridge Augmented SCM
################################################################################

#' Get Lambda Errors
#' @importFrom stats sd
#'
#' @param lambdas Vector of lambda values to compute errors for
#' @param X_c Matrix of control group pre-treatment outcomes
#' @param X_t Matrix of treatment group pre-treatment outcomes
#' @param scm_weights Synth SCM weights
#' @param holdout_length Length of conseuctive holdout period for when tuning lambdas
#' @noRd
#' @return List of lambda errors for each corresponding lambda in the lambdas parameter.
get_lambda_errors <- function(lambdas,
                              X_c,
                              X_t,
                              scm_weights,
                              holdout_length=1
) {
  
  # vector that stores the sum MSE across all CV sets
  errors <- matrix(0, nrow = ncol(X_c) - holdout_length, ncol = length(lambdas))
  lambda_errors = numeric(length(lambdas))
  lambda_errors_se = numeric(length(lambdas))
  
  for (i in 1:(ncol(X_c) - holdout_length)) {
    X_0 <- X_c[,-(i:(i + holdout_length - 1))]
    if(is.null(dim(X_0))){
      X_0 = matrix(X_0, ncol = 1)
    }
    
    X_1 <- matrix(X_t[-(i:(i + holdout_length - 1))])
    
    X_0v <- X_c[,i:(i + holdout_length - 1)]
    if(is.null(dim(X_0v))){
      X_0v = matrix(X_0v, ncol = 1)
    }
    
    X_1v <- matrix(X_t[i:(i + holdout_length - 1)], ncol = 1)
    
    syn = copy(scm_weights)
    
    for (j in 1:length(lambdas)) {
      ridge_weights <- t(X_1 - t(X_0) %*% syn) %*% solve(t(X_0) %*% X_0 + lambdas[j] * diag(ncol(X_0))) %*% t(X_0)
      aug_weights <- syn + t(ridge_weights)
      error <- X_1v - t(X_0v) %*% aug_weights
      # take sum of errors across the holdout time periods
      error <- sum(error ^ 2)
      errors[i, j] <-  error
      # lambda_errors[j] <- lambda_errors[j] + error
    }
  }
  
  lambda_errors <- apply(errors, 2, mean)
  lambda_errors_se <- apply(errors, 2, function(x) sd(x) / sqrt(length(x)))
  if(all(is.na(lambda_errors_se))){
    lambda_errors_se = rep(.0001, length(lambda_errors_se))
  }
  return(list(lambda_errors = lambda_errors,
              lambda_errors_se = lambda_errors_se))
}


