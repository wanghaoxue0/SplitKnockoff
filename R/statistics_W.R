# R function for calculating the statistics W
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)



#' W statistics generator for FDR control
#'
#' generate the split knockoff statistics W for a split LASSO 
#' path, only consider the support set estimation
#'
#' @param X the design matrix
#' @param D the linear transform
#' @param y the response vector
#' @param nu the parameter for variable splitting
#' @param option options for creating the Split Knockoff statistics; 
#' option$eta specify the choice of eta for creating the knockoff copy; 
#' option$lambda specify the choice of lambda for the split LASSO path
#'
#' @return W: the knockoff statistics
#' @return Z: feature significance
#' @return t_Z: knockoff significance
#' 
#' @usage
#' W_support(X, D, y, nu, option)
#' @export
#'
#'
W_support <-function(X, D, y, nu, option)
  # input argument:
  # X : the design matrix
  # y : the response vector
  # D : the linear transform
  # nu: the parameter for variable splitting
  # option: options for creating the Knockoff statistics
  #	option$eta : the choice of eta for creating the knockoff copy
  #	option$lambda: the choice of lambda for the path

  # output argument
  # W: the knockoff statistics
  # Z: feature significance
  # t_Z: knockoff significance
{
  m <- nrow(D)
  p <- ncol(D)

  # generate the design matrix
  creat.result  <- sk.create(X, y, D, nu, option)
  A_beta  <- creat.result$A_beta
  A_gamma <- creat.result$A_gamma
  tilde_y <- creat.result$tilde_y
  tilde_A_gamma <- creat.result$tilde_A_gamma

  ############ step 0 #############

  # set lambda
  lambda_vec <- option$lambda
  nlambda <- length(lambda_vec)

  # set penalty
  penalty <- matrix(1,m+p,1)

  for (i in 1: p) {
    penalty[i, 1] = 0
  }

  # lasso path settings for glmnet

  fit_step0 = glmnet(cbind(A_beta,A_gamma) , tilde_y, lambda =lambda_vec, penalty.factor = penalty)
  coefs <- fit_step0$beta

  # store beta(lambda)
  betas = coefs[1: p, ]

  ############ step 1 #############
  coef1 = matrix(0,m, 1)
  for (i in 1: nlambda) {
    # take beta_lambda, gamma_lambda as calculated in step 1
    y_new <- tilde_y - A_beta %*% betas[, i]
    # calculate LASSO
    # opts = struct
    lambda = lambda_vec[i]
    fit_step1 = glmnet(A_gamma, y_new, lambda =lambda)
    # coef1[, i] <- fit_step1$beta
    coef <- fit_step1$beta
    coef1 <- cbind(coef1,coef)
  }
   coef1 <- coef1[,-1]
   # calculate r and Z
   r <- matrix(0,m,1)
   Z <- matrix(0,m,1)

   for (i in 1: m) {
     hit <- hittingpoint(coef1[i, ], lambda_vec)
     Z[i] <- hit$Z
     r[i] <- hit$r
   }



  ############# step 2 ############
   coef2 <- matrix(0,m,1)
   for (i in 1: nlambda) {
     # take beta_lambda, gamma_lambda as calculated in step 1
     y_new <- tilde_y - A_beta %*% betas[, i]
     # calculate LASSO
     # opts = struct
     lambda = lambda_vec[i]
     # opts = glmnetSet(opts)
     fit_step2 <- glmnet(tilde_A_gamma, y_new, lambda=lambda)
     coef <- fit_step2$beta
     coef2 <- cbind(coef2,coef)
   }
   coef2 <- coef2[,-1]
   # calculate tilde_Z tilde_r and W
   t_Z <- matrix(0,m, 1)
   t_r <- matrix(0,m, 1)

  for (i in 1: m) {
    t_hit <- hittingpoint(coef2[i, ], lambda_vec)
    t_r[i] <- t_hit$r
    if(t_r[i] == r[i]){
      # store tilde_Z when it has the same sign
      t_Z[i] = t_hit$Z
    }
  }


  ############ W ###########
   W <- matrix(0,m, 1)
   for (i in 1: m) {
   W[i] <- max(Z[i], t_Z[i]) * sign(Z[i] - t_Z[i])
   }
   structure(list(call = match.call(),
                  W = W,
                  Z = Z,
                  t_Z = t_Z),
             class = 'path_result')
}





#' W statistics generator for directional FDR control
#'
#' generate the split knockoff statistics W for a split LASSO 
#' path, take the directional effect into account
#'
#' @param X the design matrix
#' @param D the linear transform
#' @param y the response vector
#' @param nu the parameter for variable splitting
#' @param option options for creating the Knockoff statistics
#' option$eta specify the choice of eta for creating the knockoff copy; 
#' option$lambda specify the choice of lambda for the split LASSO path
#'
#' @return W: the knockoff statistics
#' @return Z: feature significance
#' @return r: the sign estimator
#' @return t_Z: knockoff significance
#' 
#' @usage
#' W_sign(X, D, y, nu, option)
#' @export
#' 

W_sign <-function(X, D, y, nu, option)
  # input argument:
  # X : the design matrix
  # y : the response vector
  # D : the linear transform
  # nu: the parameter for variable splitting
  # option: options for creating the Knockoff statistics
  #	option$eta : the choice of eta for creating the knockoff copy
  #	option$lambda: the choice of lambda for the path

  # output argument
  # W: the knockoff statistics
# Z: feature significance
# t_Z: knockoff significance
{
  m <- nrow(D)
  p <- ncol(D)

  # generate the design matrix
  creat.result  <- sk.create(X, y, D, nu, option)
  A_beta  <- creat.result$A_beta
  A_gamma <- creat.result$A_gamma
  tilde_y <- creat.result$tilde_y
  tilde_A_gamma <- creat.result$tilde_A_gamma

  ############ step 0 #############

  # set lambda
  lambda_vec <- option$lambda
  nlambda <- length(lambda_vec)

  # set penalty
  penalty <- matrix(1,m+p,1)

  for (i in 1: p) {
    penalty[i, 1] = 0
  }

  # lasso path settings for glmnet

  fit_step0 = glmnet(cbind(A_beta,A_gamma) , tilde_y, lambda =lambda_vec, penalty.factor = penalty)
  coefs <- fit_step0$beta

  # store beta(lambda)
  betas = coefs[1: p, ]

  ############ step 1 #############
  coef1 = coefs[(p+1):(p+m),]
  # calculate r and Z
  r <- matrix(0,m,1)
  Z <- matrix(0,m,1)

  for (i in 1: m) {
    hit <- hittingpoint(coef1[i, ], lambda_vec)
    Z[i] <- hit$Z
    r[i] <- hit$r
  }
  # for (i in 1: nlambda) {
  #   # take beta_lambda, gamma_lambda as calculated in step 1
  #   y_new <- tilde_y - A_beta %*% betas[, i]
  #   # calculate LASSO
  #   # opts = struct
  #   lambda = lambda_vec[i]
  #   fit_step1 = glmnet(A_gamma, y_new, lambda =lambda)
  #   # coef1[, i] <- fit_step1$beta
  #   coef <- fit_step1$beta
  #   coef1 <- cbind(coef1,coef)
  # }
  #  coef1 <- coef1[,-1]
  #
  #
  #  for (i in 1: m) {
  #    hit <- hittingpoint(coef1[i, ], lambda_vec)
  #    Z[i] <- hit$Z
  #    r[i] <- hit$r
  #  }
  #


  ############# step 2 ############
  coef2 <- matrix(0,m,1)
  for (i in 1: nlambda) {
    # take beta_lambda, gamma_lambda as calculated in step 1
    y_new <- tilde_y - A_beta %*% betas[, i]
    # calculate LASSO
    # opts = struct
    lambda = lambda_vec[i]
    # opts = glmnetSet(opts)
    fit_step2 <- glmnet(tilde_A_gamma, y_new, lambda=lambda)
    coef <- fit_step2$beta
    coef2 <- cbind(coef2,coef)
  }
  coef2 <- coef2[,-1]
  # calculate tilde_Z tilde_r and W
  t_Z <- matrix(0,m, 1)

  for (i in 1: m) {
    t_hit <- hittingpoint(coef2[i, ], lambda_vec)
    t_Z[i] <- t_hit$Z
  }


  ############ W ###########
  W <- matrix(0,m, 1)
  for (i in 1: m) {
    W[i] <- Z[i] * sign(Z[i] - t_Z[i])
  }
  structure(list(call = match.call(),
                 W = W,
                 r = r,
                 Z = Z,
                 t_Z = t_Z),
            class = 'path_result')
}

