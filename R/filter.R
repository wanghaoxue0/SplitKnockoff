# R function for filter of structural sparsity problem
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)

#' Split Knockoff filter for structural sparsity
#'
#' the main function, Split Knockoff filter, for variable selection in structural sparsity problem.
#'
#' @param X the design matrix
#' @param D the linear transform
#' @param y the response vector
#' @param option various options for split knockoff filter, the details will be specified in the example
#'
#' @return results: a cell with the selected variable set in each cell w.r.t. nu.
#' @return Z: a cell with the feature significance Z in each cell w.r.t. nu.
#' @return t_Z: a cell with the knockoff significance tilde_Z in each cell w.r.t. nu.
#' @examples

#' option <- list(data = NA, dim = length(data), dimnames = NULL)
#' 
#' # the target (directional) FDR control
#' option$q <- 0.2    
#' 
#' # choice on threshold, the other choice is 'knockoff+'          
#' option$method <- 'knockoff' 
#' 
#' # degree of separation between original design and its split knockoff copy 
#' # in the range of [0, 2], the less the more separated
#' option$eta <- 0.1          
#'   
#' # whether to normalize the dataset
#' option$normalize <- 'true'
#' 
#' # choice on the set of regularization parameters for split LASSO path
#' option$lambda <- 10.^seq(0, -6, by=-0.01)
#' 
#' # choice of nu for split knockoffs
#' option$nu <- 10
#' 
#' # choice on whether to estimate the directional effect, 'disabled'/'enabled'
#' option$sign <- 'enabled'
#' 
#' option <- option[-1]

#'
#' # Settings on simulation parameters

#' k <- 20   # sparsity level
#' A <- 1    # magnitude
#' n <- 350  # sample size
#' p <- 100  # dimension of variables
#' c <- 0.5  # feature correlation
#' sigma <-1 # noise level

#' # generate D
#' D <- diag(p)
#' m <- nrow(D)

#' # generate X
#' Sigma = matrix(0, p, p)
#' for( i in 1: p){
#'   for(j in 1: p){
#'     Sigma[i, j] <- c^(abs(i - j))
#'   }
#' }

#' library(mvtnorm)
#' set.seed(100)
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)

#' # generate beta and gamma
#' beta_true <- matrix(0, p, 1)
#' for( i in 1: k){
#'   beta_true[i, 1] = A
#'   if ( i%%3 == 1){
#'     beta_true[i, 1] = -A
#'   }
#' }
#' gamma_true <- D %*% beta_true

#' S0 <- which(gamma_true!=0)


#' # generate varepsilon
#' set.seed(1)

#' # generate noise and y
#' varepsilon <- rnorm(n) * sqrt(sigma)
#' y <- X %*% beta_true + varepsilon
#' filter_result <- sk.filter(X, D, y, option)
#' Z_path <- filter_result$Z
#' t_Z_path <- filter_result$t_Z
#' @export
#'
sk.filter <- function(X, D, y, option)
{
  if(option$normalize == "true"){
    X <- normc(X) # normalize(X)
    y <- normc(y) # normalize(y)
  }
  nu_s = option$nu
  n_nu <- length(nu_s)
  Z <- list(data=NA,dim=c(1,n_nu))
  t_Z <- list(data=NA,dim=c(1,n_nu))
  r <- list(data=NA,dim=c(1,n_nu))

  n<- nrow(X)
  q = option$q
  method = option$method


  results <- list(data=NA,dim=c(1,n_nu))
  # if(all.equal(option$stage0, 'fixed') && all.equal(option$beta, 'ridge'))
  # {
  #   split_knockoffs.statistics.pathorder.fixed_beta(X, y, D, option)
  #   option$beta_choice <- beta_choice
  #
  # }
  for (i in 1: n_nu) {
    nu = nu_s[i]
    # filter.choose(option$stage0)
    if (all.equal(option$sign, 'disabled') == TRUE){
    path_result <- W_support(X, D, y, nu, option)
    W <- path_result$W
    Z[[i]] <- path_result$Z
    t_Z[[i]] <- path_result$t_Z
    select_results <- sk.select(W, q, option)
    result <- select_results$S
    results[[i]]<-result
    }
    else{
      path_result <- W_sign(X, D, y, nu, option)
      W <- path_result$W
      r_temp <- path_result$r
      Z[[i]] <- path_result$Z
      t_Z[[i]] <- path_result$t_Z
      select_results <- sk.select(W, q, option)
      result <- select_results$S
      results[[i]]<-result
      r[[i]] <- r_temp[result]
    }
  }
  if (all.equal(option$sign, 'disabled') == TRUE){
  structure(list(call = match.call(),
                 results = results,
                 W = W,
                 Z = Z,
                 t_Z = t_Z),
            class = 'filter_result')}
  else {
    structure(list(call = match.call(),
                   results = results,
                   W = W,
                   r = r,
                   Z = Z,
                   t_Z = t_Z),
              class = 'filter_result')}
  }








