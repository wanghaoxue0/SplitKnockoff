# R function for other uses
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)
# two functions are included





#' hitting point calculator on a given path
#'
#' calculate the hitting time and the sign of
#' respective variable in a path.
#'
#' @param coef the path for one variable
#' @param lambda_vec respective value of lambda in the path
#'
#' @return Z: the hitting time
#' @return r: the sign of respective variable at the hitting time
#' @usage
#' hittingpoint(coef, lambda_vec)
#' @export
#'
hittingpoint <-function(coef, lambda_vec)
# input argument
# coef: the path for one variable
# lambdas: respective value of lambda in the path

# output argument
# Z: the hitting time
# r: the sign of respective variable at the hitting time
{

  n_lambda <- length(lambda_vec)

  Z <- 0
  r <- 0

  # calculate Z and r
  for (j in 1: n_lambda) {
      if(abs(coef[j]) != 0){
        Z = lambda_vec[j]
        r = sign(coef[j])
        break
      }
  }
  structure(list(call = match.call(),
                 Z = Z,
                 r = r),
            class = 'hittingpoint.result')
}


#' default normalization function for matrix
#'
#' normalize columns of a matrix.
#'
#' @param X the input martix
#'
#' @return Y the output matrix
#' @examples
#' library(mvtnorm)
#' n = 350
#' p = 100
#' Sigma = matrix(0, p, p)
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)
#' X <- normc(X)
#' @export
#'
normc <- function(X){
  n = nrow(X)
  p = ncol(X)
  X  <- X  - colMeans(X)
  factors <- 1 / sqrt(colSums(X^2))
  factors <- array(factors,dim=c(1,length(factors)))
  factors <- rep(factors,n)
  factors <- array(factors,dim=c(n,p))
# Y = X %*% factors(matrix(1,1,n),)    # This is used in Knockoffs.
  Y = X * factors* sqrt(n-1)
}
