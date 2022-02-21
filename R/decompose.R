# R function for filter of structural sparsity problem
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)

# KNOCKOFFS.PRIVATE.DECOMPOSE  Decompose design matrix X for knockoff creation

#' make SVD as well as orthogonal complements
#' 
#' 
#'
#' @param X the input matrix
#' @param randomize whether to randomize
#'
#' @return U
#' @return S
#' @return V
#' @return U_perp : orthogonal complement for U
#' @examples
#' library(mvtnorm)
#' n = 350
#' p = 100
#' Sigma = matrix(0, p, p)
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)
#' decompose.result <- sk.decompose(X)
#' U_perp <- decompose.result$U_perp
#' @export
#'
sk.decompose <- function(X, randomize){
  # Check dimensions.
  n <- nrow(X)
  p <- ncol(X)
  if(n < 2*p){
    print("knockoff:DimensionError: Data matrix must have n >= 2p")
  }
  # Factorize X as X = USV' (reduced SVD).
  svd.result <- canonicalSVD(X)
  S <- svd.result$S
  S <- diag(S)
  V <- svd.result$V
  U <- svd.result$U
  UU <- cbind(U,matrix(0,n,n-p))
  qrresult <- qr(UU)
  Qreslt <- qr.Q(qrresult)
  U_perp <- Qreslt[,(p+1):n]
  structure(list(call = match.call(),
                 U = U,
                 S = S,
                 V = V,
                 U_perp =U_perp),
            class = 'decompose.result')
}

