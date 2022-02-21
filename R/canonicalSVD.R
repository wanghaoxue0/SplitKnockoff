# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)

# Reduced SVD with canonical sign choice

# Computes a reduced SVD without sign ambiguity. Our convention is that
# the sign of each vector in U is chosen such that the coefficient
# with largest absolute value is positive.

#' singular value decomposition
#'
#' Computes a reduced SVD without sign ambiguity
#'
#' @param X the input matrix
#'
#' @return S
#' @return U
#' @return V
#' @examples
#' nu = 10
#' n = 350
#' m = 100
#' A_gamma <- rbind(matrix(0,n,m),-diag(m)/sqrt(nu))
#' svd.result = canonicalSVD(A_gamma)
#' S <- svd.result$S
#' S <- diag(S)
#' V <- svd.result$V
#'
#' @export
#'
canonicalSVD <- function(X)
{
  svd<-svd(X)
  S <- svd$d
  U <- svd$u
  V <- svd$v
  for (j in min(dim(X))) {
    i <- which.max(abs(U[,j]))
    if (U[i,j] < 0){
    U[,j] = -U[,j]
    V[,j] = -V[,j]
    }
  }
  structure(list(call = match.call(),
                 S = S,
                 U = U,
                 V = V),
            class = 'svd_result')
}

