# R function for giving the variable splitting design matrix
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)


#' generate split knockoff copies
#'
#' Give the variable splitting design matrix
#' and response vector. It will also create a split
#' knockoff copy if required.
#'
#' @param X the design matrix
#' @param y the response vector
#' @param D the linear transform
#' @param nu the parameter for variable splitting
#' @param option options for creating the Knockoff copy;
#' option$copy true : create a knockoff copy;
#' option$eta the choice of eta for creating the split knockoff copy
#'
#' @return A_beta: the design matrix for beta after variable splitting
#' @return A_gamma: the design matrix for gamma after variable splitting
#' @return tilde_y: the response vector after variable splitting.
#' @return tilde_A_gamma: the knockoff copy of A_beta; will be [] if option$copy = false.
#'
#' @examples
#' option <- array(data = NA, dim = length(data), dimnames = NULL)
#' option$q <- 0.2
#' option$eta <- 0.1
#' option$method <- 'knockoff'
#' option$normalize <- 'true'
#' option$lambda <- 10.^seq(0, -6, by=-0.01)
#' option$nu <- 10
#' option$copy <- 'true'
#' option$sign <- 'enabled'
#' option <- option[-1]
#' library(mvtnorm)
#' sigma <-1
#' p <- 100
#' D <- diag(p)
#' m <- nrow(D)
#' n <- 350
#' nu = 10
#' c = 0.5
#' Sigma = matrix(0, p, p)
#' for( i in 1: p){
#'   for(j in 1: p){
#'     Sigma[i, j] <- c^(abs(i - j))
#'  }
#' }
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)
#' beta_true <- matrix(0, p, 1)
#' varepsilon <- rnorm(n) * sqrt(sigma)
#' y <- X %*% beta_true + varepsilon
#' creat.result  <- sk.create(X, y, D, nu, option)
#' A_beta  <- creat.result$A_beta
#' A_gamma <- creat.result$A_gamma
#' tilde_y <- creat.result$tilde_y
#' tilde_A_gamma <- creat.result$tilde_A_gamma
#'
#'
#'
#' @export
#'
#'
sk.create <- function(X, y, D, nu, option)
# Input Argument:
# X : the design matrix.
# y : the response vector.
# D : the linear transform.
# nu: the parameter for variable splitting.
# option: options for creating the Knockoff copy.
#	option$copy = true : create a knockoff copy.
#	option$eta : the choice of eta for creating the knockoff copy.

# Output Argument:
# A_beta: the design matrix for beta after variable splitting.
# A_gamma: the design matrix for gamma after variable splitting.
# tilde_y: the response vector after variable splitting.
# tilde_A_gamma: the knockoff copy of A_beta; will be [] if option$copy = false.
{
  n <- nrow(X)
  m <- nrow(D)

  # calculate A_beta, A_gamma
  A_beta <- rbind(X/sqrt(n),D/sqrt(nu))
  A_gamma <- rbind(matrix(0,n,m),-diag(m)/sqrt(nu))

  # calculate tifdr[is.na(fdr)] <- 0lde_y
  tilde_y <- rbind(y/sqrt(n),matrix(0,m,1))

    s_size <- 2-option$eta

    # calculte inverse for Sigma_{beta, beta}
    Sigma_bb <- t(A_beta) %*% A_beta

    if(sum(abs(eigen(Sigma_bb)$values)>1e-6) == nrow(Sigma_bb)){
      Sigma_bb_inv = solve(Sigma_bb)
    }else{Sigma_bb_inv <- ginv(Sigma_bb)}

    # calculate Sigma_{gamma, gamma}, etc
    svd.result = canonicalSVD(A_gamma)
    S <- svd.result$S
    S <- diag(S)
    V <- svd.result$V
    Sigma_gg = (V %*% Matrix(S^2)) %*% t(V)
    Sigma_gg = as.matrix(Sigma_gg)
    Sigma_gb = t(A_gamma) %*% A_beta
    Sigma_bg = t(Sigma_gb)

    # calculate C_nu
    C = Sigma_gg - Sigma_gb %*% Sigma_bb_inv %*% Sigma_bg
    C = (C + t(C))/2
    C_inv = solve(C)
    t <- min(s_size * min(eigen(C)$val), 1/nu)
    t <- matrix(1,nrow(C),1)* t
    t <- array(t)
    # generate s
    diag_s <- diag(t)


    # calculate K^T K = 2S-S C_nu^{-1} S
    KK = 2 * diag_s - diag_s %*% C_inv %*% diag_s
    KK = (KK + t(KK))/2
    See <- diag(eigen(KK)$val)
    Uee <- eigen(KK)$vec
    K = Uee %*% sqrt(See) %*% t(Uee)

    # calculate U=[U1;U_2] where U_2 = 0_m* m
    # U_1 is an orthogonal complement of X
    decompose.result <- sk.decompose(X)
    U_perp <- decompose.result$U_perp
    U_1 = U_perp[,1:m]
    U = rbind(U_1, matrix(0,m,m))

    # calculate sigma_beta beta^{-1} sigma_beta gamma
    short = Sigma_bb_inv %*% Sigma_bg


    # calculate tilde_A_gamma
    tilde_A_gamma = A_gamma %*% (diag(m) - C_inv %*% diag_s) + A_beta %*%short%*%C_inv%*%diag_s+ U %*% K
  structure(list(call = match.call(),
                 A_beta = A_beta,
                 A_gamma = A_gamma,
                 tilde_y = tilde_y,
                 tilde_A_gamma =tilde_A_gamma),
            class = 'creat.result')

}
