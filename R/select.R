# R function for calculating the statistics W
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)

#' split knockoff selector given W statistics
#'
#' @param W statistics W_j for testing null hypothesis
#' @param q target FDR
#' @param option option$method can be 'knockoff' or 'knockoff+'
#'
#' @return S array of selected variable indices
#' @usage sk.select(W, q, option)
#' @export
#'
sk.select<-function(W, q, option){

  #  Inputs:
  #     W - statistics W_j for testing null hypothesis beta_j = 0.
  #       q - target FDR
  #
  #   Outputs:
  #       S - array of selected variable indices
  plus = 0
  if (all.equal(option$method, 'knockoff+') == TRUE){
    plus = 1
  }
  W <- t(W)
  t = sort(c(0,abs(W[W!=0])))
  ratio = matrix(0, 1, length(t))
  for (i in 1:length(t)){
  ratio[i] = (plus + sum(W <= -t[i])) / max(1, sum(W >= t[i]))
  }
  nindex <- which(ratio <= q)
  if(is.null(nindex)==TRUE){
    T <- Inf
  }else{
  index <- nindex[1]
  T = t[index]
  }
  S = which(W >= T)
  structure(list(call = match.call(),
                 S = S),
            class = 'selectS_result')
}
