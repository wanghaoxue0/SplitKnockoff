% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{sk.filter}
\alias{sk.filter}
\title{Split Knockoff filter for structural sparsity}
\usage{
sk.filter(X, D, y, option)
}
\arguments{
\item{X}{the design matrix}

\item{D}{the linear transform}

\item{y}{the response vector}

\item{option}{various options for split knockoff filter, the details will be specified in the example}
}
\value{
results: a cell with the selected variable set in each cell w.r.t. nu.

Z: a cell with the feature significance Z in each cell w.r.t. nu.

t_Z: a cell with the knockoff significance tilde_Z in each cell w.r.t. nu.
}
\description{
the main function, Split Knockoff filter, for variable selection in structural sparsity problem.
}
\examples{
option <- list(data = NA, dim = length(data), dimnames = NULL)

# the target (directional) FDR control
option$q <- 0.2    

# choice on threshold, the other choice is 'knockoff+'          
option$method <- 'knockoff' 

# degree of separation between original design and its split knockoff copy 
# in the range of [0, 2], the less the more separated
option$eta <- 0.1          
  
# whether to normalize the dataset
option$normalize <- 'true'

# choice on the set of regularization parameters for split LASSO path
option$lambda <- 10.^seq(0, -6, by=-0.01)

# choice of nu for split knockoffs
option$nu <- 10

# choice on whether to estimate the directional effect, 'disabled'/'enabled'
option$sign <- 'enabled'

option <- option[-1]

# Settings on simulation parameters
k <- 20   # sparsity level
A <- 1    # magnitude
n <- 350  # sample size
p <- 100  # dimension of variables
c <- 0.5  # feature correlation
sigma <-1 # noise level
# generate D
D <- diag(p)
m <- nrow(D)
# generate X
Sigma = matrix(0, p, p)
for( i in 1: p){
  for(j in 1: p){
    Sigma[i, j] <- c^(abs(i - j))
  }
}
library(mvtnorm)
set.seed(100)
X <- rmvnorm(n,matrix(0, p, 1), Sigma)
# generate beta and gamma
beta_true <- matrix(0, p, 1)
for( i in 1: k){
  beta_true[i, 1] = A
  if ( i\%\%3 == 1){
    beta_true[i, 1] = -A
  }
}
gamma_true <- D \%*\% beta_true
S0 <- which(gamma_true!=0)
# generate varepsilon
set.seed(1)
# generate noise and y
varepsilon <- rnorm(n) * sqrt(sigma)
y <- X \%*\% beta_true + varepsilon
filter_result <- sk.filter(X, D, y, option)
Z_path <- filter_result$Z
t_Z_path <- filter_result$t_Z
}
