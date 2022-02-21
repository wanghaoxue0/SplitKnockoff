# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)

# This script reproduces Figure 1 of the paper, comparing the different
# distributions of Knockoff statistics and Split Knockoff statistics.

#root = getwd()
#setwd(sprintf('%s/R/', root))
#source("beta.r")
#source("create.r")
#source("cv.r")
#source("filter.r")
#source("magnitude.r")
#source("others.r")
#source("simu.r")
#source("statistics_W.r")
#source("canonicalSVD.R")
#source("decompose.R")
#source("select.R")
library(latex2exp)
library(ggplot2)
library(Matrix)
library(glmnet)
library(MASS)
library(SplitKnockoff)
library(testthat)


test_that('test splitknockoff.filter, the main function.', {

k <- 20   # sparsity level
A <- 1    # magnitude
n <- 350  # sample size
p <- 100  # dimension of variables
c <- 0.5  # feature correlation
sigma <-1 # noise level

option <- array(data = NA, dim = length(data), dimnames = NULL)
option$q <- 0.2
option$eta <- 0.1
option$method <- 'knockoff'
option$stage0 <- 'path'
option$normalize <- 'true'
option$cv_rule <- 'min'
option$lambda <- 10.^seq(0, -6, by=-0.01)
option$nu <- 10
option$copy <- 'true'
option$sign <- 'disabled'
option <- option[-1]


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

library(mvtnorm) # package mvtnorm needed for this generation
set.seed(100)
X <- rmvnorm(n,matrix(0, p, 1), Sigma) # generate X

# generate beta and gamma
beta_true <- matrix(0, p, 1)
for( i in 1: k){
  beta_true[i, 1] = A
  if ( i%%3 == 1){
  beta_true[i, 1] = -A
  }
}
gamma_true <- D %*% beta_true

S0 <- which(gamma_true!=0)


# generate varepsilon
set.seed(1)

# generate noise and y
varepsilon <- rnorm(n) * sqrt(sigma)
y <- X %*% beta_true + varepsilon



  ## Split Knockoff
filter_result <- sk.filter(X, D, y, option)
Z_path <- filter_result$Z
t_Z_path <- filter_result$t_Z

expect_equal(length(Z_path), 2)
expect_equal(length(t_Z_path), 2)
expect_equal(length(filter_result), 5)
})

