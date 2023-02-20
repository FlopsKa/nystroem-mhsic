# Copyright (c) 2016
# Jonas Peters  [peters@stat.math.ethz.ch]
# Niklas Pfister [pfisteni@student.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

## Level analysis


# load required libraries
library(dHSIC)
library(parallel)

# load required function
source("benchmarking.R")

# define generating function of H0 (3-variable standard normal) - Simulation 1
gaussianiid <- function(sample_size){
  return(list(X1=rnorm(sample_size),X2=rnorm(sample_size),X3=rnorm(sample_size)))
}

# define generating function of H0 (2-variable standard normal and binomial) - Simulation 2
gaussiandiscrete <- function(sample_size){
  return(list(X1=rnorm(sample_size),X2=rbinom(sample_size,20,0.2)))
}

# define generating function of H0 (d-variable standard normal) - Simulation 3
gaussianiid_d <- function(sample_size){
  d <- 10
  X <- matrix(rnorm(sample_size*d),sample_size,d)
  Xlist <- split(X, rep(1:ncol(X), each = nrow(X)))
  return(Xlist)
}

# perform benchmarking
n=1000

# Simulation 1
ptm <- proc.time()
set.seed(1)
B <- 25
sample_size=c(100,200,300,400,500,600,700,800,900,1000)
rejectionrates_level1 <- benchmarking(n,sample_size,gaussianiid,B,"gaussian",FALSE,20)
print(proc.time()-ptm)
save(rejectionrates_level1,file="rejection_rates_level1.Rda")

# Simulation 2
ptm <- proc.time()
set.seed(1)
B <- 100
sample_size=c(100,200,300,400,500,600,700,800,900,1000)
rejectionrates_level2 <- benchmarking(n,sample_size,gaussiandiscrete,B,c("gaussian","discrete"),FALSE,20)
print(proc.time()-ptm)
save(rejectionrates_level2,file="rejection_rates_level2.Rda")
