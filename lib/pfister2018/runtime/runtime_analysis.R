# Copyright (c) 2016
# Jonas Peters  [peters@stat.math.ethz.ch]
# Niklas Pfister [pfisteni@student.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 


#### Runtime Analysis

# load required libraries
library(dHSIC)
library(parallel)
library(dHSIC)


# data generating function
gendata <- function(m,d){
  X <- vector("list",d)
  for(i in 1:d){
    X[[i]] <- rnorm(m)
  }
  return(X)
}


# define runtime function
runtime <- function(m,d,n){
  time_dhsic <- vector("numeric",n)
  time_pairwise <- vector("numeric",n)
  for(i in 1:n){
    X <- gendata(m,d)
    # dhsic
    ptm1 <- proc.time()
    dhsic(X)
    time_dhsic[i] <- as.numeric((proc.time()-ptm1)[3])
    ptm2 <- proc.time()
    for(j in d:2){
      dhsic(do.call("cbind", X[1:(j - 1)]),X[[j]])
    }
    time_pairwise[i] <- as.numeric((proc.time() - ptm2)[3])
  }
  time <- list(time_dhsic,time_pairwise)
  return(time)
}

# fixed sample size
m <- 100
d <- seq(3,100,4)
n <- 100
req_vars <- list(d=d)

time.vars<- lapply(d,runtime,m=m,n=n)
dhsic <- vector("list",length(d))
pairwise <- vector("list",length(d))
for(i in 1:length(d)){
  dhsic[[i]] <- time.vars[[i]][[1]]
  pairwise[[i]] <- time.vars[[i]][[2]]
}
time.vars <- list(dhsic=dhsic,pairwise=pairwise)
save(time.vars,file="time.vars.Rda")

par(mfrow=c(1,1))
max_t <- max(unlist(time.vars$pairwise))
boxplot(time.vars$dhsic, border ="firebrick3",outline=FALSE,xaxt="n",ylim=c(0,max_t))
axis(1,d,at=1:length(d))
boxplot(time.vars$pairwise, border="chartreuse4", add=TRUE,outline=FALSE,xaxt="n",yaxt="n")


# fixed variable number
m <- seq(100,2000,100)
d <- 10
n <- 100
req_vars <- c(req_vars,list(m=m))
save(req_vars,file="req_vars.runtime.Rda")

time.sample<- lapply(m,runtime,d=d,n=n)
dhsic <- vector("list",length(d))
pairwise <- vector("list",length(d))
for(i in 1:length(m)){
  dhsic[[i]] <- time.sample[[i]][[1]]
  pairwise[[i]] <- time.sample[[i]][[2]]
  
}
time.sample <- list(dhsic=dhsic,pairwise=pairwise)
save(time.sample,file="time.sample.Rda")

max_t <- max(unlist(time.sample$pairwise))
boxplot(time.sample$dhsic, border ="firebrick3",outline=FALSE,xaxt="n",ylim=c(0,max_t))
axis(1,m,at=1:length(m))
boxplot(time.sample$pairwise, border="chartreuse4", add=TRUE,outline=FALSE,xaxt="n",yaxt="n")
