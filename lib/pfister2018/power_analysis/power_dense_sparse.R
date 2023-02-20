# Copyright (c) 2016
# Jonas Peters  [peters@stat.math.ethz.ch]
# Niklas Pfister [pfister@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 


# load library
library(dHSIC)
library(MASS)
library(parallel)


# Define nonlinearity and dimension
nonlinear <- function(x){
  return(x^2)
}


# Define dense model
dense_setting <- function(n,d,nonlinear,c){
  noise <- matrix(rnorm((d+1)*n),n,d+1)
  X <- matrix(NA,n,d)
  H <- noise[,1]
  for(i in 1:d){
    X[,i] <- nonlinear(H)+c*noise[,i+1]
  }
  return(X)
}

# Define sparse model
sparse_setting <- function(n,d,nonlinear,c){
  noise <- matrix(rnorm((d+1)*n),n,d+1)
  X <- matrix(NA,n,d)
  H <- noise[,1]
  X[,1] <- nonlinear(H)+c*noise[,2]
  X[,2] <- nonlinear(H)+c*noise[,3]
  for(i in 3:d){
      X[,i] <- noise[,i+1]
  }
  return(X)
}

##
# total variation dense model
##
diff_dense <- function(c,r){
  d <- length(r)
  # product distribution
  dist.H <- function(x){
    return(dchisq(x,1))
  }
  dist.X <- function(x){
    return(1/c*dnorm(x/c))
  }
  dist.fun <- function(x){
    return(integrate(function(z) dist.X(x-z)*dist.H(z),-Inf,Inf)$value)
  }
  dist.fun <- Vectorize(dist.fun)
  fun1 <- function(y){
    return(integrate(dist.fun,-Inf,y)$value)
  }
  prod.dist <- prod(sapply(r,fun1))
  
  # joint distribution
  fun2 <- function(h){
    return(prod(sapply(((r-h)/c),pnorm))*dchisq(h,1))
  }
  fun2 <- Vectorize(fun2)
  joint.dist <- integrate(fun2,-Inf,Inf)$value
  # return difference
  return(abs(prod.dist-joint.dist))
}

##
# total variation sparse model
##
diff_sparse <- function(c,r){
  d <- length(r)
  # explicitely simplified convolution of the density of a chi-squared and a normal distribution
  conv.fun <- function(x){
    return(1/(sqrt(pi)*gamma(0.5)*c)*integrate(function(z) ((z-x)/c^2+0.5)*sqrt(z)*exp(-z/2)*exp(-((z-x)/c)^2/2),0,Inf)$value)
  }
  conv.fun <- Vectorize(conv.fun)
  fun1 <- function(y){
    return(integrate(conv.fun,-Inf,y)$value)
  }
  prod.dist <- prod(c(sapply(r[1:2],fun1),pnorm(r[3:d])))
  
  # joint distribution
  fun2 <- function(h){
    return(prod(sapply(((r[1:2]-h)/c),pnorm))*dchisq(h,1))
  }
  fun2 <- Vectorize(fun2)
  joint.dist <- prod(c(integrate(fun2,-Inf,Inf)$value,pnorm(r[3:d])))
  # return difference
  return(abs(prod.dist-joint.dist))
}


total_variation <- function(c,d,dense,C){
  sigma <- matrix(2,d,d)
  diag(sigma) <- 2+c^2
  mean <- rep(1,d)
  sets <- mvrnorm(C,mean,sigma)
  # Calculate distance between joint distribution and product distribution
  dist <- 0
  for(i in 1:C){
    r <- sets[i,]
    if(dense){
      dist.tmp <- diff_dense(c,r)
    }
    else{
      dist.tmp <- diff_sparse(c,r)
    }
    dist <- max(dist.tmp,dist)
  }
  return(dist)
}



calc_dist <- function(c,pars){
  totvar <- vector("numeric",2)
  # compute population total variation
  ptm1 <- proc.time()
  totvar[1] <- total_variation(c[1],pars$d,dense=TRUE,pars$C)
  totvar[2] <- total_variation(c[2],pars$d,dense=FALSE,pars$C)
  print(paste("totalvariation",round(c[1],2),":"))
  print(proc.time()-ptm1)
  # emprical dHSIC
  power <- vector("numeric",2)
  pval.m1 <- vector("numeric",pars$B)
  pval.m2 <- vector("numeric",pars$B)
  ptm3 <- proc.time()
  for(i in 1:pars$B){
    ## if(i%%10==0){
    ##   print(i)
    ## }
    X.m1 <- dense_setting(pars$n.emp,pars$d,pars$nonlinear,c[1])
    X.m1 <- split(t(X.m1),rep(1:pars$d,pars$n.emp))
    X.m2 <- sparse_setting(pars$n.emp,pars$d,pars$nonlinear,c[2])
    X.m2 <- split(t(X.m2),rep(1:pars$d,pars$n.emp))
    pval.m1[i] <- dhsic.test(X.m1,kernel="gaussian",method="permutation",B=100,bandwidth=1)$p.value
    pval.m2[i] <- dhsic.test(X.m2,kernel="gaussian",method="permutation",B=100,bandwidth=1)$p.value
  }
  print(paste("emprical dHSIC",round(c[1],2),":"))
  print(proc.time()-ptm3)
  power[1] <- sum(pval.m1<0.05)/pars$B
  power[2] <- sum(pval.m2<0.05)/pars$B
  return(list(totvar=totvar,power=power))
}



##
# d=5
##
d <- 5
c.dense <- seq(0.9,3,length.out=20)
c.sparse <- seq(0.1,2,length.out=20)
# Set parameters to correct format
c.tot <- cbind(c.dense,c.sparse)
c.tot <- split(c.tot,rep(1:nrow(c.tot),ncol(c.tot)))
# Initialize variables
n.emp <- 100
C <- 10000
B <- 10000
pars <- list(n.emp=n.emp,C=C,d=d,nonlinear=nonlinear,B=B)
# Computation
result <- mclapply(c.tot,calc_dist,pars,mc.cores=20)
save(result,file="power_comparison.d5.Rda")


##
# d=10
##
d <- 10
c.dense <- seq(2,4,length.out=20)
c.sparse <- seq(0.2,2,length.out=20)
# Set parameters to correct format
c.tot <- cbind(c.dense,c.sparse)
c.tot <- split(c.tot,rep(1:nrow(c.tot),ncol(c.tot)))
# Initialize variables
n.emp <- 100
C <- 10000
B <- 10000
pars <- list(n.emp=n.emp,C=C,d=d,nonlinear=nonlinear,B=B)
# Computation
result <- mclapply(c.tot,calc_dist,pars,mc.cores=20)
save(result,file="power_comparison.d10.Rda")
