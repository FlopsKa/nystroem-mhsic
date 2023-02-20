# Copyright (c) 2016  Niklas Pfister  [pfister@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

romano.test <- function(X, alpha=0.05, method="bootstrap", B=100, C=100){

  ###
  # Prerequisites
  ###

  # convert list to matrix
  X <- do.call(cbind,X)

  # determine dimensions
  d <- ncol(X)
  len <- nrow(X)
  
  ###
  # Sample intervals
  ###

  mean <- colSums(X)/len
  sigma <- cov(X)
  sets <- mvrnorm(C,mean,sigma)

  
  ###
  # Calculate test statistic
  ###

  test.stat <- 0
  for(i in 1:C){
    tmp <- t(t(X)<sets[i,])
    emp.joint <- sum(rowSums(tmp)==d)/len
    emp.prod <- prod(colSums(tmp)/len)
    test.stat <- max(abs(emp.joint-emp.prod),test.stat)
  }
  test.stat <- sqrt(len)*test.stat

  ###
  # Bootstrap test
  ###
  
  if(method=="bootstrap"){
    romano_boot_fun <- function(l){
      X_boot <- matrix(0,len,d)
      X_boot[,1] <- X[,1]
      for (j in 2:d){
        X_boot[,j] <- matrix(X[sample(len,replace=TRUE),j],nrow=len)
      }
      # calculate test.stat
      test.stat.boot <- 0
      for(i in 1:C){
        tmp <- t(t(X_boot)<sets[i,])
        emp.joint <- sum(rowSums(tmp)==d)/len
        emp.prod <- prod(colSums(tmp)/len)
        test.stat.boot <- max(abs(emp.joint-emp.prod),test.stat.boot)
      }
      return(sqrt(len)*test.stat.boot)
    }
    # Perform bootstrapping
    romano_boot <- sapply(1:B,romano_boot_fun)
    # Compute critical value and p-value
    sortromano <- sort(romano_boot)
    Bind <- sum(test.stat==sortromano)+ceiling((1-alpha)*(B+1))
    if(Bind<=B){
      critical_value <- sortromano[Bind]
    }
    else{
      critical_value <- Inf
    }
    p_value <- (sum(romano_boot>=test.stat)+1)/(B+1)
  }

  ###
  # Permutation test
  ###
  
  else if(method=="permutation"){
    romano_perm_fun <- function(l){
      X_perm <- matrix(0,len,d)
      X_perm[,1] <- X[,1]
      for (j in 2:d){
        X_perm[,j] <- matrix(X[sample(len,replace=FALSE),j],nrow=len)
      }
      # calculate test.stat
      test.stat.perm <- 0
      for(i in 1:C){
        tmp <- t(t(X_perm)<sets[i,])
        emp.joint <- sum(rowSums(tmp)==d)/len
        emp.prod <- prod(colSums(tmp)/len)
        test.stat.perm <- max(abs(emp.joint-emp.prod),test.stat.perm)
      }
      return(sqrt(len)*test.stat.perm)
    }
    # Perform permutation
    romano_perm <- sapply(1:B,romano_perm_fun)
    # Compute critical value and p-value
    sortromano <- sort(romano_perm)
    Bind <- sum(test.stat==sortromano)+ceiling((1-alpha)*(B+1))
    if(Bind<=B){
      critical_value <- sortromano[Bind]
    }
    else{
      critical_value <- Inf
    }
    p_value <- (sum(romano_perm>=test.stat)+1)/(B+1)
  }
    
  
  ###
  # Collect result
  ###
  test=list(statistic=test.stat,
            crit.value = critical_value,
            p.value = p_value)
  
  return(test)

}
