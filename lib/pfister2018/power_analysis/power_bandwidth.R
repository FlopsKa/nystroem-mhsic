# Copyright (c) 2016
# Jonas Peters  [peters@stat.math.ethz.ch]
# Niklas Pfister [pfister@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.


# load library
library(dHSIC)
library(MASS)
library(parallel)

# source R functions
source("../extra/sampleDataFromG.R")
source("../extra/randomB.R")
source("../extra/computeCausOrder.R")
source("../extra/computeGaussKernel.R")

# load DAG file
load("allFullDagsWith4Nodes.RData")
numDAGs <- nrow(allFullDags)

# Define dependence functions
quadratic <- function(x){
  return(x^2)
}
linear <- function(x){
  return(x)
}
none <- function(x){
  return(0)
}
sinus <- function(x){
  return(sin(x))
}

# Median bandwidth
median_bandwidth <- function(x){
  len <- length(x)
  xnorm <- as.matrix(dist(x,method="euclidean",diag=TRUE,upper=TRUE))
  if(len>1000){
    sam <- sample(1:len,1000)
    xhilf <- xnorm[sam,sam]
  }
  else{
    xhilf <- xnorm
  }
  bandwidth <- sqrt(0.5*median(xhilf[lower.tri(xhilf,diag=FALSE)]^2))
  if(bandwidth==0){
      bandwidth <- 0.001
  }
  return(bandwidth)
}



# Define model
generate_data <- function(n,d,dependence){
  noise <- matrix(rnorm((d+1)*n),n,d+1)
  X <- vector("list",d)
  H <- noise[,1]
  for(i in 1:d){
    X[[i]] <- dependence(H)+2*noise[,i+1]
  }
  return(X)
}
generate_data2 <- function(n,d,dependence){
  noise <- matrix(rnorm(7*n),n,7)
  X <- vector("list",3)
  H1 <- noise[,1]
  H2 <- noise[,2]
  H3 <- noise[,7]
  X[[1]] <- H1+noise[,3]
  X[[2]] <- H1^2+noise[,4]
  X[[3]] <- H2+noise[,5]
  X[[4]] <- H2^2+noise[,6]
  return(X)
}
generate_data3 <- function(n,d,dependence){
  # sample two different random DAGs
  tmp <- sample(1:numDAGs,2,replace=FALSE)
  Gtrue <- matrix(allFullDags[tmp[1],],d,d)
  Gfalse <- matrix(allFullDags[tmp[1],],d,d)
  # generate data from causal model (with graph Gtrue)
  data <- sampleDataFromG(n,Gtrue,parsNoise=list(noiseExp=1,varMin=50,varMax=100,noiseExpVarMin=2,noiseExpVarMax=4,bound=rep(1,dim(Gtrue)[2])))
  data <- list(X1=data[,1],X2=data[,2],X3=data[,3],X4=data[,4])
  # initialize residual vector
  X <- vector("list",d)
  # regress each node on its parents and collect the residuals in X
  for(j in 1:d){
    if(sum(Gfalse[,j])==0){
      X[[j]] <- data[[j]]
    }
    else{
      formula <- paste("X",toString(j),"~",sep="")
      for(k in 1:d){
        if(Gfalse[k,j]==TRUE){
          if((sum(Gfalse[,j])-sum(Gfalse[1:k,j]))==0){
            formula <- paste(formula,"X",toString(k),sep="")
          }
          else{
            formula <- paste(formula,"X",toString(k),"+",sep="")
          }
        }
      }
      X[[j]] <- lm(as.formula(formula),data=data)$residuals
    }
  }
  return(X)
}
generate_data4 <- function(n,d,dependence){
  noise <- matrix(rnorm(4*n),n,4)
  noise <- noise/sqrt(rowSums(noise^2))
  X <- vector("list",4)
  X[[1]] <- noise[,1]
  X[[2]] <- noise[,2]
  X[[3]] <- noise[,3]
  X[[4]] <- noise[,4]
  return(X)
}
generate_data5 <- function(n,d,dependence){
  # sample a random DAGs
  tmp <- sample(1:numDAGs,2,replace=FALSE)
  G <- matrix(allFullDags[tmp[1],],d,d)
  # generate data from causal model (with graph Gtrue)
  data <- sampleDataFromG(n,G,parsNoise=list(noiseExp=1,varMin=50,varMax=100,noiseExpVarMin=2,noiseExpVarMax=4,bound=rep(1,dim(G)[2])))
  X <- list(X1=data[,1],X2=data[,2],X3=data[,3],X4=data[,4])
  return(X)
}
generate_data6 <- function(n,d,dependence){
  Gtrue <- matrix(c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),4,4)
  # generate data from causal model (with graph Gtrue)
  data <- sampleDataFromG(n,Gtrue,parsNoise=list(noiseExp=1,varMin=50,varMax=100,noiseExpVarMin=2,noiseExpVarMax=4,bound=rep(1,dim(Gtrue)[2])))
  X <- list(X1=data[,1],X2=data[,2],X3=data[,3],X4=data[,4])
  return(X)
}

generate_data7 <- function(n,d,dependence){
  X <- matrix(rnorm(3*n,0,1),n,3)
  ind <- (rowSums(X<0)==3)|(rowSums(X<0)==1)
  sum.ind <- sum(ind)
  while(sum.ind>0){    
    X[ind,] <- matrix(rnorm(3*sum.ind,0,1),sum.ind,3)
    ind <- (rowSums(X<0)==3)|(rowSums(X<0)==1)
    sum.ind <- sum(ind)
  }
  X <- list(X1=X[,1],X2=X[,2],X3=X[,3])
  return(X)
}

emp_power <- function(B,bandwidth,n,d,dependence,gen_data){
  # emprical dHSIC
  pval <- matrix(NA,B,length(bandwidth))
  pval2 <- vector("numeric",B)
  median <- matrix(NA,B,d)
  for(i in 1:B){
    if(i%%5==0){
      print(i)
    }
    X <- gen_data(n,d,dependence)
    for(k in 1:d){
      median[i,k] <- median_bandwidth(X[[k]])
    }
    pval2[i] <- dhsic.test(X,kernel="gaussian",method="permutation",B=100)$p.value
    for(j in 1:length(bandwidth)){
      pval[i,j] <- dhsic.test(X,kernel="gaussian.fixed",method="permutation",B=100,bandwidth=bandwidth[j])$p.value
    }
  }
  power <- colSums(pval<0.05)/B
  power2 <- sum(pval2<0.05)/B
  return(list(power=power,power2=power2,median=median))
}



## Run power analysis
n <- 100
d <- 4
B <- rep(50,20)
dependence <- linear

## Linear dependence example
bandwidth <- c(seq(0.1,3,0.1),seq(3.2,10,0.2))

tmp <- mclapply(B,emp_power,bandwidth,n,d,dependence,generate_data,mc.cores=20)
rejection.rate <- 0
rejection.rate2 <- 0
median <- numeric(0)
for(i in 1:length(tmp)){
  rejection.rate <- rejection.rate+tmp[[i]]$power
  rejection.rate2 <- rejection.rate2+tmp[[i]]$power2
  median <- rbind(median,tmp[[i]]$median)
}
rejection.rate <- rejection.rate/length(B)
rejection.rate2 <- rejection.rate2/length(B)

result <- list(rejection.rate=rejection.rate,
               rejection.rate2=rejection.rate2,
               median=median)

save(result,file="power.bandwidth.lin.rda")


## Complicated depedence example
bandwidth <- c(seq(0.15,1,0.05),seq(1.1,1.5,0.1),seq(1.6,20,0.5))

tmp <-  mclapply(B,emp_power,bandwidth,n,d,dependence,generate_data5,mc.cores=20)
rejection.rate <- 0
rejection.rate2 <- 0
median <- numeric(0)
for(i in 1:length(tmp)){
  rejection.rate <- rejection.rate+tmp[[i]]$power
  rejection.rate2 <- rejection.rate2+tmp[[i]]$power2
  median <- rbind(median,tmp[[i]]$median)
}
rejection.rate <- rejection.rate/length(B)
rejection.rate2 <- rejection.rate2/length(B)

result <- list(rejection.rate=rejection.rate,
               rejection.rate2=rejection.rate2,
               median=median)
save(result,file="power.bandwidth.uc.rda")


## Pairwise independent example
d <- 3
bandwidth <- c(seq(0.01,0.15,0.01),seq(0.2,10,0.1))

tmp <-  mclapply(B,emp_power,bandwidth,n,d,dependence,generate_data7,mc.cores=20)
rejection.rate <- 0
rejection.rate2 <- 0
median <- numeric(0)
for(i in 1:length(tmp)){
  rejection.rate <- rejection.rate+tmp[[i]]$power
  rejection.rate2 <- rejection.rate2+tmp[[i]]$power2
  median <- rbind(median,tmp[[i]]$median)
}
rejection.rate <- rejection.rate/length(B)
rejection.rate2 <- rejection.rate2/length(B)

result <- list(rejection.rate=rejection.rate,
               rejection.rate2=rejection.rate2,
               median=median)
save(result,file="power.bandwidth.pi.rda")
