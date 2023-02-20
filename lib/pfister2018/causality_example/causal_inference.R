# Copyright (c) 2016
# Jonas Peters  [jonas.peters@math.ku.dk]
# Niklas Pfister [pfister@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

## Causal inference example with romano method


# load required libraries
library(dHSIC)
library(mgcv)
require(Matrix)
require(MASS)
library(parallel)
library(graph)
library(RBGL)



# required R scripts
source("../extra/romano.test.R")
source("../extra/randomDAG.R")
source("../extra/sampleDataFromG.R")
source("../extra/randomB.R")
source("../extra/computeCausOrder.R")
source("../extra/computeGaussKernel.R")
source("../extra/structIntervDist.R")
source("../extra/computePathMatrix.R")
source("../extra/computePathMatrix2.R")
source("../extra/dSepAdji.R")

# load DAG file
load("allFullDagsWith4Nodes.RData")
numDAGs <- nrow(allFullDags)
d <- 4


# define function that finds best DAG
findbestDAG <- function(X, G, method){
  # initalize pvalue vector
  pval_dhsic <- vector("numeric",numDAGs)
  pval_pairwise <- vector("numeric",numDAGs)
  pval_romano <- vector("numeric",numDAGs)
  
  # find the index of the correct DAG
  correctDAG <- which(rowSums(t(matrix(rep(as.numeric(G),numDAGs),d^2,numDAGs))==allFullDags)==d^2)
  
  # compute the p-value for each DAG
  for(i in 1:numDAGs){
    # generate Adj-matrix for DAG
    A <- matrix(allFullDags[i,],d,d)
    # initialize residual vector
    res <- vector("list",d)
    ress <- matrix(NA, length(X[[1]]), length(X))
    # regress each node on its variables and collect the residuals
    for(j in 1:d){
      if(sum(A[,j])==0){
        res[[j]] <- X[[j]]
        ress[,j] <- res[[j]]
      }
      else{
        formula <- paste("X",toString(j),"~",sep="")
        for(k in 1:d){
          if(A[k,j]==TRUE){
            if((sum(A[,j])-sum(A[1:k,j]))==0){
              formula <- paste(formula,"s(X",toString(k),")",sep="")
            }
            else{
              formula <- paste(formula,"s(X",toString(k),")+",sep="")
            }
          }
        }
        res[[j]] <- gam(as.formula(formula),data=X)$residuals
        ress[,j] <- res[[j]]
      }
    }
    pval_dhsic[i] <- dhsic.test(res,method=method,kernel="gaussian",pairwise=FALSE,B=100)$p.value
    pval_pairwise[i] <- dhsic.test(res,method=method,kernel="gaussian",pairwise=TRUE,B=100)$p.value
    pval_romano[i] <- romano.test(res,method="bootstrap",B=100,C=length(res[[1]]))$p.value
#    if (pval_dhsic[i] > 0.05)
#    {
#        pairs(ress)
#    }
  }
  # compute index of best fitting DAG
  bestfitdhsic <- which.max(pval_dhsic)
  bestfitpairwise <- which.max(pval_pairwise)
  bestfitromano <- which.max(pval_romano)
  
  # compute difference in fits
  dhsic_fit <- pval_dhsic[bestfitdhsic]-pval_dhsic[correctDAG]
  pairwise_fit <- pval_pairwise[bestfitpairwise]-pval_pairwise[correctDAG]
  romano_fit <- pval_romano[bestfitromano]-pval_romano[correctDAG]
  
  # compute structural intervention distance
  sid_dhsic <- structIntervDist(matrix(allFullDags[correctDAG,],d,d),matrix(allFullDags[bestfitdhsic,],d,d))$sid
  sid_pairwise <- structIntervDist(matrix(allFullDags[correctDAG,],d,d),matrix(allFullDags[bestfitpairwise,],d,d))$sid
  sid_romano <- structIntervDist(matrix(allFullDags[correctDAG,],d,d),matrix(allFullDags[bestfitromano,],d,d))$sid
  
#  diff <- sid_dhsic
  diff2 <- cbind(pval_dhsic[bestfitdhsic],
                 pval_dhsic[correctDAG],
                 correctDAG == bestfitdhsic,
                 as.numeric(sid_dhsic),
                 pval_pairwise[bestfitpairwise],
                 pval_pairwise[correctDAG],
                 correctDAG == bestfitpairwise,
                 as.numeric(sid_pairwise),
                 pval_romano[bestfitromano],
                 pval_romano[correctDAG],
                 correctDAG == bestfitromano,
                 as.numeric(sid_romano))
  
  return(diff2)
  
}


# define simulation function
simulate_causalex <- function(sample_size,n){
  result <- matrix(nrow=n,ncol=12)
  for(i in 1:n){
    G <- randomDAG(4,1)
    X <- sampleDataFromG(sample_size,G)
    X <- list(X1=X[,1],X2=X[,2],X3=X[,3], X4=X[,4])
    result[i,] <- findbestDAG(X, G, "gamma")
  }
  return(result)
}


# Simulation
n <- 10
sample_size <- c(100,200,300,400,500)

set.seed(1)
ptm <- proc.time()
results <- mclapply(sample_size,simulate_causalex, n=n, mc.cores=10)
print(proc.time()-ptm)

save(results,file='causal_inference.Rda')
