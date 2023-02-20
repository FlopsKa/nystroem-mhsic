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
load("../power_analysis/allFullDagsWith4Nodes.RData")
numDAGs <- nrow(allFullDags)
d <- 4


# define function that finds best DAG
generate_residuals <- function(X, G, run){
  
  # find the index of the correct DAG
  correctDAG <- which(rowSums(t(matrix(rep(as.numeric(G),numDAGs),d^2,numDAGs))==allFullDags)==d^2)
  
  # compute the residuals for each DAG
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
        print(formula)
        res[[j]] <- gam(as.formula(formula),data=X)$residuals
        ress[,j] <- res[[j]]
      }
    }
    p <- paste("../../../../data/causal_sim/", s	, "/", run, "/" , i,".csv",sep="")
    print(p)
    write.csv(res, p)
  

  }
  write.csv(correctDAG, paste("../../../../data/causal_sim/", s, "/", run, "/" , "correct.csv",sep="")) 				
 
  
}

n <- 100
sample_size <- c(500,1000,1500)

for(s in sample_size) {
print(s)
	for(i in 1:n){
	    G <- randomDAG(4,1)
	    X <- sampleDataFromG(s,G)
	    X <- list(X1=X[,1],X2=X[,2],X3=X[,3], X4=X[,4])
	    generate_residuals(X, G, i)
	  }
  }
