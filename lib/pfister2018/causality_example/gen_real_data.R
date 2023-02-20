# Copyright (c) 2016
# Jonas Peters  [jonas.peters@math.ku.dk]
# Niklas Pfister [pfister@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.

## Real World Example


# load required libraries
library(dHSIC)
library(MASS)
library(mgcv)

# load romano.test
source("../extra/romano.test.R")


# read data
example_data <- read.table("data.csv",sep=" ",header = TRUE)
example_data <- as.matrix(example_data)
example_data <- example_data[sample(1:349,349,replace=FALSE),]
X <- list(X1=example_data[,1],X2=example_data[,2],X3=example_data[,3])
# load DAG file
load("allDagsWith3Nodes.RData")

d <- 3
numDAGs <- nrow(allDags)

# initalize pvalue vector
pval_dhsic <- vector("numeric",numDAGs)
pval_pairwise <- vector("numeric",numDAGs)
pval_romano <- vector("numeric",numDAGs)

# compute the p-value for each DAG
for(i in 1:numDAGs){
  print(i)
  # generate Adj-matrix for DAG
  A <- matrix(allDags[i,],d,d)
  # initialize residual vector
  res <- vector("list",d)
  # regress each node on its variables and collect the residuals
  for(j in 1:d){
    if(sum(A[,j])==0){
      res[[j]] <- X[[j]]
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
    }
  }
  p <- paste("../../../../data/causal_real/" , i,".csv",sep="")
  print(p)
  write.csv(res, p)
}

pval <- list(dhsic=pval_dhsic,pairwise=pval_pairwise,romano=pval_romano)

save(pval,file="pval_realworld.Rda")



