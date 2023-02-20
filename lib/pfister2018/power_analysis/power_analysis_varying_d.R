# Copyright (c) 2016
# Jonas Peters  [peters@stat.math.ethz.ch]
# Niklas Pfister [pfister@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

## Power analysis


# load required libraries
library(dHSIC)
require(Matrix)
require(MASS)
library(parallel)


# required R scripts
source("../extra/romano.test.R")
source("../extra/randomDAG.R")
source("../extra/sampleDataFromG.R")
source("../extra/randomB.R")
source("../extra/computeCausOrder.R")
source("../extra/computeGaussKernel.R")



# define function that generates a DAG G with 10 nodes:
DAG_data <- function(sample_size,d){
  G <- as.matrix(randomDAG(d,1))
  X <- sampleDataFromG(sample_size,G,parsNoise=list(noiseExp=1,varMin=50,varMax=100,noiseExpVarMin=2,noiseExpVarMax=4,bound=rep(1,dim(G)[2])))
  return(list(G=G,X=X))
}

# define function that computes the rejection rate when testing for dhsic and pairwise 
compute_rejectionrate <- function(d, method, reps, sample_size){
  # initalize pvalue vectors
  # computation
  pval1 <- vector("numeric",reps)
  pval2 <- vector("numeric",reps)
  pval3 <- vector("numeric",reps)
  pval4 <- vector("numeric",reps)
  for(i in 1:reps){
    GX <- DAG_data(sample_size,d)
    Xlist <-  split(GX$X, rep(1:ncol(GX$X), each = nrow(GX$X)))
    pval1[i] <- dhsic.test(Xlist,method=method,pairwise=FALSE,B=100)$p.value
    pval2[i] <- dhsic.test(Xlist,method=method,pairwise=TRUE,B=100)$p.value
    pval3[i] <- romano.test(Xlist,method=method,B=100,C=100)$p.value
    pval4[i] <- romano.test(Xlist,method=method,B=100,C=1000)$p.value
  }
  rej_dhsic <- sum(pval1<0.05)/reps
  rej_pairwise <- sum(pval2<0.05)/reps
  rej_romano <- sum(pval3<0.05)/reps
  rej_romano2 <- sum(pval4<0.05)/reps
  pval_dhsic <- pval1
  pval_pairwise <- pval2
  pval_romano <- pval3
  pval_romano2 <- pval4
  result <- list(rejection_rate=list(dhsic=rej_dhsic,pairwise=rej_pairwise,romano100=rej_romano,romano1000=rej_romano2),
                 pval=list(dhsic=pval_dhsic,pairwise=pval_pairwise,romano100=pval_romano,romano1000=pval_romano2))
  return(result)
}


# Simulation
reps <- 1000
sample_size <- 100
d <- c(4,6,8,10)
set.seed(1)

ptm <- proc.time()
power_analysis_varying_d <- mclapply(d,compute_rejectionrate,sample_size=sample_size,method="bootstrap",reps=reps,mc.cores=10)
print(proc.time()-ptm)
save(power_analysis_varying_d,file="power_analysis_varying_d.Rda")

