# Copyright (c) 2016
# Jonas Peters  [peters@stat.math.ethz.ch]
# Niklas Pfister [pfister@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

## Power analysis (including Romano test)


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



# define function that generates a DAG G with 10 nodes from the cases:
# 1) edge b/w 1 and 2
# 2) random full DAG
# 3) no edge
# 4) one random edge
# and samples sample_size data observations from this DAG
DAG_data <- function(case,sample_size,d){
  if(case==1){
    G <- matrix(0,d,d)
    G[1,2] <- 1
  }
  else if(case==2){
    G <- as.matrix(randomDAG(d,1))
  }
  else if(case==3){
    G <- matrix(0,d,d)
  }
  else if(case==4){
    G <- matrix(0,d,d)
    ind <- 1:d
    i <- sample(ind,1)
    j <- sample(ind[-i],1)
    G[i,j] <- 1
  }
  X <- sampleDataFromG(sample_size,G,parsNoise=list(noiseExp=1,varMin=50,varMax=100,noiseExpVarMin=2,noiseExpVarMax=4,bound=rep(1,dim(G)[2])))
  return(list(G=G,X=X))
}

# define function that computes the rejection rate when testing for dhsic and pairwise 
compute_rejectionrate <- function(method, reps, sample_size){
  # initalize vectors
  rej_dhsic <- vector("numeric",2)
  rej_pairwise <- vector("numeric",2)
  rej_romano <- vector("numeric",2)
  pval_dhsic <- vector("list",2)
  pval_pairwise <- vector("list",2)
  pval_romano <- vector("list",2)
  time_dhsic <- 0
  time_pairwise <- 0
  time_romano <- 0
  # computation
  for(case in 1:2){
    pval1 <- vector("numeric",reps)
    pval2 <- vector("numeric",reps)
    pval3 <- vector("numeric",reps)
    for(i in 1:reps){
      GX <- DAG_data(case,sample_size,4)
      Xlist <-  split(GX$X, rep(1:ncol(GX$X), each = nrow(GX$X)))
      ptm <- proc.time()
      pval1[i] <- dhsic.test(Xlist,method=method,pairwise=FALSE)$p.value
      time_dhsic <- time_dhsic+(proc.time()-ptm)[3]
      ptm <- proc.time()
      pval2[i] <- dhsic.test(Xlist,method=method,pairwise=TRUE)$p.value
      time_pairwise <- time_pairwise+(proc.time()-ptm)[3]
      if(method!="gamma"){
        ptm <- proc.time()
        pval3[i] <- romano.test(Xlist,method=method,C=1000)$p.value
        time_romano <- time_romano+(proc.time()-ptm)[3]
      }
    }
    rej_dhsic[case] <- sum(pval1<0.05)/reps
    rej_pairwise[case] <- sum(pval2<0.05)/reps
    rej_romano[case] <- sum(pval3<0.05)/reps
    pval_dhsic[[case]] <- pval1
    pval_pairwise[[case]] <- pval2
    pval_romano[[case]] <- pval3
  }
  result <- list(rejection_rate=list(dhsic=rej_dhsic,pairwise=rej_pairwise,romano=rej_romano),
                 pval=list(dhsic=pval_dhsic,pairwise=pval_pairwise,romano=pval_romano),
                 time=list(dhsic=time_dhsic,pairwise=time_pairwise,romano=time_romano))
  return(result)
}


# Simulation
reps <- 1
sample_size <- c(50,100,150,200)
set.seed(1)

ptm <- proc.time()
gamma_result <- mclapply(sample_size,compute_rejectionrate,method="gamma",reps=reps,mc.cores=10)
print(proc.time()-ptm)
save(gamma_result,file='gamma_result.Rda')

ptm <- proc.time()
permutation_result <- mclapply(sample_size,compute_rejectionrate,method="permutation",reps=reps,mc.cores=10)
print(proc.time()-ptm)
save(permutation_result,file='permutation_result.Rda')

ptm <- proc.time()
bootstrap_result <- mclapply(sample_size,compute_rejectionrate,method="bootstrap",reps=reps,mc.cores=10)
print(proc.time()-ptm)
save(bootstrap_result,file='bootstrap_result.Rda')
