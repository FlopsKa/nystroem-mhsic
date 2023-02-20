# Copyright (c) 2016
# Jonas Peters  [peters@stat.math.ethz.ch]
# Niklas Pfister [pfister@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 


benchmarking <- function(n,sample_size,simulation_fun,Bset,kernels,pairwise,cores){
  # define function that computes rejection rate
  rejectionrate_fun <- function(sample_size,method){
    crit <- vector("numeric",n)
    stat <- vector("numeric",n)
    for(i in 1:n){
      X <- simulation_fun(sample_size)
      temp <- dhsic.test(X,method=method,kernel=kernels,B=Bset,pairwise=pairwise)
      crit[i] <- temp$crit.val
      stat[i] <- temp$statistic
    }
    return(sum(stat>crit)/n)
  }
  # compute rejection rates
  print("computing permutation level")
  rej_rate_perm <- mclapply(sample_size,rejectionrate_fun,"permutation",mc.cores=cores)
  print("computing bootstrap level")
  rej_rate_boot <- mclapply(sample_size,rejectionrate_fun,"bootstrap",mc.cores=cores)
  print("computing gamma level")
  rej_rate_gamma <- mclapply(sample_size,rejectionrate_fun,"gamma",mc.cores=cores)
  # collect results
  rejectionrate <- list(permutation=rej_rate_perm,
                        bootstrap=rej_rate_boot,
                        gamma=rej_rate_gamma)
  return(rejectionrate)  
}
