import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6

import power
import nystroemmhsic.helpers as helpers

if __name__=="__main__":
    low, high = 200,2000
    steps = 5
    num_runs = 5
    data_func = helpers.data_MSD
    num_threads = 7
    alpha = 0.01

    power.common(low, high, steps, num_runs, data_func, num_threads, "../results/fig4_msd.csv", alpha=alpha, include_nfsic=True)
    
