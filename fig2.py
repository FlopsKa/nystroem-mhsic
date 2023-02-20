import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6

import power
import nystroemmhsic.helpers as helpers

if __name__=="__main__":
    low, high = 20, 500
    steps = 10
    num_runs = 5
    data_func = helpers.gen_data
    num_threads = 7 # not for nfsic

    power.common(low, high, steps, num_runs, data_func, num_threads, "../results/fig2.csv", include_nfsic=True)

    

