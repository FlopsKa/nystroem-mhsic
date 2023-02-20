import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6


import numpy as np
import pandas as pd
import multiprocess as mp

import nystroemmhsic.estimators as estimators, nystroemmhsic.data as data, nystroemmhsic.kernels as kernels, nystroemmhsic.helpers as helpers
from nystroemmhsic.helpers import Testbench as tb



def get_alg(name, gamma_x, gamma_y, n_prime,alpha,seed):
    kX = kernels.RBF(gamma_x)
    kY = kernels.RBF(gamma_y)
    algs = {"N-MHSIC" : estimators.HSICNy(kernels=[kX,kY],num_nystrom_samples=n_prime,seed=seed),
            "N-MHSIC2" : estimators.HSICNy(kernels=[kX,kY],num_nystrom_samples=2*n_prime,seed=seed),
            "N-MHSIC4" : estimators.HSICNy(kernels=[kX,kY],num_nystrom_samples=4*n_prime,seed=seed),
            "N-MHSIC6" : estimators.HSICNy(kernels=[kX,kY],num_nystrom_samples=6*n_prime,seed=seed),
            "N-MHSIC8" : estimators.HSICNy(kernels=[kX,kY],num_nystrom_samples=8*n_prime,seed=seed),
            "N-MHSIC16" : estimators.HSICNy(kernels=[kX,kY],num_nystrom_samples=16*n_prime,seed=seed),
            "N-MHSIC32" : estimators.HSICNy(kernels=[kX,kY],num_nystrom_samples=32*n_prime,seed=seed),
            "N-HSIC" :  estimators.LargeScaleHSICNy( kernelX=kX,kernelY=kY, num_nystrom_samples=n_prime,),
            "RFF-HSIC" : estimators.LargeScaleHSICRFF( kernelX=kX,kernelY=kY, num_rff_samples=n_prime,),
            "NFSIC" :  estimators.FSICAdapter( kernelX=kX,kernelY=kY, seed=seed, J=5,alpha=alpha),
            "HSIC" : estimators.HSIC([kX,kY]),
    }
    return algs[name]

def runner(config):
    print(config)
    name = config["name"]
    n = config["n"]
    data_func = config["data_func"]
    gamma_x = config["gamma_x"]
    gamma_y = config["gamma_y"]
    alpha = config["alpha"]
    seed = config["seed"]
    alg = get_alg(name, gamma_x, gamma_y,n_prime=int(np.sqrt(n)), alpha=alpha,seed=seed)
    p, r = tb(seed=seed, data_func=data_func, alpha=alpha).get_power_runtime(alg.hsic, n=n)
    result = pd.DataFrame({
        "alg" : [name],
        "n" : [n],
        "power" : [p],
        "runtime" : [r]
    })
    return result

def common(low, high, steps, num_runs, data_func, num_threads, out_file, include_nfsic=True, alpha=0.05):
    rng= np.random.default_rng()
    X,Y=data_func(1000,rng)

    gamma_x = helpers.g(helpers.est_sigma(X))
    gamma_y = helpers.g(helpers.est_sigma(Y)) 
    
    df_nfsic = pd.DataFrame()
    if include_nfsic:
        for _ in range(num_runs):
            for n in np.linspace(low,high,num=steps,dtype=int):
                print(n)
                alg = get_alg("NFSIC", gamma_x, gamma_y,n_prime=int(np.sqrt(n)),alpha=alpha,seed=rng.integers(100000))
                
                power = []
                with helpers.ContextTimer() as t:
                    n_samples = 100
                    for s in range(n_samples):
                        print(f"{s+1}/{n_samples}")
                        X,Y = data_func(n, rng)
                        power += [alg.test(X,Y)]
                    p = np.mean(power)
                r = t.secs / n_samples
                df_nfsic = pd.concat([df_nfsic, pd.DataFrame({
                    "alg" : ["NFSIC"],
                    "n" : [n],
                    "power" : [p],
                    "runtime" : [r]
                }
                )])
    print(df_nfsic)
    
    args = []
    rng = np.random.default_rng(1234)
    for _ in range(num_runs):
        for n in np.linspace(low,high,num=steps,dtype=int):        
            for a in ["N-MHSIC4", "N-MHSIC6", "N-MHSIC8", "N-HSIC", "RFF-HSIC"]:# "HSIC"]:
                t1 = {"name" : a,
                      "gamma_x" : gamma_x,
                      "gamma_y" : gamma_y,
                      "n" : n,
                      "data_func" : data_func,
                      "alpha" : alpha,
                      "seed" : rng.integers(100000)
                }
                args += [(t1,)]

    df=pd.DataFrame()
   
    with helpers.ContextTimer() as t:
        with mp.Pool(num_threads) as p:
            df = pd.concat(p.starmap(runner, args)).reset_index(drop="True")
    res = pd.concat([df_nfsic, df])
    print(res)
   
    print("Runtime: %.2f seconds." % t.secs)
    res.to_csv(out_file)
