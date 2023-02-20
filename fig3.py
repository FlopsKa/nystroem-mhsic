import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6


import numpy as np
import pandas as pd
from pathlib import Path
from nystroemmhsic.helpers import g, est_sigma
import time

import multiprocessing as mp

import nystroemmhsic.estimators as estimators, nystroemmhsic.data as data, nystroemmhsic.kernels as kernels, nystroemmhsic.helpers as helpers

rng = np.random.default_rng(12)


def get_alg(name, gammas, sample_size):
    ks = [kernels.RBF(g) for g in gammas]
    
    algs = {"1" : estimators.HSICNy(kernels=ks,num_nystrom_samples=int(1*np.sqrt(sample_size))),
            "2" : estimators.HSICNy(kernels=ks,num_nystrom_samples=int(2*np.sqrt(sample_size))), 
            "4" : estimators.HSICNy(kernels=ks,num_nystrom_samples=int(4*np.sqrt(sample_size))),
            "8" : estimators.HSICNy(kernels=ks,num_nystrom_samples=int(8*np.sqrt(sample_size))), 
            "HSIC" : estimators.HSIC(ks),
    }
    return algs[name]


def r(config):
    print(config)
    alg = config["alg"]
    sample_id = config["sample_id"]
    num_samples = config["num_samples"]
    null_runs = config["null_runs"]
    data_dir = config["data_dir"]

    start = time.time()
    correct = finds_correct_dag(alg, num_samples, null_runs, sample_id, data_dir)
    runtime = time.time() - start

    return pd.DataFrame({
        "alg" : [alg],
        "sample_id" : [sample_id],
        "num_samples": [num_samples],
        "correct" : [correct],
        "runtime" : [runtime]
        })

def p_val(est, Xs, null_runs=100, rep=False):
    stat = np.mean([est.mw_hsic(Xs) for _ in range(1)])
    n = len(Xs[0])
    null_stats = []
    for _ in range(null_runs):
        idxs = [rng.choice(np.arange(n), size=n, replace=rep) for _ in Xs]
        Xs_null_dist = [X[idx] for X, idx in zip(Xs, idxs)]
        null_stats += [est.mw_hsic(Xs_null_dist)]

    return (np.count_nonzero(null_stats >= stat)+1) / (null_runs+1)

def finds_correct_dag(alg, sample_size, null_runs, draw, data_dir):
    res = pd.DataFrame()
    for dag in range(1,25):
        df = pd.read_csv(data_dir / str(sample_size) / str(draw) / (str(dag) + ".csv"), 
                index_col=0, 
                names=["X1","X2","X3", "X4"], 
                skiprows=1)
        Xs = [df[col].values.reshape(-1,1) for col in df.columns]
        gammas = [g(est_sigma(X,max_len=100)) for X in Xs] ## hier nochmal prüfen, ob anderer Schätzer besser ist.


        est = get_alg(alg, gammas, sample_size)
        p = p_val(est, Xs,null_runs=null_runs)
        res = pd.concat([res, pd.DataFrame({"dag" : [dag], "pval" : [p]})])
    correct_dag = pd.read_csv(data_dir / str(sample_size) / str(draw) / "correct.csv",index_col=0)["x"].values[0]

    res.reset_index(inplace=True,drop=True)

    return res.iloc[res["pval"].idxmax()]["dag"] == correct_dag


if __name__ == "__main__":
    null_runs = 250
    
    args = []
    for n in [1500]:
        for i in range(1,101):
            for a in ["1", "2", "4", "8", "HSIC"]:
                t1 = {"alg" : a,
                      "sample_id" : i,
                      "num_samples" : n,
                      "null_runs" : null_runs,
                      "data_dir" : Path("../data/causal_weak1")
                }
                args += [(t1,)]
    with helpers.ContextTimer() as t:
        with mp.Pool(30) as p:
            df = pd.concat(p.starmap(r, args)).reset_index(drop="True")
    print(df)
    print("Runtime: %.2f seconds." % t.secs)
    df.to_csv("../results/fig3.csv")
    
        
    
