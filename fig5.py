import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from pathlib import Path
from helpers import g, est_sigma
#import seaborn as sns
import multiprocessing as mp

import estimators, kernels

rng = np.random.default_rng(12)
data_dir = Path("../../data/causal_real/")


def p_val(est, Xs, null_runs=100, rep=False):
    stat = np.mean([est.mw_hsic(Xs) for _ in range(100)])
    n = len(Xs[0])
    null_stats = []
    for _ in range(null_runs):
        idxs = [rng.choice(np.arange(n), size=n, replace=rep) for _ in Xs]
        Xs_null_dist = [X[idx] for X, idx in zip(Xs, idxs)]
        null_stats += [est.mw_hsic(Xs_null_dist)]

    return (np.count_nonzero(null_stats >= stat)+1) / (null_runs+1)

def p_val_for_dag(dag, null_runs, data_dir):
    print(dag)
    df = pd.read_csv(data_dir / (str(dag) + ".csv"), 
                     index_col=0)
    Xs = [df[col].values.reshape(-1,1) for col in df.columns]
    ks = [kernels.RBF(g(est_sigma(X))) for X in Xs]

    dhsic_est = estimators.HSIC(ks) 
    dhsic_p = p_val(dhsic_est, Xs,null_runs=null_runs)
    ny_est = estimators.HSICNy(ks,num_nystrom_samples_func=lambda _ : 100)
    ny_est_p = p_val(ny_est, Xs,null_runs=null_runs)
    return pd.DataFrame({"alg" : ["dHSIC", "Nyström HSIC"], "dag" : [dag, dag], "pval" : [dhsic_p, ny_est_p]})
    

def primed(dag):
    return p_val_for_dag(dag=dag,null_runs=500,data_dir=data_dir)

if __name__ == "__main__":
    start = time.time()
    with mp.pool.Pool(6) as p:
        is_correct = p.map(primed, range(1,26))

    print("runtime: ", time.time()- start)

    res = pd.concat(is_correct)
    print(res)

    res.to_csv("../../results/causal_real_world.csv")
