import numpy as np
import pandas as pd
import time

import nystroemmhsic.estimators as estimators, nystroemmhsic.data as data, nystroemmhsic.kernels as kernels
from nystroemmhsic.helpers import g, est_sigma

if __name__=="__main__":
    # config
    low, high = 100, 1000
    steps = 10
    runs = 100

    
    X, Y = data.indep(n=500, d=2)
    rbf = kernels.RBF(g(est_sigma(X)))
    results = pd.DataFrame()
    for i in range(runs):
        for n in np.linspace(low,high,num=steps,dtype=int):
            print(n)
            algs = [
                ("N-MHSIC", estimators.HSICNy(kernels=[rbf, rbf], num_nystrom_samples=2*int(np.sqrt(n)))),
                ("N-HSIC", estimators.LargeScaleHSICNy(kernelX=rbf, kernelY=rbf, num_nystrom_samples=2*int(np.sqrt(n)))),
                ("RFF-HSIC", estimators.LargeScaleHSICRFF(kernelX=rbf, kernelY=rbf, num_rff_samples=2*int(np.sqrt(n)))),
            ]
            X, Y = data.indep(n=n, d=2, seed=n+i)
            for (name, impl) in algs:
                start = time.time()
                estimate = impl.hsic(X,Y)
                runtime = time.time() - start
                results = pd.concat([results, pd.DataFrame({"n" : [n], "Algorithm" : [name], "Estimate" : [estimate], "Runtime" : [runtime]})])

        print(results)

    results.to_csv("../results/fig1.csv")
