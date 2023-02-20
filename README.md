Instructions for Reproducing the Experiments
============================================

# Data

The data must be obtained from the sources within the article. For the causality experiments, see below.

## Causality Data

We generate the data to run the causality experiments with R, using code from [Pfister, 2018], save it as `.csv` and continue with python. The code is in `lib/pfister2018`.

### Simulated Causality Data

Note that we change line 214 of the original file `sampleDataFromG.R` to obtain samples with larger variance:

    noisetmp <- (1*sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)


Create the folder structure with

    cd data
    mkdir -p causal_sim/500/{1..100} causal_sim/1000/{1..100} causal_sim/1500/{1..100}


And generate the data (this may take some time) by

    cd causality_example
    Rscript gen_sim_data.R


### Weather Data

The commands to obtain the residuals are

    cd data
    mkdir causal_real
    cd causality_example
    Rscript gen_real_data.R

# Experiments

The code to reproduce the experiments is in `figxx.py`, respectively. The plots can be created with `plot_figxx.ipynb`.

For NFSIC, we use the code from [fsic-test](https://github.com/wittawatj/fsic-test) with commit id `1bc6318` [Jitkrittum, 2017].


[Pfister, 2018] Niklas Pfister, Peter Bühlmann, Bernhard Schölkopf, and Jonas Peters. Kernel-based tests for joint independence. Journal of the Royal Statistical Society. Series B. Statistical Methodology, pages 5–31, 2018.

[Jitkrittum, 2017] Jitkrittum, Wittawat, Zoltán Szabó, and Arthur Gretton. An adaptive test of independence with analytic kernel embeddings. International Conference on Machine Learning. PMLR, 2017.