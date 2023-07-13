Each file in this directory represents the output of one simulation study. The
title tells you which experiment was run. For example, in the file:

```
z_bi0-m_bi0-ay_e0-uz_eNA-uw_eNA-zmw_NA-n_bo64-dgp1-data1-samples4000.csv
```

This can be interpreted as:
- `z_bi0`: Is Z binary? 0 means no (Gaussian), 1 means yes. This is relevant to Table 3.
- `m_bi0`: Is M binary? Same as above.
- `ay_e0`: the direct effect of A on Y is set to be 0. This effect is unconstrained ("NA") for all experiments except those in Table 2.
- `uz_eNA`: The direct effect of U on Z is unconstrained (NA). This effect is unconstrained for all experiments except those in Table 4.
- `uw_eNA`: The direct effect of U on W is unconstrained (NA). Same as above.
- `zmw_NA`: The direct effects of Z on M and of M on W are unconstrained. These effects are unconstrained for all experiments except those in Table 3.
- `n_bo64`: The experiment uses 64 bootstrap resamplings. Table 1 does not use bootstrap resampling; all other tables use 64 resamplings. 
- `dgp1`: The DGP's parameters were randomly chosen with seed `1`. We consider four DGPs with different parameters chosen this way.
- `data1`: The dataset was randomly sampled (independently of the DGPs parameters) with seed `1`. We sample 64 different datasets for each of our four DGPs.
- `samples4000`: The dataset contains 4000 samples. Table 1 varies the sample size; all other tables use 4000 samples.

Within each file, the first two lines are `//` comments in JSON format. The
first line contains the arguments of the experiments; the second line shows the
parameters of the DGP that was used to sample the dataset. The subsequent lines
of the file contain valid CSV format showing the experimental results of each
method's estimate of the average causal effect of A on Y. For experiments with
bootstrap resampling, "width" indicates the width of the 95% bootstrap
confidence interval, and "coverage" is a binary indicator of whether the true
causal effect was contained within that interval.
