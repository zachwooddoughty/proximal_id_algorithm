This directory is intentionally left empty. If you want to work with the data,
please email Zach. Then, if you want to reproduce the MTX experiments from the paper,
put the `mtx.csv.gz` file in this directory.

Then, run `python scripts/preprocess_mtx_data.py --delta_outcome` to produce
`mtx_cont.m6_w7_y12.ydelta.csv.gz` which is used in Figure 5 and run 
`python scripts/preprocess_mtx_data.py --delta_outcome --binary_treat` to produce
`mtx_bin.m6_w7_y12.ydelta.csv.gz` which is used in Table 5.

Then, you can run `bash scripts/run_mtx_analysis.sh` to rerun the experiments.
