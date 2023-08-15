# Proximal ID Algorithm Code

Code for the paper:

Shpitser, Ilya, Zach Wood-Doughty, and Eric J. Tchetgen Tchetgen. "The proximal
ID algorithm." *Journal of Machine Learning Research* 24.188 (2023): 1-46.

If you use this code in published work, please cite the paper as:

```
@article{JMLR:v24:21-0950,
  author  = {Ilya Shpitser and Zach Wood-Doughty and Eric J. Tchetgen Tchetgen},
  title   = {The Proximal ID Algorithm},
  journal = {Journal of Machine Learning Research},
  year    = {2023},
  volume  = {24},
  number  = {188},
  pages   = {1--46},
  url     = {http://jmlr.org/papers/v24/21-0950.html}
}
```

## Notes

For questions about the code, please open an issue here or email zach [at] northwestern [dot] edu.

To request access to the MTX data, please email Zach.

Some code (GMM and MRF functions for simple proximal) was adapted from the
appendix released by:

Miao, Wang, Xu Shi, and Eric Tchetgen Tchetgen. "A confounding bridge approach
for double negative control inference on causal effects." (2018).
https://arxiv.org/pdf/1808.04945.pdf

## How to use this repository

### Installation

1. Create a (conda) virtual environment (e.g., `conda create -n pid python
   r-base r-essentials`).
2. To preprocess MTX data and run or compile experiments, install Python
   requirements from `requirements.txt` (e.g., `conda install -y --file 
   requirements.txt -c conda-forge`).
3. To be able to run the estimation code, install R requirements with `Rscript
   requirements.R`.

### Rebuild the sim study tables from the paper

1. See experimental results in `sim_results/` directory.
2. Run `bash scripts/build_sim_tables.sh` to produce LaTeX for Tables 1-4.
   This will output additional tables that are not included in the paper.

### Rebuild the MTX plots from the paper

1. See experimental results in `mtx_results/` directory.
2. To produce Table 5, run `python scripts/postprocess_mtx_results.py`.
3. To produce Figure 5, compile `mtx_results/figure5.tex`.

### Rerun the sim studies from the paper

1. Run `bash scripts/tableX_experiments.sh` for tables 1 through 4.
2. This will take a while and can be parallelized.
3. When done, use `bash scripts/build_sim_tables.sh` to rebuild the tables.

### Rerun the MTX analysis from the paper

1. Email Zach to request access to the data and copy it into `mtx_data/`.
2. Run `python scripts/preprocess_mtx_data.py --delta_outcome` and
   `python scripts/preprocess_mtx_data.py --delta_outcome --binary_treat` to
   produce the two preprocessed datasets for Figure 5 and Table 5 respectively.
3. Use `bash scripts/run_mtx_analysis.sh` to run the analysis on the
   preprocessed data.
