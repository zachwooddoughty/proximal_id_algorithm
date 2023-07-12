## MTX Analysis Results

This directory contains the results of the analysis contained in Section 8. 

* To produce the `*.txt` files, use `bash scripts/run_mtx_analysis.sh`

* To produce Table 5, use `python scripts/postprocess_mtx_results.py`.

* To produce Figure 5, compile `figure5.tex`; note that this uses the `*.dat`
  files which may have been created/updated by the previous step.

The proximal frontdoor estimator is much slower than the naive frontdoor and
simple proximal estimators, which is particularly noticable when using 256
bootstrap resamples.
