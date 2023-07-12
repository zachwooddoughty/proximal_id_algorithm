agg="meanabs"
latex="--latex"  # Uncomment to print LaTeX table
# latex=""  # Uncomment to print plaintext tables

echo "Table 1: Varied Sample Size"
python scripts/compile_sims.py sim_results/ 1 --n_bootstrap=1 --agg=$agg $latex
echo ""

echo "Table 2: Varied A->Y Effect"
python scripts/compile_sims.py sim_results/ 2 --n_bootstrap=64 --agg=$agg $latex
echo ""

echo "Table 3a: Varied Z->M->W Effect, Gaussian Z, M"
python scripts/compile_sims.py sim_results/ 3 --n_bootstrap=64 --agg=$agg $latex
echo ""

echo "Table 3b: Varied Z->M->W Effect, Binary Z, M"
python scripts/compile_sims.py sim_results/ 3 --n_bootstrap=64 --agg=$agg $latex --zm_bin
echo ""

echo "Table 4: Varied Z<-U->W Effect"
python scripts/compile_sims.py sim_results/ 4 --n_bootstrap=64 --agg=$agg $latex
