set -e

outdir="mtx_results/"
mkdir -p $outdir
n_boot=256

dtype="binary"
for model in "naive_frontdoor" "simple_proximal" "proximal_frontdoor"; do 
    echo "$dtype $model $n_boot"
    time Rscript mtx_analysis.R --binary_treat=TRUE --delta_outcome=TRUE --a_ints=FALSE \
      --n_bootstrap=$n_boot --models="$model" &> $outdir/${model}_${dtype}_${n_boot}.txt
done

dtype="cat"
for model in "naive_frontdoor" "simple_proximal" "proximal_frontdoor"; do 
    echo "$dtype $model $n_boot"
    time Rscript mtx_analysis.R --binary_treat=FALSE --delta_outcome=TRUE --a_ints=FALSE \
      --n_bootstrap=$n_boot --models="$model" &> $outdir/${model}_${dtype}_${n_boot}.txt
done
