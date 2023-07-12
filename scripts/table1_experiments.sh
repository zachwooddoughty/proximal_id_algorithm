sample_sizes=4000,16000,64000
n_bootstrap=1
n_dgps=4
dgp_start=1
n_datasets=64
dataset_start=1
outdir="sim_results"

time Rscript proximal_front_door.R \
  --sample_sizes=$sample_sizes \
  --save=TRUE --n_bootstrap=$n_bootstrap \
  --n_datasets=$n_datasets --n_dgps=$n_dgps \
  --dgp_start=$dgp_start --dataset_start=$dataset_start \
  --outdir=$outdir
