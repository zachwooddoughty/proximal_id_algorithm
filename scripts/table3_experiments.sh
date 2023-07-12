sample_sizes=4000
n_bootstrap=64
n_dgps=1
dgp_start=1
n_datasets=64
dataset_start=1
outdir="sim_results"

for zmw_effect in 0 0.2 0.4 0.8; do

    # Gaussian Z, M
    time Rscript front_door.R \
      --zmw_effect=$zmw_effect --sample_sizes=$sample_sizes \
      --save=TRUE --outdir=$outdir --n_bootstrap=$n_bootstrap \
      --n_datasets=$n_datasets --n_dgps=$n_dgps \
      --dgp_start=$dgp_start --dataset_start=$dataset_start

    # Binary Z, M
    time Rscript front_door.R \
      --zmw_effect=$zmw_effect --sample_sizes=$sample_sizes \
      --save=TRUE --outdir=$outdir --n_bootstrap=$n_bootstrap \
      --n_datasets=$n_datasets --n_dgps=$n_dgps \
      --dgp_start=$dgp_start --dataset_start=$dataset_start \
      --z_bin=TRUE --m_bin=TRUE
done
