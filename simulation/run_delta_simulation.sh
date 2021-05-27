#!/bin/bash

echo "Start delta simulation run!"

mkdir -p delta_results

python delta_simulation.py \
    --p_i_range 0.04 0.1 \
    --n_trials 3 \
    --n_seq 1000 \
    --n_pairs 3 \
    --n_randoms 2 \
    --seed 7 \
    --p_ji_range 0.05 1 \
    --num_steps 77 \
    --n_proc 4 \
    > delta_results/delta_results.csv

python plot_delta_results.py \
    delta_results/delta_results.csv \
    delta_results/delta_acc_plot.png

echo "Run done!"
