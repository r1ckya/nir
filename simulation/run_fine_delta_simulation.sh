#!/bin/bash

echo "Start delta simulation run!"

mkdir -p delta_results

python delta_simulation.py \
    --p_i_range 0.04 0.2 \
    --n_trials 10 \
    --n_seq 1000 \
    --n_pairs 3 \
    --n_randoms 2 \
    --seed 7 \
    --p_ji_range 0.3 0.5 \
    --num_steps 21 \
    --n_proc 4 \
    > delta_results/fine_delta_results.csv

python plot_delta_results.py \
    delta_results/fine_delta_results.csv \
    delta_results/fine_delta_acc_plot.png

echo "Run done!"
