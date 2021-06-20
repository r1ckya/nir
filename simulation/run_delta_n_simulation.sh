#!/bin/bash

echo "Start delta n simulation run!"

mkdir -p delta_results

# for n in $(seq 1000 500 10000)
# do
#     python delta_simulation.py \
#         --p_i_range 0.04 0.2 \
#         --n_trials 3 \
#         --n_seq $n \
#         --n_pairs 3 \
#         --n_randoms 2 \
#         --seed 7 \
#         --p_ji_range 0.05 1 \
#         --num_steps 39 \
#         --n_proc 4 \
#         > delta_results/delta_results_$n.csv
# done

python plot_delta_n_results.py delta_results/delta_q_point_plot.png

echo "Run done!"
