#!/bin/bash

echo "Start simulation run!"

mkdir -p results

echo "Generation data..."
# default p_i is 1/20 = 0.05
python simulate.py \
    --p_i_range 0.04 0.07 \
    --p_ji_range 0.4 0.6 \
    --n_seq 1000 \
    --n_pairs 3 \
    --n_randoms 2 \
    > data.txt \
    2> /dev/null
echo "Data generated!"

echo "Processing alignment..."
python ../process_aln.py \
    --out_dir results \
    --c_ij_cutoff 10 \
    --p_ij_cutoff auto \
    --p_j_cond_i_cutoff auto \
    < data.txt \
    > results.csv \
    2> /dev/null
echo "Processing done!"

echo "Run done!"
