#!/bin/bash

echo "Running pipeline for cdr1 regions..."
./run.sh \
    -f data/cdr1_filter_seq.txt \
    -i data/cdr1_aa_len500.txt \
    -o data/results.csv \
    -d results \
    -a data/cdr1_clustalo.txt
echo "cdr1 pipeline done!"
