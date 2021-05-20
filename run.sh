#!/bin/bash

python filter.py < data/cdr1_aa_len500.txt > data/cdr1_aa_unique.txt
# clustlao cdr1_aa_unique.txt to cdr1_aa_len500aln.txt in fasta format
python process_aln.py < data/cdr1_aa_len500aln.txt > results/ratios_uncommon.txt
