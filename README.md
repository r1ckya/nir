# NIR

## Quickstart

To run pipeline for cdr1 seqs run:

```[bash]
bash run_cdr1.sh
```

To generate simulation data and run dependency search run:

```[bash]
cd simulatin && bash run_simulation.sh
```

## `run.sh`

Runs pipeline on your data

Flags:

* `-f` file with stopwords
* `-i` input file with fasta seqs
* `-o` output csv file
* `-d` dir to save images
* `-a` file to safe alignment to

see example at `run_cdr1.sh`

## `simulate.py`

In simulation dependent positions are generated in pairs (1 depends on 0, 2 on 3...), also we can add positions with random letters at the end

Each 2 positions (0 1, 2 3, ...) have single unique pair of dependent letters

Usage:

```[bash]
# generates strings of 8 letters, 2 last are random
python simulate.py \
    --p_i_range 0.04 0.07 \
    --p_ji_range 0.4 0.6 \
    --n_seq 5000 \
    --n_pairs 3 \
    --n_randoms 2 \
    --seed 123 \
    > data.txt
```

## `filter.py`

Leaves unique sequences, removes stopwords and filters by sequence length

Usage:

```[bash]
python filter.py \
    --filter_seq_file data/cdr1_filter_seq.txt \
    --min_len 8 \
    --max_len 8 \
    < cdr1_aa_len500.txt \
    > cdr1_aa_len500_filtered.txt
```

## `process_aln.py`

Runs search of dependent positions

Usage:

```[bash]
python process_aln.py \
    --out_dir results \
    --ratio_cutoff 0 \
    --c_ij_cutoff 20 \
    --p_j_cond_i_cutoff auto \
    --p_ij_cutoff auto \
    --uncommon \
    < data/cdr1_clustalo.txt \
    > results.csv
```

* `--uncommon` flag removes item from results if dependent letter is most common at it's position
* `--*_cutoff` params are cutoffs for different values, some cutoffs can be set to `auto`, `auto` means that cutoff will be calculated dynamically using normal distribution assumption and 3-sigma rule

## `delta_simulation.py`

Runs simulation with different `p_ji` values

Usage:

* `--n_trials` number of runs per `p_ji` value
* `--p_ji_range` min/max values of `p_ji`
* `--num_steps` number `p_ji` points from `p_ji_range`
* `--n_proc` number of processes

```[bash]
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
```

## `extract_seqs.sh`

Extracts some seqs with highly dependent positions

Usage:

```[bash]
bash extract_seqs.sh {pipeline_results.csv} {alignment_file}
```

Example:

```[bash]
bash extract_seqs.sh data/results.csv data/cdr1_clustalo.txt > data/extracted_seqs.txt
```

Output format:

```[txt]
_ before codition position ^ before dependent position

Position (5 L) depends on position (8 L) with p=0.74
--GYT^LTE_LS-
--GHT^LSE_LS-
--GYT^LTE_LP-
```
