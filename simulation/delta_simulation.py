import argparse
import csv
import itertools
import os
import sys
import tempfile
from multiprocessing import Pool

import numpy as np
from tqdm import tqdm

from simulate import AminoAcids


def chunks(iterable, size):
    it = iter(iterable)
    while True:
        group = tuple(itertools.islice(it, size))
        if not group:
            break
        yield group


def is_dependent(i: int, j: int, aa_i: str, aa_j: str, n_dep_pos: int) -> bool:
    return (
        i % 2 == 0
        and j % 2 == 1
        and j - i == 1
        and i < n_dep_pos
        and j < n_dep_pos
        and i == AminoAcids[aa_i].value - 1
        and j == AminoAcids[aa_j].value - 1
    )


def get_acc(res_file, n_dep_pos):
    with open(res_file, "r") as f:
        reader = csv.DictReader(f)
        n_rows = 0
        n_drop = 0
        drop = True

        for row in reader:
            i, j, aa_i, aa_j = (
                int(row["i"]),
                int(row["j"]),
                row["aa_i"],
                row["aa_j"],
            )
            if drop:
                if is_dependent(i, j, aa_i, aa_j, n_dep_pos):
                    drop = False
                else:
                    n_drop += 1
            n_rows += 1
        return 0 if n_rows == n_drop else n_dep_pos / (n_rows - n_drop)


def run(
    p_ji: float,
    p_i_min: float,
    p_i_max: float,
    n_seq: int,
    n_pairs: int,
    n_randoms: int,
    seed: int,
) -> float:

    data_fd, data_file = tempfile.mkstemp()
    res_fd, res_file = tempfile.mkstemp()

    os.system(
        f"""python simulate.py \
            --p_i_range {p_i_min} {p_i_max} \
            --p_ji_range {p_ji} {p_ji} \
            --n_seq {n_seq} \
            --n_pairs {n_pairs} \
            --n_randoms {n_randoms} \
            --seed {seed} \
            > {data_file} \
            2> /dev/null
            """
    )
    os.system(
        f"""python ../process_aln.py \
            --out_dir /tmp \
            --c_ij_cutoff 0 \
            --p_ij_cutoff auto \
            --p_j_cond_i_cutoff auto \
            < {data_file} \
            > {res_file} \
            2> /dev/null
            """
    )

    acc = get_acc(res_file, 2 * n_pairs)

    os.close(data_fd)
    os.close(res_fd)
    os.unlink(res_file)
    os.unlink(data_file)

    return acc


def run_unpack(kwargs):
    return run(**kwargs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p_i_range", nargs=2, type=float, required=True)
    parser.add_argument("--n_trials", type=int, required=True)
    parser.add_argument("--n_seq", type=int, required=True)
    parser.add_argument("--n_pairs", type=int, required=True)
    parser.add_argument("--n_randoms", type=int, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--p_ji_range", nargs=2, type=float, required=True)
    parser.add_argument("--num_steps", type=int, required=True)
    parser.add_argument("--n_proc", type=int, default=1)

    args = parser.parse_args()

    p_i_min, p_i_max = args.p_i_range

    np.random.seed(args.seed)

    p_ji_space = np.linspace(*args.p_ji_range, num=args.num_steps)

    common_params = dict(
        p_i_min=p_i_min,
        p_i_max=p_i_max,
        n_seq=args.n_seq,
        n_pairs=args.n_pairs,
        n_randoms=args.n_randoms,
    )

    with Pool(args.n_proc) as pool:
        runs_params = (
            dict(
                p_ji=p_ji,
                seed=args.seed + k * args.n_trials + trial,
                **common_params,
            )
            for k, p_ji in enumerate(p_ji_space)
            for trial in range(args.n_trials)
        )
        results = pool.imap(run_unpack, runs_params)
        results = tqdm(results, total=args.num_steps * args.n_trials)

        writer = csv.DictWriter(sys.stdout, fieldnames=["p_ji", "acc"])
        writer.writeheader()
        for p_ji, trials_acc in zip(
            p_ji_space, chunks(results, args.n_trials)
        ):
            acc = sum(trials_acc) / args.n_trials
            writer.writerow(dict(p_ji=p_ji, acc=acc))


if __name__ == "__main__":
    main()
