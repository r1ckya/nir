import argparse
import csv
import sys
from collections import defaultdict
from os.path import join
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def get_counts_single(
    records: List[Seq],
) -> Tuple[Dict[Tuple[int, str], int], np.ndarray]:
    n = len(records[0])
    counts = defaultdict(int)
    counts_pos = np.zeros(n, dtype=int)
    for record in records:
        for i, aa in enumerate(record.seq):
            if aa != "-":
                counts[(i, aa)] += 1
                counts_pos[i] += 1
    return counts, counts_pos


def get_counts_double(
    records: List[Seq],
) -> Tuple[Dict[Tuple[int, str], int], np.ndarray]:
    counts = defaultdict(int)
    n = len(records[0])
    counts_pos = np.zeros((n, n), dtype=int)
    for record in records:
        for i, aa_i in enumerate(record.seq):
            for j, aa_j in enumerate(record.seq):
                if i == j or aa_i == "-" or aa_j == "-":
                    continue
                counts[(i, j, aa_i, aa_j)] += 1
                counts_pos[i, j] += 1
    return counts, counts_pos


def is_most_common(
    i: int, aa_i: str, cnt_single: Dict[Tuple[int, str], int]
) -> bool:
    c_i = cnt_single[(i, aa_i)]
    for (j, aa_j), c_j in cnt_single.items():
        if j == i and aa_j != aa_i and c_i < c_j:
            return False
    return True


def n_aa_at_pos(i: int, cnt_single: Dict[Tuple[int, str], int]) -> int:
    res = set()
    for j, aa_j in cnt_single.keys():
        if j == i and aa_j != "-":
            res.add(aa_j)
    return len(res)


def show(info: List[Tuple], out_dir: str) -> None:
    ratios = list(map(lambda x: x[1], info))
    conds = list(map(lambda x: x[0], info))
    c_ij = list(map(lambda x: x[-1], info))

    plt.scatter(ratios, conds)
    plt.xlabel("ratio")
    plt.ylabel("cond")
    plt.savefig(join(out_dir, "scatter.png"))
    plt.close()

    plt.hist(ratios, bins=100)
    plt.title("ratios hist")
    plt.savefig(join(out_dir, "ratios_hist.png"))
    plt.close()

    plt.hist(c_ij, bins=100)
    plt.title("c_ij hist")
    plt.savefig(join(out_dir, "c_ij_hist.png"))
    plt.close()

    plt.hist(conds, bins=100)
    plt.xticks(np.linspace(0, 1, 11))
    plt.tight_layout()
    plt.title("conds hist")
    plt.savefig(join(out_dir, "conds_hist.png"))
    plt.close()

    plt.plot(sorted(ratios))
    plt.title("ratios_log")
    plt.yscale("log")
    plt.savefig(join(out_dir, "ratios_log.png"))
    plt.close()

    plt.plot(sorted(ratios))
    plt.title("ratios")
    plt.savefig(join(out_dir, "ratios.png"))
    plt.close()

    plt.plot(sorted(conds))
    plt.title("conds")
    plt.hlines(1 / 20, 0, len(conds), "r", label=f"y=1/20")
    plt.legend()
    plt.savefig("results/conds.png")
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", type=str, required=True)
    parser.add_argument(
        "--p_ij_cutoff", type=str, default="auto", help="auto or float(0.006)"
    )
    parser.add_argument(
        "--p_j_cond_i_cutoff",
        type=str,
        default="auto",
        help="auto or float(0.15)",
    )
    parser.add_argument("--ratio_cutoff", type=float, default=0.0)
    parser.add_argument("--c_ij_cutoff", type=int, default=0)

    parser.add_argument(
        "--uncommon",
        action="store_true",
        help="discard dependent letters from result if it is most common at position",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    records = list(SeqIO.parse(sys.stdin, "fasta"))

    cnt_single, cnt_pos_single = get_counts_single(records)
    cnt_double, cnt_pos_double = get_counts_double(records)

    info = []

    for (i, j, aa_i, aa_j), c_ij in cnt_double.items():

        c_i = cnt_single[(i, aa_i)]
        c_j = cnt_single[(j, aa_j)]

        p_ij = c_ij / cnt_pos_double[i, j]
        p_i = c_i / cnt_pos_single[i]
        p_j = c_j / cnt_pos_single[j]

        p_j_cond_i = c_ij / c_i

        ratio = p_ij / (p_i * p_j)

        if args.p_ij_cutoff == "auto":
            n_aa_i = n_aa_at_pos(i, cnt_single)
            n_aa_j = n_aa_at_pos(j, cnt_single)
            mean = (1 / n_aa_i) * (1 / n_aa_j)
            std = (mean * (1 - mean) / cnt_pos_double[i, j]) ** 0.5
            p_ij_cutoff = mean + 3 * std
        else:
            p_ij_cutoff = float(args.p_ij_cutoff)

        if args.p_j_cond_i_cutoff == "auto":
            n_aa_j = n_aa_at_pos(j, cnt_single)
            mean = 1 / n_aa_j
            std = (mean * (1 - mean) / c_i) ** 0.5
            p_j_cond_i_cutoff = mean + 3 * std
        else:
            p_j_cond_i_cutoff = float(args.p_ij_cutoff)

        if (
            not (args.uncommon and is_most_common(j, aa_j, cnt_single))
            and p_ij > p_ij_cutoff
            and p_j_cond_i > p_j_cond_i_cutoff
            and ratio > args.ratio_cutoff
            and c_ij > args.c_ij_cutoff
        ):
            info.append(
                (
                    p_j_cond_i,
                    ratio,
                    i,
                    aa_i,
                    j,
                    aa_j,
                    p_i,
                    p_j,
                    p_ij,
                    c_i,
                    c_j,
                    c_ij,
                )
            )

    if args.verbose:
        print("number of different amino acids per position")
        for i in range(len(records[0])):
            print(f"{i:>2} {n_aa_at_pos(i, cnt_single),:>2}", file=sys.stderr)

    writer = csv.writer(sys.stdout)
    writer.writerow(
        [
            "p_j_cond_i",
            "ratio",
            "i",
            "aa_i",
            "j",
            "aa_j",
            "p_i",
            "p_j",
            "p_ij",
            "c_i",
            "c_j",
            "c_ij",
        ]
    )
    writer.writerows(sorted(info))

    show(info, args.out_dir)


if __name__ == "__main__":
    main()
