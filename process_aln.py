import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO


def get_counts_single(records):
    counts_single = defaultdict(int)
    for record in records:
        for i, aa in enumerate(record.seq):
            counts_single[(i, aa)] += 1
    return counts_single


def get_counts_double(records):
    counts_double = defaultdict(int)
    for record in records:
        for i, aa_i in enumerate(record.seq):
            for j, aa_j in enumerate(record.seq):
                if i == j:
                    continue
                counts_double[(i, j, aa_i, aa_j)] += 1
    return counts_double


def is_most_common(i, aa_i, cnt_single):
    c_i = cnt_single[(i, aa_i)]
    for (j, aa_j), c_j in cnt_single.items():
        if j == i and aa_j != aa_i and c_i < c_j:
            return False
    return True


def n_aa_at_pos(i, cnt_single):
    res = set()
    for j, aa_j in cnt_single.keys():
        if j == i:
            res.add(aa_j)
    return len(res)


def show(info):
    ratios = list(map(lambda x: x[1], info))
    conds = list(map(lambda x: x[0], info))

    plt.scatter(ratios, conds)
    plt.xlabel("ratio")
    plt.ylabel("cond")
    plt.savefig("results/scatter.png")
    plt.close()

    plt.hist(ratios)
    plt.title("ratios hist")
    plt.savefig("results/ratios_hist.png")
    plt.close()

    plt.hist(conds, bins=20)
    plt.xticks(np.linspace(0, 1, 11))
    plt.tight_layout()
    plt.title("conds hist")
    plt.savefig("results/conds_hist.png")
    plt.close()

    plt.plot(sorted(ratios))
    plt.title("ratios_log")
    plt.yscale("log")
    plt.savefig("results/ratios_log.png")
    plt.close()

    plt.plot(sorted(ratios))
    plt.title("ratios")
    plt.savefig("results/ratios.png")
    plt.close()

    plt.plot(sorted(conds))
    plt.title("conds")
    plt.hlines(1 / 20, 0, len(conds), "r", label=f"y=1/20")
    plt.legend()
    plt.savefig("results/conds.png")
    plt.close()


def main():
    records = []

    for record in SeqIO.parse(sys.stdin, "fasta"):
        if (
            record.seq.startswith("--")
            and record.seq.startswith("-")
            and len(record.seq.strip("-")) == 8
        ):
            record.seq = record.seq[2:-1]
            records.append(record)

    cnt_single = get_counts_single(records)
    cnt_double = get_counts_double(records)

    n_records = len(records)

    info = []

    for (i, j, aa_i, aa_j), c_ij in cnt_double.items():

        # if i > j:
        #     continue
        # if is_most_common(i, aa_i, cnt_single):
        #     continue
        # if is_most_common(j, aa_j, cnt_single):
        #     continue

        c_i = cnt_single[(i, aa_i)]
        c_j = cnt_single[(j, aa_j)]

        if c_i < 3 or c_j < 3:
            continue

        p_ij = c_ij / n_records
        p_i = c_i / n_records
        p_j = c_j / n_records

        p_j_cond_i = c_ij / c_i

        ratio = p_ij / (p_i * p_j)

        info.append((p_j_cond_i, ratio, i, j, aa_i, aa_j, c_i, c_j, c_ij))

    for i in range(8):
        print(n_aa_at_pos(i, cnt_single), file=sys.stderr)

    print(*sorted(info), sep="\n")

    show(info)


if __name__ == "__main__":
    main()
