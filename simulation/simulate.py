import argparse
import enum
import sys
from enum import Enum, auto
from typing import List, Tuple

import numpy as np
from Bio import Seq, SeqIO


class AminoAcids(Enum):
    A = auto()  # Alanine (Ala)
    C = auto()  # Cysteine (Cys)
    D = auto()  # Aspartic Acid (Asp)
    E = auto()  # Glutamic Acid (Glu)
    F = auto()  # Phenylalanine (Phe)
    G = auto()  # Glycine (Gly)
    H = auto()  # Histidine (His)
    I = auto()  # Isoleucine (Ile)
    K = auto()  # Lysine (Lys)
    L = auto()  # Leucine (Leu)
    M = auto()  # Methionine (Met)
    N = auto()  # Asparagine (Asn)
    P = auto()  # Proline (Pro)
    Q = auto()  # Glutamine (Gln)
    R = auto()  # Arginine (Arg)
    S = auto()  # Serine (Ser)
    T = auto()  # Threonine (Thr)
    V = auto()  # Valine (Val)
    W = auto()  # Tryptophan (Trp)
    Y = auto()  # Tyrosine (Tyr)

    @staticmethod
    def to_tuple():
        return tuple(c.name for c in AminoAcids)

    @staticmethod
    def wo_aa(aa):
        return tuple(c.name for c in AminoAcids if c != aa)


def generate_pairs(
    aa_i: AminoAcids,
    aa_j: AminoAcids,
    p_i: float,  # condition freq
    p_ji: float,  # dependency
    n_seq: int,  # number of samples
) -> List[str]:

    """Generates sample of two letters, first dependent on second p(j|i)"""

    n_aa_i = np.random.binomial(n_seq, p_i)
    n_aa_ji = np.random.binomial(n_aa_i, p_ji)

    pos1 = np.random.choice(AminoAcids.wo_aa(aa_i), size=n_seq, replace=True)
    pos1[:n_aa_i] = aa_i.name

    pos2 = np.random.choice(AminoAcids.wo_aa(aa_j), size=n_seq, replace=True)
    pos2[:n_aa_ji] = aa_j.name

    print(n_aa_i, file=sys.stderr)
    print(n_aa_ji, file=sys.stderr)
    # print(file=sys.stderr)

    # random shuffle to prevent dependency on different pair
    res = [x + y for x, y in zip(pos1, pos2)]
    np.random.shuffle(res)
    cnt = 0
    for s in res:
        if s == aa_i.name + aa_j.name:
            cnt += 1

    print(cnt, file=sys.stderr)
    print(n_aa_i / n_seq, file=sys.stderr)
    print(n_aa_ji / n_seq, file=sys.stderr)
    print(n_aa_ji / n_aa_i, file=sys.stderr)
    print(aa_i.name, aa_j.name, file=sys.stderr)
    print(file=sys.stderr)

    return res


def generate_sample(
    p_i_range: Tuple[float, float],
    p_ji_range: Tuple[float, float],
    n_seq: int,
    n_pairs: int,
) -> List[str]:

    """Generates samples with dependent positions (1 depends on 0, 3 depends on 2...)"""

    assert 2 * n_pairs < len(
        AminoAcids
    ), f"Too many pairs, max value is {len(AminoAcids) // 2}"

    aa = iter(AminoAcids)
    pairs = []
    for _ in range(n_pairs):
        aa_i = next(aa)
        aa_j = next(aa)
        p_i = np.random.uniform(*p_i_range)
        p_ji = np.random.uniform(*p_ji_range)

        pairs.append(generate_pairs(aa_i, aa_j, p_i, p_ji, n_seq))

    return ["".join(seq) for seq in zip(*pairs)]


def random_sample(n_seq: int, n_aa: int) -> List[str]:

    """Generate random strings"""

    return [
        "".join(s)
        for s in np.random.choice(AminoAcids.to_tuple(), size=(n_seq, n_aa))
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p_i_range", nargs=2, type=float, required=True)
    parser.add_argument("--p_ji_range", nargs=2, type=float, required=True)
    parser.add_argument("--n_seq", type=int, required=True)
    parser.add_argument("--n_pairs", type=int, required=True)
    parser.add_argument("--n_randoms", type=int, required=True)
    parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()

    np.random.seed(args.seed)

    dependent = generate_sample(
        args.p_i_range, args.p_ji_range, args.n_seq, args.n_pairs
    )
    independent = random_sample(args.n_seq, args.n_randoms)

    seqs = [x + y for x, y in zip(dependent, independent)]

    records = [
        SeqIO.SeqRecord(
            seq=Seq.Seq(seq),
            id=str(i),
            description="",
        )
        for i, seq in enumerate(seqs)
    ]

    SeqIO.write(records, sys.stdout, "fasta")


if __name__ == "__main__":
    main()
