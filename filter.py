import argparse
import sys

from Bio import SeqIO


def main() -> None:

    parser = argparse.ArgumentParser()
    parser.add_argument("--filter_seq_file", type=str, nargs="?", default=None)
    parser.add_argument("--min_len", type=int, default=0)
    parser.add_argument("--max_len", type=int, default=sys.maxsize)

    args = parser.parse_args()

    if args.filter_seq_file is not None:
        with open(args.filter_seq_file, "r") as f:
            seq_set = set(map(lambda s: s.strip(), f.readlines()))
    else:
        seq_set = set()

    for record in SeqIO.parse(sys.stdin, "fasta"):
        if (
            record.seq in seq_set
            or not args.min_len <= len(record.seq) <= args.max_len
        ):
            print(f"skipping {record.id}", file=sys.stderr)
        else:
            SeqIO.write([record], sys.stdout, "fasta")
            seq_set.add(record.seq)
    print("Number of filtered seqs:", len(seq_set), file=sys.stderr)


if __name__ == "__main__":
    main()
