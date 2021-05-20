import sys

from Bio import SeqIO


def main():
    seq_set = set(
        [
            "EDTFKHLG",
            "SDYASNHA",
            "DNFASNYA",
            "EDDDWSPHWVNPAPEHY",
            "RVSVWTSE",
            "LS",
            "EGIYHW",
            "GYNIHYHG",
        ]
        + [
            "WIQVYRQL",
            "IHLTRHY",
            "WIHLTNYV",
            "GYY",
            "WNTFTSYY",
            "WIHLTGYY",
            "SYA",
            "SSYTI",
            "FSSYA",
            "GYTLPTT",
            "TGYY",
            "SYD",
            "VFTSYA",
        ]
        + ["VTPLPVMV", "GDSFKINT"]
    )

    # seq_set = set()

    for record in SeqIO.parse(sys.stdin, "fasta"):
        if record.seq in seq_set or len(record.seq) != 8:
            print(f"skipping {record.id}", file=sys.stderr)
        else:
            SeqIO.write([record], sys.stdout, "fasta")
            seq_set.add(record.seq)
    print(len(seq_set), file=sys.stderr)


if __name__ == "__main__":
    main()
