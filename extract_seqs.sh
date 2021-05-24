#!/bin/bash
set -u
export LC_NUMERIC="en_US.UTF-8"

lines=$(tail -n 10 "$1" | tac)

while IFS=, read -r p i aa_i j aa_j
do
    if ((i > j))
        then
            ndots=$((i - j - 1))
            pattern=".\{$j\}$aa_j.\{$ndots\}$aa_i"
        else
            ndots=$((j - i - 1))
            pattern=".\{$i\}$aa_i.\{$ndots\}$aa_j"
    fi
    printf "Position (%d %s) depends on position (%d %s) with p=%.2f\n" \
        $j $aa_j $i $aa_i $p
    grep "$pattern" --no-group-separator -B 1 "$2" | head -n $((2 * 3))
    echo


done <<< $(cut -f 1,3,4,5,6 -d ',' <<< "$lines")
