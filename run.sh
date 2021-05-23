#!/bin/bash
set -u

while getopts "f:i:o:d:a:" opt
do
    case "$opt" in
        f) filter_seq_file="$OPTARG";;
        i) input="$OPTARG";;
        o) output="$OPTARG";;
        d) result_dir="$OPTARG";;
        a) alignment="$OPTARG";;
        *) echo "Unknown parameter passed: $opt"; exit 1 ;;
    esac
done

filtered="$input".filtered

echo "Filtering..."

base="$(echo "$input" | cut -f 1 -d '.')"
ext="$(echo "$input" | cut -f 2 -d '.')"
filtered="$base"_filtered."$ext"

python filter.py \
    --filter_seq_file "$filter_seq_file" \
    --min_len 8 \
    --max_len 8 \
    < "$input" \
    > "$filtered" \
    2> /dev/null
echo "Filtering done"

# docker or use webform on https://www.ebi.ac.uk/Tools/msa/clustalo/
# docker image size ~400Mb

echo "Aligning..."
docker run -it --rm \
    -v $(pwd):/data \
    biocontainers/clustal-omega:v1.2.1_cv5 \
    clustalo \
        -i "$filtered" \
        -o "$alignment" \
        --force \
        --output-order input-order
echo "Aligned done!"


echo "Processing alignment..."
python process_aln.py \
    --out_dir "$result_dir" \
    --p_j_cond_i_cutoff 0 \
    < "$alignment" \
    > "$output"
    2> /dev/null
echo "Processing done!"

echo "Run done!"
