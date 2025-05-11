#!/bin/bash

# output file
OUTPUT="all_results.csv"

# write csv header
echo "Threshold,MCC,E-value source,Test set" > "$OUTPUT"

# function to process a classification set (either set1 or set2)
process_set() {
    local setfile=$1
    local setname=$2

    for i in $(seq 1 9); do
        evalue="1e-$i"
        python3 performance.py "$setfile" "$evalue"
    done | grep 'threshold' | while read -r line; do

        threshold=$(echo "$line" | grep -oP 'threshold=\s*\K\S+')
        mcc=$(echo "$line" | grep -oP 'MCC=\s*\K\S+')
        fullseq=$(echo "$line" | grep -oP 'fullseq=\s*\K\S+')

        if [[ "$fullseq" == "True" ]]; then
            evalue_type="full-seq"
        else
            evalue_type="single-dom"
        fi

        formatted_threshold=$(printf "%.0E" "$threshold")

    # print each result row to stdout (later sorted and appended to file)
        echo "$formatted_threshold,$mcc,$evalue_type,$setname"
    done
}

# process both sets and collect all results
#change these to the correct set files name that you have and that you want to analyze (MSA or MStA alignments results)
# _strali.class is the MStA alignment results
# _seqali.class is the MSA alignment results
results=$( 
    # process_set set_1_strali.class set1
    # process_set set_2_strali.class set2
    process_set set_1_seqali.class set1
    process_set set_1_seqali.class set2
)

# sort results by test set, then by e-value source, then by threshold
echo "$results" | sort -t',' -k4,4 -k3,3 -k1,1 >> "$OUTPUT"
