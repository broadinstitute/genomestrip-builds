#!/bin/bash

partitionName=$1
vcfFile=$2
referenceFile=$3

vcfextract -fields CHROM,POS,INFO:END ${vcfFile} \
| tail -n +2 \
| awk -v partitionName=${partitionName} -v OFS='\t' '{
    if (start[$1] == 0 || $2 < start[$1]) start[$1] = $2
    if (end[$1] == 0 || $3 > end[$1]) end[$1] = $3
} END {
    for (s in start) {
        print s, start[s], end[s], partitionName
    }
}' \
| joinfiles \
    <(cut -f1 ${referenceFile}.fai) 1 \
    - 1 \
| awk 'NF > 1' | cut -f2- > call_spans.dat || exit 1

