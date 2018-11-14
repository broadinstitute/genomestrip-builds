#!/bin/bash

partitionFile=$1
vcfSpansFile=$2
vcfFileList=$3
nullVcfPath=$4
outputFile=$5

bedtools intersect -a ${partitionFile} -b ${vcfSpansFile} -wa -wb -loj \
| joinfiles \
    - 8 \
    <(awk '{printf("%s\tP%04d\n", $0, NR)}' ${vcfFileList}) 2 \
| sort -u | sort -k4,4n \
| awk -v nullVcfPath=${nullVcfPath} 'BEGIN {
    part = 0
    isEmpty = 1
} {
    if ($4 != part) {
        if (part > 0) {
            if (isEmpty) {
                print nullVcfPath
            } else {
                printf("\n")
            }
        }
        part = $4
        isEmpty = 1
    }
    if ($9 != "") {
        if (!isEmpty) {
            printf("\t")
        }
        printf("%s", $9)
        isEmpty = 0
    }
} END {
    if (isEmpty) {
        print nullVcfPath
    } else {
        printf("\n")
    }
}' > ${outputFile}

