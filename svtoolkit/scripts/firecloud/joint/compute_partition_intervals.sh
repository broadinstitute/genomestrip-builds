#!/bin/bash

partitionFile=$1
outputFile=$2

awk 'BEGIN {
    part = 0
} {
    if ($4 != part) {
        if (part > 0) printf("\n")
        part = $4
    } else {
        printf("\t")
    }
    printf("%s:%d-%d", $1, $2, $3)
} END {
    printf("\n")
}' ${partitionFile} > ${outputFile}

