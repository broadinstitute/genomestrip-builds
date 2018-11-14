#!/bin/bash

# It's assumed that inputFile has been sorted by groupingField
inputFile=$1
targetField=$2
groupingField=$3

awk -F '\t' -v targetField=${targetField} -v groupingField=${groupingField} 'BEGIN {
    group = ""
} {
    if ($groupingField != group) {
        if (group != "") printf("\n")
        group = $groupingField
    } else {
        printf("\t")
    }
    targetValue = $targetField
    printf("%s", targetValue)
} END {
    printf("\n")
}' ${inputFile}

