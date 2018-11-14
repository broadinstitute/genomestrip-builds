#!/bin/bash

partitionBatchVcfList=$1

cut -f 1 ${partitionBatchVcfList} | sort -u > partitions.list

cat /dev/null > vcf_files.tsv
cut -f1 ${partitionBatchVcfList} | sort -u | sort -k1n \
| while read partition; do
    awk -v partition=$partition ' BEGIN {count=0} {
        if ($1 == partition) {
            if (count > 0) printf("\t")
            printf($3)
            count = count + 1
        }
    } END {printf("\n")}' ${partitionBatchVcfList} >> vcf_files.tsv
done

