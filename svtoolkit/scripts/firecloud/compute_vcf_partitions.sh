#!/bin/bash

vcfFile=$1
parallelRecords=$2
partitionsFile=$3

if [ -z "${partitionsFile}" ]; then
    echo "Usage: compute_vcf_partitions.sh <vcfFile> <parallelRecords> <partitionsFile>" && exit 1
fi

java -Xmx2g \
    -cp $SV_CLASSPATH \
    org.broadinstitute.sv.apps.ComputeVCFPartitions \
    -I ${vcfFile} \
    -parallelRecords ${parallelRecords} \
    -O ${partitionsFile} \
    || exit 1

