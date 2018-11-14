#!/bin/bash

chromsListFile=$1
windowSize=$2
intervalsListFile=$3
referenceBundle=$4

if [ -z "${referenceBundle}" ]; then
    echo "Usage: create_genome_partitions.sh <chromsListFile> <windowSize> <intervalsListFile> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

joinfiles -select \
    ${referenceFile}.fai 1 \
    ${chromsListFile} 1 \
| awk -v windowSize=${windowSize} '{
    startPos = 1
    endPos = 0
    while (endPos < $2) {
        endPos = startPos + windowSize - 1
        if (endPos > $2) endPos = $2
        print $1":"startPos"-"endPos

        startPos = endPos + 1
    }
}' > ${intervalsListFile}

