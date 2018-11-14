#!/bin/bash

pairsType=$1
pairsFileList=$2
mergedFile=$3

if [ "${mergedFile}" == "" ]; then
    echo "Usage: <pairsType> <pairsFileList> <mergedFile>" && exit 1
fi

if [ "${pairsType}" == "discovery" ]; then
    sortFields="-k2,2 -k8,8n -k13,13n"
    tabixArgs="-s 2 -b 8 -e 13"
else
    sortFields="-k8,8 -k9,9n -k13,13n"
    tabixArgs="-s 8 -b 9 -e 13"
fi

(head -1 $(head -1 ${pairsFileList});
xargs -i sh -c 'zcat -f {} | tail -n +2' < ${pairsFileList} \
| sort -u -t$'\t' $sortFields) \
| awk -v FS="\t" -v OFS="\t" '{
    if (NR == 1) {
        $1 = "#" $1
    }
    print
}' \
| bgzip -c > ${mergedFile} || exit 1

tabix ${tabixArgs} ${mergedFile}

