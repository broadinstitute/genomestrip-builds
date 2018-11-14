#!/bin/bash

homologyFileList=$1
mergedFile=$2

if [ "${mergedFile}" == "" ]; then
    echo "Usage: <homologyFileList> <mergedFile>" && exit 1
fi

startField=$(zcat -f $(head -1 ${homologyFileList}) | head -1 | sed 's/\t/\n/g' | grep -n LEFTCHR | cut -d ':' -f1)

(zcat -f $(head -1 ${homologyFileList}) | head -1;
xargs -i sh -c 'zcat -f {} | tail -n +2' < ${homologyFileList}) \
| cut -f ${startField}- \
| awk -v FS="\t" -v OFS="\t" '{
    if (NR == 1) {
        if (substr($1, 1, 1) != "#") {
            $1 = "#" $1
        }
        print
        next
    }
    print | "sort -u -k1,1 -k2,2n -k3,6n"
}' \
| bgzip -c > ${mergedFile} || exit 1

tabix -s 1 -b 2 -e 6 ${mergedFile} || exit 1

