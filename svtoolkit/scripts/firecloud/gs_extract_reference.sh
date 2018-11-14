#!/bin/bash

referenceBundle=$1

if [ -z "${referenceBundle}" ]; then
    echo "Usage: gs_extract_reference.sh <referenceBundle>" && exit 1
fi

echo $(date +"[%b %d %H:%M:%S] Extracting Genome STRiP reference bundle")
referenceName=$(tar -tf ${referenceBundle} | head -1 | cut -d / -f1)

tar -xvzf ${referenceBundle} || exit 1
rm -f ${referenceBundle}

referenceFile="${referenceName}/${referenceName}.fasta"
echo "Reference FASTA: "$referenceFile
 
source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1
