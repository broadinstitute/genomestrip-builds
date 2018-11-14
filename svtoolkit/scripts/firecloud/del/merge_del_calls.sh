#!/bin/bash

vcfFileList=$1
referenceBundle=$2

if [ "${referenceBundle}" == "" ]; then
    echo "Usage: <vcfFileList> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

mergedSitesVcf="gs_dels.sites.vcf.gz"

java -cp ${SV_CLASSPATH} -Xmx3g \
    org.broadinstitute.sv.apps.MergeDeletionSites \
    -R ${referenceFile} \
    -vcf ${vcfFileList} \
    -O ${mergedSitesVcf} \
    || exit 1
