#!/bin/bash

batchGenotypesVcfList=$1
referenceBundle=$2

if [ -z "${referenceBundle}" ]; then
    echo "Usage: create_final_callset.sh <batchGenotypesVcfList> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

vcfArg=$(awk '{print "-vcf " $0}' ${batchGenotypesVcfList} | xargs)

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFMerge \
    -R ${referenceFile} \
    ${vcfArg} \
    -includeInfoTag END \
    -includeInfoTag GSELENGTH \
    -includeInfoTag SVTYPE \
    -includeInfoTag GSRECLASSIFIEDDEL \
    -O reclassified.genotypes.vcf.gz \
    || exit 1

