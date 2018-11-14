#!/bin/bash

delVcf=$1
cnvVcf=$2
genderMapFile=$3
reclassifiedSitesVcf=$4
referenceBundle=$5

if [ "${referenceBundle}" == "" ]; then
    echo "Usage: <delVcf> <cnvVcf> <genderMapFile> <reclassifiedSitesVcf> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

# Index is needed to be able to query the vcf inside ReclassifyDeletions
echo "Indexing the cnv vcf file..."
bcftools index -t ${cnvVcf} || exit 1

java -cp ${SV_CLASSPATH} -Xmx3g \
    org.broadinstitute.sv.apps.ReclassifyDeletions \
    -R ${referenceFile} \
    -ploidyMapFile ${ploidyMapFile} \
    -delVcf ${delVcf} \
    -cnvVcf ${cnvVcf} \
    -genderMapFile ${genderMapFile} \
    -O ${reclassifiedSitesVcf} \
    || exit 1

