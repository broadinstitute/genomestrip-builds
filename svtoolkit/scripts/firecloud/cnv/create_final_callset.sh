#!/bin/bash

partitionName=$1
interval=$2
vcfFileList=$3
referenceBundle=$4

if [ -z "${referenceBundle}" ]; then
    echo "Usage: create_final_callset.sh <partitionName> <interval> <vcfFileList> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

runDir=cnv_callset
mkdir -p ${runDir}/merged/eval
mkdir -p ${runDir}/final

mergedVcf="${runDir}/merged/gs_cnv.genotypes.vcf.gz"

java -cp $SV_CLASSPATH -Xmx4g \
    org.broadinstitute.sv.apps.VCFMerge \
    -R $referenceFile \
    -vcf ${vcfFileList} \
    -O ${mergedVcf} \
    || exit 1

xargs rm < ${vcfFileList}

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -A GCContent \
    -A CopyNumberClass \
    -A NonVariant \
    -R ${referenceFile} \
    -vcf ${mergedVcf} \
    -writeReport true \
    -writeSummary true \
    -reportDirectory ${runDir}/merged/eval \
    || exit 1

awk 'NR > 1 && $8 != "NA" {print $1}' ${runDir}/merged/eval/CopyNumberClass.report.dat \
> ${runDir}/final/SelectedSites.list || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFFilter \
    -R ${referenceFile} \
    -vcf ${mergedVcf} \
    -includeSite ${runDir}/final/SelectedSites.list \
    -O ${runDir}/final/cnv.genotypes.vcf.gz \
    || exit 1

rm ${mergedVcf}

find ${runDir}/final | tar -cvzf cnv_callset.tar.gz --files-from -

