#!/bin/bash

intervalFile=$1
vcfFileList=$2
sampleListFile=$3
filteredVcf=$4
referenceBundle=$5

if [ -z "${referenceBundle}" ]; then
    echo "Usage: filter_redundant_calls.sh <intervalFile> <vcfFileList> <sampleListFile> <filteredVcf> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

mkdir -p merged
mkdir -p deduped/eval

vcfFileName=$(basename ${filteredVcf})
mergedVcf=merged/${vcfFileName}
dedupedVcf=deduped/${vcfFileName}

echo $(date +"[%b %d %H:%M:%S] Merging/filtering vcf file(s)")
bcftools concat -a -Oz -R ${intervalFile} -o ${mergedVcf} $(cat ${vcfFileList} | xargs) || exit 1
bcftools index -t ${mergedVcf} || exit 1
echo "# variants to dedupe: "$(zcat ${mergedVcf} | grep -v '^#' | wc -l)

xargs rm < ${vcfFileList}

sampleArg=
if [ ! -z "${sampleListFile}" ]; then
    sampleArg="-sample ${sampleListFile}"
fi

echo $(date +"[%b %d %H:%M:%S] Filtering redundant sites")
java -cp $SV_CLASSPATH -Xmx12g \
    org.broadinstitute.sv.apps.FilterRedundantSites \
    -R $referenceFile \
    -vcf ${mergedVcf} \
    ${sampleArg} \
    -O ${dedupedVcf} \
    || exit 1
rm ${mergedVcf}

echo "# deduped variants: "$(zcat ${dedupedVcf} | grep -v '^#' | wc -l)

echo $(date +"[%b %d %H:%M:%S] Annotating the deduped vcf file")
java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -A GCContent \
    -A CopyNumberClass \
    -R ${referenceFile} \
    -vcf ${dedupedVcf} \
    -writeReport true \
    -writeSummary true \
    -reportDirectory deduped/eval \
    || exit 1

awk 'NR > 1 && $8 != "NA" {print $1}' deduped/eval/CopyNumberClass.report.dat \
> SelectedSites.list

echo $(date +"[%b %d %H:%M:%S] Creating the filtered vcf file")
bcftools view -i 'ID=@SelectedSites.list' ${dedupedVcf} -Oz -o ${filteredVcf}
bcftools index -t ${filteredVcf} || exit 1

echo "# filtered variants: "$(cat SelectedSites.list | wc -l)

echo $(date +"[%b %d %H:%M:%S] Done successfully")
