#!/bin/bash

intervalFile=$1
delVcf=$2
cnvVcf=$3
reclassifiedVcfFileList=$4
jointVcf=$5
referenceBundle=$6

if [ "${referenceBundle}" == "" ]; then
    echo "Usage: <intervalFile> <delVcf> <cnvVcf> <reclassifiedDelVcfList> <jointVcf> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

vcfFileName=$(basename ${jointVcf})
mergedVcf=merged/${vcfFileName}
dedupedVcf=deduped/${vcfFileName}

reclassifiedVcf="reclassified.genotypes.vcf.gz"

echo $(date +"[%b %d %H:%M:%S] Filtering reclassified vcf file(s)")
bcftools concat -a -Oz -R ${intervalFile} -o ${reclassifiedVcf} $(cat ${reclassifiedVcfFileList} | xargs) || exit 1
bcftools index -t ${reclassifiedVcf} || exit 1

xargs rm < ${reclassifiedVcfFileList}

vcfextract -fields ID ${delVcf} | tail -n +2 > del_sites.list
vcfextract -fields ID ${cnvVcf} | tail -n +2 > cnv_sites.list

echo "#dels: "$(cat del_sites.list | wc -l)
echo "#cnvs: "$(cat cnv_sites.list | wc -l)
echo "#reclassified: "$(zcat ${reclassifiedVcf} | grep -v '^#' | wc -l)

mkdir -p reclassified/eval
java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -A GCContent \
    -A CopyNumberClass \
    -R ${referenceFile} \
    -vcf ${reclassifiedVcf} \
    -writeReport true \
    -writeSummary true \
    -reportDirectory reclassified/eval \
    || exit 1

awk '$8 == "DUP" || $8 == "MIXED" {print $1}' reclassified/eval/CopyNumberClass.report.dat \
| joinfiles -select \
    <(vcfextract -fields ID,INFO:GSRECLASSIFIEDDEL ${reclassifiedVcf}) 1 \
    - 1 \
| tee >(cut -f1 > reclassified/include_cnvs.list) >(cut -f2 > reclassified/exclude_dels.list) > /dev/null

joinfiles \
    del_sites.list 1 \
    reclassified/exclude_dels.list 1 \
| awk 'NF == 1' \
| cat \
    - \
    cnv_sites.list \
    reclassified/include_cnvs.list \
> SelectedSites.list

echo $(date +"[%b %d %H:%M:%S] Merging del/cnv/reclassified vcf files")
java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.MergeVCFFiles \
    -R ${referenceFile} \
    -vcf ${delVcf} \
    -vcf ${cnvVcf} \
    -vcf ${reclassifiedVcf} \
    -includeSite SelectedSites.list \
    -O merged.genotypes.vcf.gz \
    || exit 1

#echo $(date +"[%b %d %H:%M:%S] Filtering redundant sites")
#echo "# variants to dedupe: "$(cat SelectedSites.list | wc -l)
#java -cp $SV_CLASSPATH -Xmx4g \
#    org.broadinstitute.sv.apps.FilterRedundantSites \
#    -R $referenceFile \
#    -vcf merged.genotypes.vcf.gz \
#    -O ${jointVcf} \
#    || exit 1

mv merged.genotypes.vcf.gz ${jointVcf}

echo "# variants: "$(zcat ${jointVcf} | grep -v '^#' | wc -l)
echo $(date +"[%b %d %H:%M:%S] Done successfully")

