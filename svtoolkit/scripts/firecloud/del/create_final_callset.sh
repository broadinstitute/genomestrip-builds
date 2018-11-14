#!/bin/bash

vcfFileList=$1
referenceBundle=$2

if [ -z "${referenceBundle}" ]; then
    echo "Usage: create_final_callset.sh <vcfFileList> <referenceBundle>" && exit 1
fi

if [[ $referenceBundle =~ .fasta$ ]]; then
    referenceFile=${referenceBundle}
    source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1
else
    source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
fi

runDir=del_callset
mkdir -p ${runDir}/merged/eval
mkdir -p ${runDir}/final/eval

mergedVcf="${runDir}/merged/del.genotypes.vcf.gz"

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFMerge \
    -R ${referenceFile} \
    -vcf ${vcfFileList} \
    -includeInfoTag END \
    -includeInfoTag GSELENGTH \
    -includeInfoTag SVTYPE \
    -includeInfoTag SVLEN \
    -includeInfoTag CIPOS \
    -includeInfoTag CIEND \
    -includeInfoTag GSBKPT \
    -O ${mergedVcf} \
    || exit 1

#xargs rm < ${vcfFileList}

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
> ${runDir}/final/SelectedSites.list

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFFilter \
    -R ${referenceFile} \
    -vcf ${mergedVcf} \
    -includeSite ${runDir}/final/SelectedSites.list \
    -O ${runDir}/final/del.genotypes.vcf.gz \
    || exit 1

rm ${mergedVcf}

tar -cvzf del_callset.tar.gz ${runDir}/final

