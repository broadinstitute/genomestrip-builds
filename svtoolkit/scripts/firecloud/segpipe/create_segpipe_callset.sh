#!/bin/bash

callType=$1
genotypesVcf=$2
callsBedFile=$3
sampleListFile=$4
genderMapFile=$5
referenceBundle=$6

if [ -z "${referenceBundle}" ]; then
    echo "Usage: <callType> <genotypesVcf> <callsBedFile> <sampleListFile> <genderMapFile> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

mkdir metadata || exit 1
mkdir -p corrected/eval || exit 1

java -cp $SV_CLASSPATH -Xmx4g \
    org.broadinstitute.sv.apps.CorrectUncalledGenotypes \
    -R ${referenceFile} \
    -ploidyMapFile ${ploidyMapFile} \
    -genderMapFile ${genderMapFile} \
    -vcf ${genotypesVcf} \
    -callsFile ${callsBedFile} \
    -O corrected/segpipe_${callType}.genotypes.vcf.gz \
    -reportFile corrected/corrections.report.dat \
    || exit 1

sampleArg=
if [ ! -z "${sampleListFile}" ]; then
    sampleArg="-sample ${sampleListFile}"
fi

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -A GCContent \
    -A CopyNumberClass \
    -R ${referenceFile} \
    -md metadata \
    -genderMapFile ${genderMapFile} \
    -vcf corrected/segpipe_${callType}.genotypes.vcf.gz \
    ${sampleArg} \
    -writeReport true \
    -writeSummary true \
    -reportDirectory corrected/eval \
    || exit 1

awk 'NR > 1 && $8 != "NA" {print $1}' corrected/eval/CopyNumberClass.report.dat > SelectedSites.list

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFFilter \
    -R ${referenceFile} \
    -vcf corrected/segpipe_${callType}.genotypes.vcf.gz \
    -includeSite SelectedSites.list \
    -O segpipe_${callType}.genotypes.vcf.gz \
    || exit 1

