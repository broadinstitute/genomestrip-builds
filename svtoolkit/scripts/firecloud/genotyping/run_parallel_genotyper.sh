#!/bin/bash

vcfFile=$1
partitionName=$2
partitionArg=$3
mdPath=$4
bamFileList=$5
referenceBundle=$6
useLcMask=$7

if [ -z "${useLcMask}" ]; then
    echo "Usage: run_parallel_genotyper.sh <vcfFile> <partitionName> <partitionArg> <mdPath> <bamFileList> <referenceBundle> <useLcMask>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}
runDir=gtrun
mkdir -p ${runDir}/logs || exit 1

for attempt in {1..3}; do
    echo $(date +"[%b %d %H:%M:%S] Running ZipExtract, attempt=${attempt}")

java -Xmx2g -cp ${SV_CLASSPATH} \
    org.broadinstitute.sv.apps.ZipExtract \
    -I ${mdPath} \
    -file headers.bam \
    -file headers.bam.bai \
    -file sample_gender.report.txt

    rc=$?
    echo $(date +"[%b %d %H:%M:%S] ZipExtract completed, rc=${rc}")

    if [ "${rc}" -eq 0 ]; then
        break
    fi
done

if [ "${rc}" -ne 0 ]; then
    echo "Error running ZipExtract" && exit 1
fi

extraArgs=
if [ "${bamFileList}" == "NULL" ]; then
    bamFileList="headers.bam"
    extraArgs="-P genotyping.modules:depth"
fi

genomeMaskArg="-genomeMaskFile ${svMaskFile}"
if [ "${useLcMask}" == "true" ]; then
    genomeMaskArg="${genomeMaskArg} -genomeMaskFile ${lcMaskFile}"
fi

for attempt in {1..3}; do
    echo $(date +"[%b %d %H:%M:%S] Running genotyping for partition ${partitionName}, attempt=${attempt}")

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVGenotyper \
    -T SVGenotyperWalker \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R ${referenceFile} \
    -md ${mdPath} \
    ${genomeMaskArg} \
    -genderMapFile sample_gender.report.txt \
    -I ${bamFileList} \
    -disableGATKTraversal true \
    -vcf ${vcfFile} \
    -partitionName ${partitionName} \
    -partition ${partitionArg} \
    ${extraArgs} \
    -runDirectory ${runDir} \
    -O ${runDir}/${partitionName}.genotypes.vcf.gz

    rc=$?
    echo $(date +"[%b %d %H:%M:%S] Genotyping completed, rc=${rc}")

    if [ "${rc}" -eq 0 ]; then
        break
    fi
done

if [ "${rc}" -ne 0 ]; then
    echo "Error running genotyping" && exit 1
fi

echo "Archiving output..."

tar -cvzf gtrun.tar.gz ${runDir}
