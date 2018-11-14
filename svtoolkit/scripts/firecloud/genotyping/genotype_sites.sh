#!/bin/bash

vcfFile=$1
mdPath=$2
bamFileList=$3
useLcMask=$4
runDir=$5
referenceFile=$6

if [ "${referenceFile}" == "" ]; then
    echo "Usage: genotype_sites.sh <vcfFile> <mdPath> <bamFileList> <useLcMask> <runDir> <referenceFile>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}

genomeMaskArg="-genomeMaskFile ${svMaskFile}"
if [ "${useLcMask}" == "true" ]; then
    genomeMaskArg="${genomeMaskArg} -genomeMaskFile ${lcMaskFile}"
fi

outFile="${runDir}/"`echo ${vcfFile} | awk -F / '{ print $NF}' | sed 's/.sites.vcf.gz$/.vcf.gz/' | sed 's/.genotypes.vcf.gz$/.vcf.gz/' | sed 's/.vcf.gz$/.genotypes.vcf.gz/'`

for attempt in {1..5}; do

    echo $(date +"[%b %d %H:%M:%S] Running SVGenotyper, attempt=${attempt}")

    java -cp ${SV_CLASSPATH} -Xmx4g \
        org.broadinstitute.gatk.queue.QCommandLine \
        -S ${SV_DIR}/qscript/SVGenotyper2.q \
        -S ${SV_DIR}/qscript/SVQScript.q \
        -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
        -cp ${SV_CLASSPATH} \
        -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
        -jobLogDir ${runDir}/logs \
        -tempDir ${SV_TMPDIR} \
        -R ${referenceFile} \
        -md ${mdPath} \
        ${genomeMaskArg} \
        -genderMapFile sample_gender.report.txt \
        -I ${bamFileList} \
        -parallelRecords 1000 \
        ${extraArgs} \
        -runDirectory ${runDir} \
        -vcf ${vcfFile} \
        -skipAnnotator CNQuality \
        -O ${outFile} \
        -jobRunner ParallelShell \
        -gatkJobRunner ParallelShell \
        -run

    rc=$?
    if [ "${rc}" -eq 0 ]; then
        break
    fi
done

if [ "${rc}" -ne 0 ]; then
    echo "Error running SVGenotyper" && exit ${rc}
fi

echo "SVGenotyper completed successfully"
