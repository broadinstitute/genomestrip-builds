#!/bin/bash

vcfFile=$1
mdPath=$2
bamFileList=$3
outDir=$4
referenceBundle=$5
credentialsKeyFile=$6

if [ -z "${credentialsKeyFile}" ]; then
    echo "Usage: run_parallel_genotyper.sh <<vcfFile> <partitionName> <partitionArg> <mdPath> <bamFileList> <referenceBundle> <useLcMask> <credentialsKeyFile>" && exit 1
fi

mkdir -p ${outDir}/logs || exit 1

genomeMaskArg="-genomeMaskFile ${svMaskFile} -genomeMaskFile ${lcMaskFile}"

extraArgs=
if [ "${depthOnly}" == "true" ]; then
    extraArgs="-P genotyping.modules:depth"
fi

outFile="${outDir}/"`echo ${vcfFile} | awk -F / '{ print $NF}' | sed 's/.sites.vcf.gz$/.vcf.gz/' | sed 's/.genotypes.vcf.gz$/.vcf.gz/' | sed 's/.vcf.gz$/.genotypes.vcf.gz/'`

for attempt in {1..5}; do

    echo $(date +"[%b %d %H:%M:%S] Running SVGenotyper, attempt=${attempt}")

    java -cp ${SV_CLASSPATH} -Xmx4g \
        org.broadinstitute.gatk.queue.QCommandLine \
        -S ${SV_DIR}/qscript/SVGenotyper2.q \
        -S ${SV_DIR}/qscript/SVQScript.q \
        -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
        -cp ${SV_CLASSPATH} \
        -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
        -jobLogDir ${outDir}/logs \
        -tempDir ${SV_TMPDIR} \
        -R ${referenceFile} \
        -md ${mdPath} \
        ${genomeMaskArg} \
        -genderMapFile sample_gender.report.txt \
        -I ${bamFileList} \
        -parallelRecords 1000 \
        ${extraArgs} \
        -outDirectory ${outDir} \
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
