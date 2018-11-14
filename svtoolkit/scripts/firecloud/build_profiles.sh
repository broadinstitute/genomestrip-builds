#!/bin/bash

chrom=$1
binSize=$2
mdPath=$3
referenceBundle=$4
useLcMask=$5

if [ -z "${referenceBundle}" ]; then
    echo "Usage: build_profiles.sh <chrom> <binSize> <mdPath> <referenceBundle> <useLcMask>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

mkdir tmp
profilesDir=profiles_${binSize}
mkdir -p ${profilesDir}/logs

genomeMaskArg="-genomeMaskFile ${svMaskFile}"
if [ "${useLcMask}" == "true" ]; then
    genomeMaskArg="${genomeMaskArg} -genomeMaskFile ${lcMaskFile}"
fi

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/profiles/GenerateDepthProfiles.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -jobLogDir ${profilesDir}/logs \
    -R ${referenceFile} \
    -L ${chrom} \
    -profileBinSize ${binSize} \
    -md ${mdPath} \
    ${genomeMaskArg} \
    -runDirectory ${profilesDir} \
    -jobRunner ParallelShell \
    -run \
    || exit 1

echo $(date +"[%b %d %H:%M:%S] Archiving profiles directory")
tar -cvzf profiles_${binSize}.tar.gz ${profilesDir} || exit 1
echo $(date +"[%b %d %H:%M:%S] Building profiles completed successfully")

