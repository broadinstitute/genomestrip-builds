#!/bin/bash

archiveList=$1
referenceBundle=$2

if [ -z "${referenceBundle}" ]; then
    echo "Usage: gs_merge_metadata.sh <archiveList> <referenceBundle>" && exit 1
fi

source $(dirname $0)/gs_extract_reference.sh ${referenceBundle} || exit 1

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}
sampleMdDir=sample_metadata
mkdir -p ${sampleMdDir}

echo $(date +"[%b %d %H:%M:%S] Unarchiving per-sample metadata")

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/scripts/firecloud/UnarchiveData.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -R ${referenceFile} \
    -archive ${archiveList} \
    -deleteArchives true \
    -O ${sampleMdDir} \
    -tempDir ${SV_TMPDIR} \
    -jobLogDir ${SV_TMPDIR} \
    -jobRunner ParallelShell \
    -run \
    || exit 1

echo $(date +"[%b %d %H:%M:%S] Merging metadata")
echo "Disk space statistics"
df -h .

mdInputLocations=$(find ${sampleMdDir} -mindepth 1 -maxdepth 1 | awk '{print "-mdi " $0}' | xargs)
echo "mdInputLocations=${mdInputLocations}"

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/preprocess/SVMergeMetadata.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -jobLogDir ${profilesDir}/logs \
    -R ${referenceFile} \
    ${mdInputLocations} \
    -md metadata \
    -useParallelMerge true \
    -jobRunner ParallelShell \
    -run \
    -l DEBUG \
    || exit 1

echo "Disk space statistics"
df -h .

echo $(date +"[%b %d %H:%M:%S] Deleting per-sample metadata to reclaim space")
rm -rf ${mdInputLocations}

echo "Disk space statistics"
df -h .

echo $(date +"[%b %d %H:%M:%S] Zipping metadata directory")

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.ZipMetadata \
    -md metadata \
    -O metadata.zip \
    || exit 1

echo $(date +"[%b %d %H:%M:%S] Merging per-sample metadata completed successfully")

