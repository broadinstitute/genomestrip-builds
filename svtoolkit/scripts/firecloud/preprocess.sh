#!/bin/bash

sampleId=$1
bamFileList=$2
referenceBundle=$3

if [ -z "${referenceBundle}" ]; then
    echo "Usage: gs_preprocess.sh <sampleId> <bamFileList>> <referenceBundle>" && exit 1
fi

echo $(date +"[%b %d %H:%M:%S] Extracting Genome STRiP reference bundle")
referenceName=$(tar -tf ${referenceBundle} | head -1 | cut -d / -f1)
tar -xvzf ${referenceBundle}
referenceFile="${referenceName}/${referenceName}.fasta"
echo "Reference FASTA: "$referenceFile
        
mkdir tmp
mkdir metadata
mkdir metadata/logs

echo $(date +"[%b %d %H:%M:%S] Starting preprocessing")

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVPreprocess.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -jobLogDir metadata/logs \
    -R ${referenceFile} \
    -md metadata \
    -I ${bamFileList} \
    -deleteIntermediateDirs true \
    -jobRunner ParallelShell \
    -run \
    || exit 1

echo $(date +"[%b %d %H:%M:%S] Zipping metadata directory")
java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.ZipMetadata \
    -md metadata \
    -O metadata.zip \
    || exit 1

echo $(date +"[%b %d %H:%M:%S] Preprocessing completed successfully")

