#!/bin/bash

mdPath=$1
segdupFile=$2
sequenceMapFile=$3
referenceBundle=$4

if [ -z "${referenceBundle}" ]; then
    echo "Usage: run_segdup_pipeline.sh <mdPath> <segdupFile> <referenceBundle> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1

runDir=segdup_output
mkdir -p ${runDir}
SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}

java -Xmx2g -cp ${SV_CLASSPATH} \
    org.broadinstitute.sv.apps.ZipExtract \
    -I ${mdPath} \
    -file headers.bam \
    -file headers.bam.bai \
    -file sample_gender.report.txt \
    || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    -Djava.io.tmpdir=${SV_TMPDIR} \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/discovery/segdup/SegDupDiscoveryPipeline.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -R ${referenceFile} \
    -md ${mdPath} \
    -genderMapFile sample_gender.report.txt \
    -segdupFile ${segdupFile} \
    -sequenceMapFile ${sequenceMapFile} \
    -I headers.bam \
    -runDirectory ${runDir} \
    -jobLogDir ${runDir}/logs \
    -jobRunner Shell \
    -gatkJobRunner Shell \
    -run

tar -cvzf segdup_output.tar.gz ${runDir}

