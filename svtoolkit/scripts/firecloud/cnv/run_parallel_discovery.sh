#!/bin/bash

partitionName=$1
mdPath=$2
intervalListFile=$3
referenceBundle=$4
credentialsKeyFile=$5

if [ -z "${credentialsKeyFile}" ]; then
    echo "Usage: run_parallel_discovery.sh <partitionName> <mdPath> <intervalListFile> <referenceBundle> <credentialsKeyFile>" && exit 1
fi

echo "Exporting GOOGLE_APPLICATION_CREDENTIALS..."
export GOOGLE_APPLICATION_CREDENTIALS=${credentialsKeyFile}

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1

genomeMaskArg="-genomeMaskFile ${svMaskFile} -genomeMaskFile ${lcMaskFile}"

bamFileList=headers.bam

java -Xmx2g -cp ${SV_CLASSPATH} \
    org.broadinstitute.sv.apps.ZipExtract \
    -I ${mdPath} \
    -file headers.bam \
    -file headers.bam.bai \
    -file sample_gender.report.txt \
    || exit 1

runDir=cnv_output
mkdir -p ${runDir}
SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}

tilingWindowSize=1000
tilingWindowOverlap=500
maximumReferenceGapLength=1000
boundaryPrecision=100
minimumRefinedLength=500

for attempt in {1..2}; do

    echo $(date +"[%b %d %H:%M:%S] Running CNV pipeline, attempt=${attempt}")

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/discovery/cnv/CNVDiscoveryPipeline.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -tempDir ${SV_TMPDIR} \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -R ${referenceFile} \
    -md ${mdPath} \
    -genderMapFile sample_gender.report.txt \
    -I headers.bam \
    -intervalList ${intervalListFile} \
    -tilingWindowSize ${tilingWindowSize} \
    -tilingWindowOverlap ${tilingWindowOverlap} \
    -maximumReferenceGapLength ${maximumReferenceGapLength} \
    -boundaryPrecision ${boundaryPrecision} \
    -minimumRefinedLength ${minimumRefinedLength} \
    -useAllSamplesInDiscovery \
    -lastStage 11 \
    --disableJobReport \
    -runDirectory ${runDir} \
    -jobLogDir ${runDir}/logs \
    -jobRunner ParallelShell \
    -gatkJobRunner ParallelShell \
    -run

    rc=$?
    echo $(date +"[%b %d %H:%M:%S] Deletion pipeline completed, rc=${rc}")

    if [ "${rc}" -eq 0 ]; then
        break
    fi
done

tar -cvzf cnv_output.tar.gz ${runDir}

if [ "${rc}" -eq 0 ]; then
    echo "CNV pipeline completed successfully"
else
    echo "Error running CNV pipeline"
fi
exit ${rc}

#tar -C ${runDir} -cvzf cnv_output.tar.gz cnv_stage11 || exit 1
