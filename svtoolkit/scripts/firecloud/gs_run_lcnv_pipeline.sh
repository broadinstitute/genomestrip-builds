#!/bin/bash

sampleList=$1
profilesArchive=$2
referenceBundle=$3

if [ -z "${referenceBundle}" ]; then
    echo "Usage: gs_run_lcnv_pipeline.sh <sampleList> <profilesArchive> <referenceBundle>" && exit 1
fi

source $(dirname $0)/gs_extract_reference.sh ${referenceBundle} || exit 1

sampleListArg=
fileLength=$(cat $sampleList | wc -c)
if [ "${fileLength}" -gt 1 ]; then
    sampleListArg="-sample ${sampleList}"
fi

profilesDir=profiles
mkdir -p ${profilesDir} || exit 1

tar -xvzf ${profilesArchive} --strip-components 1 -C ${profilesDir} || exit 1

genderMapFile=${profilesDir}/sample_gender.report.txt
runDir=lcnv_calls
logDir=${runDir}/logs
mkdir -p ${logDir} || exit 1

# Progressive scan
maxDepth=50
binSize=10000

echo $(date +"[%b %d %H:%M:%S] Running LCNV pipeline")

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/discovery/lcnv/LCNVDiscoveryPipeline.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -jobLogDir ${logDir} \
    -R ${referenceFile} \
    -profilesDir ${profilesDir} \
    -maxDepth ${maxDepth} \
    -perSampleScan true \
    ${sampleListArg} \
    -genderMapFile ${genderMapFile} \
    -runDirectory ${runDir} \
    -jobRunner ParallelShell \
    -run \
    || exit 1

echo $(date +"[%b %d %H:%M:%S] Archiving lcnv_calls directory")
tar -cvzf lcnv_calls.tar.gz ${runDir}

echo $(date +"[%b %d %H:%M:%S] LCNV pipeline completed successfully")
