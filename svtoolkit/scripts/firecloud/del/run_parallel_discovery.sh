#!/bin/bash

partitionName=$1
partitionArgs="$2"
storeReadPairFile=$3
mdPath=$4
bamFileList=$5
referenceBundle=$6
repeatTrackFile=$7

if [ -z "${repeatTrackFile}" ]; then
    echo "Usage: run_parallel_discovery.sh <partitionName> <partitionArgs> <storeReadPairFile> <mdPath> <bamFileList> <referenceBundle> <repeatTrackFile>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}
runDir=del_output
mkdir -p ${runDir}/logs || exit 1

java -Xmx2g -cp ${SV_CLASSPATH} \
    org.broadinstitute.sv.apps.ZipExtract \
    -I ${mdPath} \
    -file sample_gender.report.txt \
    || exit 1

# Avoid problem with certain read pairs - this should be debugged
extraArgs="-P select.validateReadPairs:false"

for attempt in {1..3}; do

    echo $(date +"[%b %d %H:%M:%S] Running Deletion discovery for partition ${partitionName}, attempt=${attempt}")

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVDiscovery \
    -T SVDiscoveryWalker \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R ${referenceFile} \
    -md ${mdPath} \
    -genomeMaskFile ${svMaskFile} \
    -genderMapFile sample_gender.report.txt \
    -I ${bamFileList} \
    -disableGATKTraversal true \
    -partitionName ${partitionName} \
    ${partitionArgs} \
    -runFilePrefix ${partitionName} \
    -storeReadPairFile ${storeReadPairFile} \
    ${extraArgs} \
    -runDirectory ${runDir} \
    -O ${runDir}/${partitionName}.discovery.vcf.gz

    rc=$?
    echo $(date +"[%b %d %H:%M:%S] Deletion discovery completed, rc=${rc}")

    if [ "${rc}" -eq 0 ]; then
        break
    fi
done

if [ "${rc}" -ne 0 ]; then
    echo "Error running Deletion discovery" && exit 1
fi

echo "Archiving output..."
tar -cvzf del_output.tar.gz del_output

