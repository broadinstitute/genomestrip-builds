#!/bin/bash

mdPath=$1
chrom=$2
profileFile=$3

if [ -z "${profileFile}" ]; then
    echo "Usage: <mdPath> <chrom> <profileFile>" && exit 1
fi

runDir=segpipe_output/
mkdir -p ${runDir} || exit 1

echo "Generating the index file..."
tabix -f -S 1 -s 2 -b 3 -e 4 ${profileFile} || exit 1
ls -l ${profileFile}*

SV_CLASSPATH=${SV_CLASSPATH}:${SV_DIR}/lib/SVToolkit-cnseg.jar
echo $SV_CLASSPATH
java -Xmx3g -cp java:${SV_CLASSPATH} \
    org.broadinstitute.sv.cnseg.CNVSegmentationPipeline \
    -md ${mdPath} \
    -L ${chrom} \
    -profileFile ${profileFile} \
    -block all \
    -referenceId HG38 \
    -OD ${runDir} \
    -P cnv.genotypingWindowSize:11 \
    || exit 1

# Delete all the currently unused files
find ${runDir} -type f ! -name cnSeg.tbl.gz -and ! -name dpSeg.tbl.gz -and ! -name .done -exec rm {} \;

tar -cvzf segpipe_output.tar.gz ${runDir} || exit 1

