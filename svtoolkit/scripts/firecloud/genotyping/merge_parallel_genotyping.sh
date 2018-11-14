#!/bin/bash

filePrefix=$1
partitionsArchiveList=$2
referenceBundle=$3

if [ -z "${referenceBundle}" ]; then
    echo "Usage: merge_parallel_genotyping.sh <filePrefix> <partitionsArchiveList> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}
runDir=gtrun
mkdir -p ${runDir}/logs || exit 1

awk '{print "tar -xvzf " $0 " -C gtrun --strip-components 1"}' ${partitionsArchiveList} | sh

java -cp $SV_CLASSPATH -Xmx4g \
    org.broadinstitute.sv.apps.MergeGenotyperOutput \
    -R $referenceFile \
    -runDirectory gtrun \
    -O ${runDir}/${filePrefix}.genotypes.vcf.gz \
    || exit 1

find ${runDir} -maxdepth 1 -type f -name 'P[0-9]*.genotypes.*' -exec rm {} \;

tar -cvzf gtrun.tar.gz ${runDir} || exit 1

