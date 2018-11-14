#!/bin/bash

filePrefix=$1
partitionsArchiveList=$2
storeReadPairFile=$3
referenceBundle=$4
repeatTrackFile=$5

if [ -z "${repeatTrackFile}" ]; then
    echo "Usage: merge_parallel_discovery.sh <filePrefix> <partitionsArchiveList> <storeReadPairFile> <referenceBundle> <repeatTrackFile>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
rm -f ${referenceBundle}

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}
runDir=del_output
mkdir -p del_input
mkdir -p ${runDir}/logs || exit 1

awk '{print "tar -xvzf " $0 " -C del_input --strip-components 1; rm -f " $0}' ${partitionsArchiveList} | sh

# Temporary fix needed to substitute spaces with underscores in the sample names
find del_input -name '*.vcf.gz' -exec $(dirname $0)/fix_vcfs_with_spaces.sh {} \;

java -cp $SV_CLASSPATH -Xmx4g \
    org.broadinstitute.sv.apps.MergeDiscoveryOutput \
    -R $referenceFile \
    -runDirectory del_input \
    -O ${runDir}/${filePrefix}.sites.unfiltered.vcf.gz \
    || exit 1

$(dirname $0)/apply_default_filters.sh \
    ${runDir}/${filePrefix}.sites.unfiltered.vcf.gz \
    ${runDir}/${filePrefix}.sites.vcf.gz \
    ${referenceFile} \
    || exit 1

if [ "${storeReadPairFile}" == "true" ]; then
    samtools merge -c -p ${runDir}/gs_dels.pairs.bam `find del_input -name '*.pairs.bam' | xargs` || exit 1
    samtools index ${runDir}/gs_dels.pairs.bam || exit 1
fi

rm -rf del_input

$(dirname $0)/filter_discovered_sites.sh \
    ${filePrefix} \
    ${runDir} \
    ${referenceFile} \
    ${repeatTrackFile} \
    || exit 1

gtDir=${runDir}/genotyping
mkdir -p ${gtDir} || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFFilter \
    -R ${referenceFile} \
    -vcf ${runDir}/${filePrefix}.sites.vcf.gz \
    -includeSite ${runDir}/filtering/${filePrefix}.sites.list \
    -O ${gtDir}/${filePrefix}.sites.vcf.gz \
    || exit 1

tar -cvzf del_output.tar.gz ${runDir} || exit 1

