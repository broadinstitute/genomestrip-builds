#!/bin/bash

filePrefix=$1
runDir=$2
referenceFile=$3
repeatTrackFile=$4

# Filter calls in alpha-satellite regions
echo "Filtering the discovered variants..."

filtDir=${runDir}/filtering
mkdir -p ${filtDir} || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -R ${referenceFile} \
    -A MobileElements \
    -repeatTrackFile ${repeatTrackFile} \
    -vcf ${runDir}/${filePrefix}.sites.vcf.gz \
    -O ${filtDir}/${filePrefix}.annotated.sites.vcf \
    || exit 1

# Filter to select passing variants with acceptably low alpha-satellite fraction

vcfextract \
    -fields ID,FILTER,INFO:GSALPHASATFRACTION \
    ${filtDir}/${filePrefix}.annotated.sites.vcf \
    > ${filtDir}/${filePrefix}.info.dat || exit 1
cat ${filtDir}/${filePrefix}.info.dat | awk '$2 == "PASS" && $3 < 0.5 { print $1 }' > ${filtDir}/${filePrefix}.sites.list || exit 1


