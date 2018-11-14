#!/bin/bash

mdPath=$1
bamFileList=$2
filePrefix=$3
runDir=$4
referenceFile=$5
repeatTrackFile=$6

if [ -z "${repeatTrackFile}" ]; then
    echo "Usage: post_deletion.sh <mdPath> <bamFileList> <filePrefix> <runDir> <referenceFile> <repeatTrackFile>" && exit 1
fi

SV_TMPDIR=tmp

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

java -Xmx1g -cp ${SV_CLASSPATH} org.broadinstitute.sv.apps.VCFExtract \
    -fields ID,FILTER,INFO:GSALPHASATFRACTION \
    ${filtDir}/${filePrefix}.annotated.sites.vcf \
    > ${filtDir}/${filePrefix}.info.dat || exit 1
cat ${filtDir}/${filePrefix}.info.dat | awk '$2 == "PASS" && $3 < 0.5 { print $1 }' > ${filtDir}/${filePrefix}.sites.list || exit 1

gtDir=${runDir}/genotyping
mkdir -p ${gtDir} || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFFilter \
    -R ${referenceFile} \
    -vcf ${runDir}/${filePrefix}.sites.vcf.gz \
    -includeSite ${filtDir}/${filePrefix}.sites.list \
    -O ${gtDir}/${filePrefix}.sites.vcf.gz \
    || exit 1

numSites=$(cat ${filtDir}/${filePrefix}.sites.list | wc -l)
echo "Genotyping ${numSites} selected variants..."

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVGenotyper.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R ${referenceFile} \
    -md ${mdPath} \
    -I ${bamFileList} \
    -runDirectory ${gtDir} \
    -jobLogDir ${runDir}/logs \
    -jobRunner ParallelShell \
    -vcf ${gtDir}/${filePrefix}.sites.vcf.gz \
    -O ${gtDir}/${filePrefix}.genotypes.vcf.gz \
    --disableJobReport \
    -run \
    -retry 3 \
    || exit 1


