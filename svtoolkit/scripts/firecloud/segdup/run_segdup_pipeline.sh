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
mkdir -p ${runDir}/logs
mkdir -p ${runDir}/vcfs
SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}

java -Xmx2g -cp ${SV_CLASSPATH} \
    org.broadinstitute.sv.apps.ZipExtract \
    -I ${mdPath} \
    -file headers.bam \
    -file headers.bam.bai \
    -file sample_gender.report.txt \
    || exit 1

# Stage1
java -cp ${SV_CLASSPATH} -Xmx2g \
    org.broadinstitute.sv.discovery.GenerateSegdupSites \
    -R ${referenceFile} \
    -segdupFile ${segdupFile} \
    -sequenceMapFile ${sequenceMapFile} \
    -O ${runDir}/vcfs/segdup_scan.sites.vcf \
    || exit 1

# Stage2
java -cp ${SV_CLASSPATH} -Xmx2g \
    -Djava.io.tmpdir=${SV_TMPDIR} \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVGenotyper2.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -R ${referenceFile} \
    -md ${mdPath} \
    -genderMapFile sample_gender.report.txt \
    -I headers.bam \
    -vcf ${runDir}/vcfs/segdup_scan.sites.vcf \
    -parallelRecords 1000 \
    -P genotyping.modules:depth \
    -O ${runDir}/segdup_scan.genotypes.vcf \
    -runDirectory ${runDir} \
    -jobLogDir ${runDir}/logs \
    -jobRunner ParallelShell \
    -availableProcessorsMultiplier 0.625 \
    -disableGATKTraversal \
    -run

rc=$?
if [ "${rc}" -eq 0 ]; then
    find ${runDir} -maxdepth 1 -type f -name 'P[0-9]*.genotypes.*' -exec rm {} \;
    find ${runDir} -maxdepth 1 -type f -name '[.]*P[0-9]*.genotypes.*' -exec rm {} \;
fi

tar -cvzf segdup_output.tar.gz ${runDir}

exit ${rc}
