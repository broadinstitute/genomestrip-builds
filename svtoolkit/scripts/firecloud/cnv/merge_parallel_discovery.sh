#!/bin/bash

mdPath=$1
partitionsArchiveList=$2
referenceBundle=$3
credentialsKeyFile=$4

if [ -z "${credentialsKeyFile}" ]; then
    echo "Usage: merge_parallel_discovery.sh <mdPath> <partitionsArchiveList> <referenceBundle> <credentialsKeyFile>" && exit 1
fi

echo "Exporting GOOGLE_APPLICATION_CREDENTIALS..."
export GOOGLE_APPLICATION_CREDENTIALS=${credentialsKeyFile}

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}
runDir=cnv_output
mkdir -p ${runDir}/eval || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/scripts/firecloud/UnarchiveData.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -R ${referenceFile} \
    -archive ${partitionsArchiveList} \
    -deleteArchives true \
    -O cnv_partitions \
    -tempDir ${SV_TMPDIR} \
    -jobLogDir ${SV_TMPDIR} \
    -jobRunner ParallelShell \
    -run \
    || exit 1

find cnv_partitions/*/cnv_stage11/ -name '*.filtered.genotypes.vcf.gz' > vcf_files.list || exit 1

java -cp $SV_CLASSPATH -Xmx4g \
    org.broadinstitute.sv.apps.MergeVCFFiles \
    -R ${referenceFile} \
    -vcf vcf_files.list \
    -O cnv_partitions/gs_cnv.genotypes.vcf.gz \
    || exit 1

java -Xmx2g -cp ${SV_CLASSPATH} \
    org.broadinstitute.sv.apps.ZipExtract \
    -I ${mdPath} \
    -file sample_gender.report.txt \
    || exit 1

duplicateScoreThreshold=0.0
duplicateOverlapThreshold=0.5

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -A Redundancy \
    -R ${referenceFile} \
    -md ${mdPath} \
    -vcf cnv_partitions/gs_cnv.genotypes.vcf.gz \
    -comparisonFile cnv_partitions/gs_cnv.genotypes.vcf.gz \
    -filterVariants false \
    -duplicateOverlapThreshold ${duplicateOverlapThreshold} \
    -duplicateScoreThreshold ${duplicateScoreThreshold} \
    -O ${runDir}/gs_cnv.annotated.genotypes.vcf.gz \
    || exit 1

vcfextract -fields ID,INFO:GSDUPLICATEOVERLAP ${runDir}/gs_cnv.annotated.genotypes.vcf.gz \
| awk '$2 == "NA" {print $1}' > ${runDir}/SelectedSites.list

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFFilter \
    -R ${referenceFile} \
    -vcf ${runDir}/gs_cnv.annotated.genotypes.vcf.gz \
    -includeSite ${runDir}/SelectedSites.list \
    -O ${runDir}/gs_cnv.genotypes.vcf.gz \
    || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -A GCContent \
    -A CopyNumberClass \
    -R ${referenceFile} \
    -md ${mdPath} \
    -genderMapFile sample_gender.report.txt \
    -ploidyMapFile ${ploidyMapFile} \
    -vcf ${runDir}/gs_cnv.genotypes.vcf.gz \
    -filterVariants false \
    -writeReport true \
    -writeSummary true \
    -reportDirectory ${runDir}"/eval" \
    || exit 1

tar -cvzf cnv_output.tar.gz ${runDir} || exit 1

