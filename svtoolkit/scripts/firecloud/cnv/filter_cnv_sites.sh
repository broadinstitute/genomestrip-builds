#!/bin/bash

cnvOutput=$1
mdPath=$2
madRange=$3
referenceBundle=$4
credentialsKeyFile=$5

if [ -z "${credentialsKeyFile}" ]; then
    echo "Usage: filter_cnv_sites.sh <cnvOutput> <mdPath> <madRange> <referenceBundle> <credentialsKeyFile>" && exit 1
fi

echo "Exporting GOOGLE_APPLICATION_CREDENTIALS..."
export GOOGLE_APPLICATION_CREDENTIALS=${credentialsKeyFile}

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1

tar -xvzf ${cnvOutput} || exit 1
rm ${cnvOutput}

java -Xmx2g -cp ${SV_CLASSPATH} \
    org.broadinstitute.sv.apps.ZipExtract \
    -I ${mdPath} \
    -file sample_gender.report.txt \
    || exit 1

awk -F '\t' 'NR > 1 && $2 != "chrX" && $2 != "chrY" {print $1}' cnv_output/eval/GCContent.report.dat > autosomal_sites.list
echo "There are "$(cat autosomal_sites.list | wc -l)" autosomal sites"

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFFilter \
    -R ${referenceFile} \
    -vcf cnv_output/gs_cnv.genotypes.vcf.gz \
    -includeSite autosomal_sites.list \
    -O autosomal.genotypes.vcf.gz \
    || exit 1

mkdir eval
java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -A VariantsPerSample \
    -R ${referenceFile} \
    -md ${mdPath} \
    -genderMapFile sample_gender.report.txt \
    -ploidyMapFile ${ploidyMapFile} \
    -vcf autosomal.genotypes.vcf.gz \
    -filterVariants false \
    -writeReport true \
    -writeSummary true \
    -reportDirectory eval \
    || exit 1

Rscript ${SV_DIR}/scripts/firecloud/cnv/compute_discovery_samples.R "eval/VariantsPerSample.report.dat" ${madRange} "eval/DiscoverySamples.list" || exit 1
echo "There are "$(cat eval/DiscoverySamples.list | wc -l)" discovery samples"

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -A GCContent \
    -A CopyNumberClass \
    -R ${referenceFile} \
    -md ${mdPath} \
    -genderMapFile sample_gender.report.txt \
    -ploidyMapFile ${ploidyMapFile} \
    -vcf cnv_output/gs_cnv.genotypes.vcf.gz \
    -sample eval/DiscoverySamples.list \
    -filterVariants false \
    -writeReport true \
    -writeSummary true \
    -reportDirectory eval \
    || exit 1

joinfiles \
    eval/CopyNumberClass.report.dat 1 \
    eval/GCContent.report.dat 1 \
| awk '($8 == "DEL" && $15 >= 1000) || ($8 == "DUP" || $8 == "MIXED") && $15 >= 2000' | cut -f 10-13 \
> eval/SelectedSites.dat

