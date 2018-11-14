#!/bin/bash

partitionName=$1
mdPath=$2
intervalListFile=$3
referenceBundle=$4
credentialsKeyFile=$5

if [ -z "${credentialsKeyFile}" ]; then
    echo "Usage: run_parallel_discovery.sh <partitionName> <mdPath> <intervalListFile> <referenceBundle> <credentialsKeyFile>" && exit 1
fi

echo "Exporting GOOGLE_APPLICATION_CREDENTIALS..."
export GOOGLE_APPLICATION_CREDENTIALS=${credentialsKeyFile}

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}

genomeMaskArg="-genomeMaskFile ${svMaskFile} -genomeMaskFile ${lcMaskFile}"

bamFileList=headers.bam

java -Xmx2g -cp ${SV_CLASSPATH} \
    org.broadinstitute.sv.apps.ZipExtract \
    -I ${mdPath} \
    -file headers.bam \
    -file headers.bam.bai \
    -file sample_gender.report.txt \
    || exit 1

runDir=cnv_output
mkdir -p ${runDir}
for stage in {1..11}; do
    mkdir ${runDir}/cnv_stage${stage}
done

# Stage1
echo "Running Stage1"

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.discovery.SVDepthScanner \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -R ${referenceFile} \
    ${genomeMaskArg} \
    -L ${intervalListFile} \
    -tilingWindowSize 1000 \
    -tilingWindowOverlap 500 \
    -maximumReferenceGapLength 1000 \
    -O ${runDir}/cnv_stage1/cnv.sites.vcf.gz \
    || exit 1

# Stage2
echo "Running Stage2"

source ${SV_DIR}/scripts/firecloud/cnv/run_cnv_genotyper.sh \
    ${runDir}/cnv_stage1/cnv.sites.vcf.gz \
    ${mdPath} \
    ${bamFileList} \
    ${runDir}/cnv_stage2 \
    ${referenceFile} \
    ${credentialsKeyFile} \
    || exit 1

# Stage3
echo "Running Stage3"

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.discovery.MergeDepthScannerSites \
    -R ${referenceFile} \
    ${genomeMaskArg} \
    -vcf ${runDir}/cnv_stage2/cnv.genotypes.vcf.gz \
    -reportFile ${runDir}/cnv_stage3/MergeScannerSites.report.dat \
    -outputVariants ALL \
    -duplicateScoreThresholdMax 0.0 \
    -O ${runDir}/cnv_stage3/merged.sites.vcf.gz \
    || exit 1

# Stage4
echo "Running Stage4"

source ${SV_DIR}/scripts/firecloud/cnv/run_cnv_genotyper.sh \
    ${runDir}/cnv_stage3/merged.sites.vcf.gz \
    ${mdPath} \
    ${bamFileList} \
    ${runDir}/cnv_stage4 \
    ${referenceFile} \
    ${credentialsKeyFile} \
    || exit 1

# Stage6
echo "Running Stage6"

awk '$2 == "PASS"' ${runDir}/cnv_stage4/eval/GenotypeSiteFilters.report.dat \
| cut -f1 \
| joinfiles -select \
    ${runDir}/cnv_stage4/eval/CopyNumberClass.report.dat 1 - 1 \
| awk '$8 != "NA"' | cut -f 1 \
> ${runDir}/cnv_stage6/SelectedVariants.list

# Stage7
echo "Running Stage7"

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.genotyping.RefineCNVBoundaries \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    ${parameterFileArg} \
    -R ${referenceFile} \
    -ploidyMapFile ${ploidyMapFile} \
    -md ${mdPath} \
    ${genomeMaskArg} \
    -genderMapFile sample_gender.report.txt \
    -I ${bamFileList} \
    -vcf ${runDir}/cnv_stage4/merged.genotypes.vcf.gz \
    -site cnv_output/cnv_stage6/SelectedVariants.list \
    -boundaryPrecision 100 \
    -minimumRefinedLength 500 \
    -maximumReferenceGapLength 1000 \
    -O ${runDir}/cnv_stage7/merged.brig.vcf \
    || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.discovery.MergeBrigVcfFiles \
    -R ${referenceFile} \
    -vcfFile ${runDir}/cnv_stage7/merged.brig.vcf \
    -mergedVcfFile ${runDir}/cnv_stage7/brig.sites.vcf.gz \
    || exit 1

# Stage8
echo "Running Stage8"

source ${SV_DIR}/scripts/firecloud/cnv/run_cnv_genotyper.sh \
    ${runDir}/cnv_stage7/brig.sites.vcf.gz \
    ${mdPath} \
    ${bamFileList} \
    ${runDir}/cnv_stage8 \
    ${referenceFile} \
    ${credentialsKeyFile} \
    || exit 1

# Stage9
echo "Running Stage9"
mkdir ${runDir}/cnv_stage9/eval

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.discovery.MergeDepthScannerSites \
    -R ${referenceFile} \
    ${genomeMaskArg} \
    -vcf ${runDir}/cnv_stage8/brig.genotypes.vcf.gz \
    -reportFile ${runDir}/cnv_stage9/eval/MergeAdjacentSites.report.dat \
    -outputVariants ALL \
    -adjacentSiteDistance 1000000 \
    -duplicateScoreThresholdMax 0 \
    -O ${runDir}/cnv_stage9/adjacent_merged.sites.vcf.gz \
    || exit 1

# Stage10
echo "Running Stage10"

source ${SV_DIR}/scripts/firecloud/cnv/run_cnv_genotyper.sh \
    ${runDir}/cnv_stage9/adjacent_merged.sites.vcf.gz \
    ${mdPath} \
    ${bamFileList} \
    ${runDir}/cnv_stage10 \
    ${referenceFile} \
    ${credentialsKeyFile} \
    || exit 1

# Stage11
echo "Running Stage11"

java -Xmx4g -jar $GATK_JAR \
    -T VariantFiltration \
    -V ${runDir}/cnv_stage10/adjacent_merged.genotypes.vcf.gz \
    -o ${runDir}/cnv_stage11/filtered.genotypes.vcf.gz \
    -R $referenceFile \
    -filterName LENGTH -filter "GSCNCATEGORY == \"NA\" || (GSCNCATEGORY == \"DEL\" || GSCNCATEGORY == \"MIXED\") && GCLENGTH < 1000 || GSCNCATEGORY == \"DUP\" && GCLENGTH < 2000" \
    -filterName CALLRATE -filter "GSCALLRATE == \"NA\" || GSCALLRATE < 0.9" \
    -filterName DENSITY -filter "(1.0 * GSELENGTH / GCLENGTH) < 0.5" \
    -filterName CLUSTERSEP -filter "GSCLUSTERSEP == \"NA\" || GSCLUSTERSEP < 5.0" \
    -filterName VDJREGION -filter "GSVDJFRACTION > 0.0" \
    || exit 1

tar -cvzf cnv_output.tar.gz ${runDir}
#tar -C ${runDir} -cvzf cnv_output.tar.gz cnv_stage11 || exit 1
