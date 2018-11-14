#!/bin/bash

filePrefix=$1
searchInterval=$2
mdPath=$3
bamFileList=$4
parameterArgs="$5"
referenceBundle=$6
repeatTrackFile=$7

if [ -z "${repeatTrackFile}" ]; then
    echo "Usage: run_complete_discovery.sh <filePrefix> <searchInterval> <mdPath> <bamFileList> <parameterArgs> <referenceBundle> <repeatTrackFile>" && exit 1
fi

if [[ $referenceBundle =~ .fasta$ ]]; then
    referenceFile=${referenceBundle}
    source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1
else
    source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
fi

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}
runDir=del_output
mkdir -p ${runDir}/logs || exit 1

storeReadPairFile="true"

$(dirname $0)/compute_partitions.sh ${searchInterval} partition.dat ${referenceFile} || exit 1

numPartitions=$(cat partition.dat | wc -l)
if [ "${numPartitions}" -ne 1 ]; then
    echo "Error: this method can be run on one partition only, numPartitions=${numPartitions}" && exit 1
fi

# Avoid problem with certain read pairs - this should be debugged
extraArgs="-P select.validateReadPairs:false"

for attempt in {1..3}; do

    echo $(date +"[%b %d %H:%M:%S] Running Deletion discovery, attempt=${attempt}")

# TODO
#    -genderMapFile sample_gender.report.txt \

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVDiscovery \
    -T SVDiscoveryWalker \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R ${referenceFile} \
    -md ${mdPath} \
    -genomeMaskFile ${svMaskFile} \
    -I ${bamFileList} \
    -disableGATKTraversal true \
    $(cut -f 2- partition.dat) \
    -runFilePrefix ${filePrefix} \
    -storeReadPairFile ${storeReadPairFile} \
    ${extraArgs} \
    ${parameterArgs} \
    -runDirectory ${runDir} \
    -O ${runDir}/${filePrefix}.unfiltered.sites.vcf.gz

    rc=$?
    echo $(date +"[%b %d %H:%M:%S] Deletion discovery completed, rc=${rc}")

    if [ "${rc}" -eq 0 ]; then
        break
    fi
done

if [ "${rc}" -ne 0 ]; then
    echo "Error running Deletion discovery" && exit 1
fi

samtools index ${runDir}/${filePrefix}.discovery.pairs.bam || exit 1

$(dirname $0)/apply_default_filters.sh \
    ${runDir}/${filePrefix}.unfiltered.sites.vcf.gz \
    ${runDir}/${filePrefix}.sites.vcf.gz \
    ${referenceFile} \
    || exit 1

$(dirname $0)/filter_discovered_sites.sh \
    ${filePrefix} \
    ${runDir} \
    ${referenceFile} \
    ${repeatTrackFile} \
    || exit 1

gtDir=${runDir}/genotyping
mkdir -p ${gtDir}/eval || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.apps.VCFFilter \
    -R ${referenceFile} \
    -vcf ${runDir}/${filePrefix}.sites.vcf.gz \
    -includeSite ${runDir}/filtering/${filePrefix}.sites.list \
    -O ${gtDir}/${filePrefix}.sites.vcf.gz \
    || exit 1

# TODO
#    -genderMapFile sample_gender.report.txt \

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVGenotyper \
    -T SVGenotyperWalker \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R ${referenceFile} \
    -md ${mdPath} \
    -genomeMaskFile ${svMaskFile} \
    -I ${runDir}/${filePrefix}.discovery.pairs.bam \
    -disableGATKTraversal true \
    -vcf ${gtDir}/${filePrefix}.sites.vcf.gz \
    ${parameterArgs} \
    -runDirectory ${gtDir} \
    -O ${runDir}/${filePrefix}.genotypes.vcf.gz \
    || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -R ${referenceFile} \
    -vcf ${runDir}/${filePrefix}.genotypes.vcf.gz \
    -auxFilePrefix ${gtDir}/${filePrefix}.genotypes \
    -A GCContent \
    -A ClusterSeparation \
    -writeReport true \
    -writeSummary true \
    -reportDirectory ${gtDir}/eval \
    || exit 1

gcContentReport=${gtDir}/eval/GCContent.report.dat
clusterSepReport=${gtDir}/eval/ClusterSeparation.report.dat

cat ${clusterSepReport} \
    | awk 'NR > 1 && $4 != "NA" && $4 > 3.0 { print $1 }' \
    > ${gtDir}/SelectedCalls.list || exit 1

java -cp ${SV_CLASSPATH} -Xmx3g \
    org.broadinstitute.sv.apps.VCFFilter \
    -R ${referenceFile} \
    -vcf ${gtDir}/${filePrefix}.sites.vcf.gz \
    -includeSite ${gtDir}/SelectedCalls.list \
    -includeInfoTag END \
    -includeInfoTag CIPOS \
    -includeInfoTag CIEND \
    -includeInfoTag GSBKPT \
    -includeInfoTag GSNPAIRS \
    -includeInfoTag GSNSAMPLES \
    -O ${filePrefix}.sites.vcf.gz \
    || exit 1

