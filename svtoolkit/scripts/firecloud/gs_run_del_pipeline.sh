#!/bin/bash

mdPath=$1
bamFileList=$2
intervalList="$3"
referenceBundle=$4
repeatTrackFile=$5
credentialsKeyFile=$6

if [ -z "${credentialsKeyFile}" ]; then
    echo "Usage: gs_run_del_pipeline.sh <mdPath> <bamFileList> <intervalList> <referenceBundle> <repeatTrackFile> <credentialsKeyFile>" && exit 1
fi

echo "SV_CLASSPATH=${SV_CLASSPATH}"

echo "Exporting GOOGLE_APPLICATION_CREDENTIALS..."
export GOOGLE_APPLICATION_CREDENTIALS=${credentialsKeyFile}

source $(dirname $0)/gs_extract_reference.sh ${referenceBundle} || exit 1

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}
runDir=del_output
mkdir -p ${runDir}/logs || exit 1

filePrefix=gs_dels

# Avoid problem with certain read pairs - this should be debugged
extraArgs="-P select.validateReadPairs:false"

intervalListArg=
if [ "${intervalList}" != "" ]; then
    intervalListArg="-L "$(echo ${intervalList} | sed 's/ /-L /g')
fi

echo "Running Deletion discovery pipeline..."

for attempt in {1..2}; do

    echo $(date +"[%b %d %H:%M:%S] Running Deletion pipeline, attempt=${attempt}")

java -cp ${SV_CLASSPATH} -Xmx2g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVDiscovery.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    --disableJobReport \
    -cp ${SV_CLASSPATH} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R ${referenceFile} \
    -md ${mdPath} \
    -runDirectory ${runDir} \
    -jobLogDir ${runDir}/logs \
    -jobRunner ParallelShell \
    -windowSize 5000000 \
    -windowPadding 10000 \
    -minimumSize 100 \
    -maximumSize 100000 \
    ${intervalListArg} \
    -I ${bamFileList} \
    -O ${runDir}/${filePrefix}.sites.vcf \
    ${extraArgs} \
    -run \
    -retry 3

    rc=$?
    echo $(date +"[%b %d %H:%M:%S] Deletion pipeline completed, rc=${rc}")

    if [ "${rc}" -eq 0 ]; then
        break
    fi
done

# For now, produce intermediate output FireCloud can save in Google Storage
tar -cvzf del_output.tar.gz ${runDir} || exit 1

if [ "${rc}" -ne 0 ]; then
    echo "Error running Deletion pipeline" && exit 1
fi

# Filter calls in alpha-satellite regions
echo "Filtering the discovered variants..."

filtDir=${runDir}/filtering
mkdir -p ${filtDir} || exit 1

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -R ${referenceFile} \
    -A MobileElements \
    -repeatTrackFile ${repeatTrackFile} \
    -vcf ${runDir}/${filePrefix}.sites.vcf \
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
    -vcf ${runDir}/${filePrefix}.sites.vcf \
    -includeSite ${filtDir}/${filePrefix}.sites.list \
    -O ${gtDir}/${filePrefix}.sites.vcf.gz \
    || exit 1

echo "Genotyping the selected variants..."

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
    -parallelRecords 10 \
    -vcf ${gtDir}/${filePrefix}.sites.vcf.gz \
    -O ${gtDir}/${filePrefix}.genotypes.vcf.gz \
    --disableJobReport \
    -run \
    -retry 3 \
    || exit 1

evalDir=${gtDir}/eval
mkdir -p ${evalDir} || exit 1

# Generate useful report files
echo "Generating reports..."

java -Xmx1g -cp ${SV_CLASSPATH} org.broadinstitute.sv.apps.VCFExtract \
    -fields ID,CHROM,POS,INFO:END,INFO:GSELENGTH \
    ${gtDir}/${filePrefix}.genotypes.vcf.gz \
    | awk -v OFS="\t" '{ if (NR == 1) { $6 = "LENGTH"; $7 = "DENSITY" } else { $6 = $4-$3+1; $7 = sprintf("%1.4f", $5/$6) }; print }' \
    > ${evalDir}/VariantLength.report.dat || exit 1

auxFilePrefix="${gtDir}/${filePrefix}.genotypes"

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
        -R ${referenceFile} \
        -vcf ${gtDir}/${filePrefix}.genotypes.vcf.gz \
        -auxFilePrefix ${auxFilePrefix} \
        -A CopyNumberClass \
        -A ClusterSeparation \
        -A VariantsPerSample \
        -A GCContent \
        -writeReport true \
        -writeSummary true \
        -reportDirectory ${evalDir} \
        || exit 1

cat ${evalDir}/VariantsPerSample.report.dat | cut -f 1 | tail -n +2 > ${runDir}/samples.list || exit 1

echo "Deletion discovery pipeline completed successfully"

tar -cvzf del_output.tar.gz ${runDir} || exit 1


