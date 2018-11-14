#!/bin/bash

archiveList=$1
binSize=$2
referenceBundle=$3

if [ -z "${referenceBundle}" ]; then
    echo "Usage: gs_build_profiles.sh <archiveList> <binSize> <referenceBundle>" && exit 1
fi

source $(dirname $0)/gs_extract_reference.sh ${referenceBundle} || exit 1

mkdir -p tmp
mdDir=metadata
profilesDir=profiles_${binSize}
mkdir -p ${mdDir}
mkdir -p ${profilesDir}/logs

if [ -z "${archiveList}" ]; then
    echo "Archive list is empty" && exit 1
fi

# Support tar.gz format of per-sample archives
if [[ $(head -1 ${archiveList}) =~ tar.gz$ ]]; then
    echo $(date +"[%b %d %H:%M:%S] Unarchiving metadata")
    java -cp ${SV_CLASSPATH} -Xmx4g \
        org.broadinstitute.gatk.queue.QCommandLine \
        -S ${SV_DIR}/scripts/firecloud/UnarchiveData.q \
        -S ${SV_DIR}/qscript/SVQScript.q \
        -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
        -cp ${SV_CLASSPATH} \
        -R ${referenceFile} \
        -archive ${archiveList} \
        -deleteArchives true \
        -O ${mdDir} \
        -tempDir tmp \
        -jobLogDir tmp \
        -jobRunner ParallelShell \
        -run \
        || exit 1

    mdDirArg="-md "$(find ${mdDir} -name headers.bam -exec dirname {} \; | xargs | sed 's/ / -md /g')

else
    idx=0
    while read mdPath; do
        idx=$((idx+1))
        mkdir -p ${mdDir}/${idx}

        cd ${mdDir}/${idx}
        java -Xmx2g -cp ${SV_CLASSPATH} \
            org.broadinstitute.sv.apps.ZipExtract \
            -I ${mdPath} \
            -file headers.bam \
            -file headers.bam.bai \
            -file sample_gender.report.txt \
            || exit 1
        cd -
    done < ${archiveList}

    mdDirArg=$(awk '{print "-md " $0}' ${archiveList} | xargs)

fi

echo "MD_DIR=${mdDirArg}"
find ${mdDir} -name headers.bam > bam_files.list
cat bam_files.list

# Store the merged genderMapFile within the profiles directory, as this gender map is required to run the LCNV pipeline
genderMapFile=${profilesDir}/sample_gender.report.txt
find ${mdDir} -name sample_gender.report.txt | head -1 | xargs head -1 > ${genderMapFile}
find ${mdDir} -name sample_gender.report.txt -exec tail -n +2 {} \; >> ${genderMapFile}

echo $(date +"[%b %d %H:%M:%S] Building profiles")
java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/profiles/GenerateDepthProfiles.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -cp ${SV_CLASSPATH} \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir tmp \
    -jobLogDir ${profilesDir}/logs \
    -R ${referenceFile} \
    ${mdDirArg} \
    -I bam_files.list \
    -profileBinSize ${binSize} \
    -runDirectory ${profilesDir} \
    -jobRunner ParallelShell \
    -run \
    || exit 1

echo $(date +"[%b %d %H:%M:%S] Archiving profiles directory")
tar -cvzf profiles_${binSize}.tar.gz ${profilesDir} || exit 1

echo $(date +"[%b %d %H:%M:%S] Building profiles completed successfully")

