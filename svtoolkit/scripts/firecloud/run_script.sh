#!/bin/bash

scriptPath=$1
scriptArgs="$2"
resourceArchive=$3
gsOutputPath=$4
outputArchive=$5
referenceBundle=$6
credentialsKeyFile=$7

if [ "${resourceArchive}" == "" ]; then
    echo "Usage: run_script.sh <scriptPath> <scriptArgs> <resourceArchive> <gsOutputPath> <outputArchive> <referenceBundle> <credentialsKeyFile>" && exit 1
fi

if [ ! -z "${credentialsKeyFile}" ]; then
    source ${SV_DIR}/scripts/firecloud/gs_activate_service_account.sh ${credentialsKeyFile}
fi

if [ ! -z "${referenceBundle}" ]; then
    source $(dirname $0)/gs_extract_reference.sh ${referenceBundle} || exit 1
fi

tar -xvzf ${resourceArchive} || exit 1

source ${scriptPath} ${scriptArgs} || exit 1
echo "Script ${scriptPath} completed successfully"

if [ -e ${outputArchive} ]; then
    echo "Copying ${outputArchive} to ${gsOutputPath}..."
    gsutil cp ${outputArchive} ${gsOutputPath}/${outputArchive} || exit 1
else
    echo "${outputArchive} does not exist"
fi
