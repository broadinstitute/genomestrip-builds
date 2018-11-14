#!/bin/bash

credentialsKeyFile=$1
if [ "${credentialsKeyFile}" == "" ]; then
    echo "Usage: gs_activate_service_account.sh <credentialsKeyFile>" && exit 1
fi

echo "Activating google service account ${credentialsKeyFile}"
gcloud auth activate-service-account --key-file=${credentialsKeyFile}

echo "Exporting GOOGLE_APPLICATION_CREDENTIALS..."
export GOOGLE_APPLICATION_CREDENTIALS=${credentialsKeyFile}

