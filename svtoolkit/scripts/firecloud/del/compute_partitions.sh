#!/bin/bash

searchIntervals="$1"
partitionsFile=$2
referenceBundle=$3

if [ -z "${referenceBundle}" ]; then
    echo "Usage: compute_partitions.sh <searchIntervals> <partitionsFile> <referenceBundle>" && exit 1
fi

if [[ $referenceBundle =~ .fasta$ ]]; then
    referenceFile=${referenceBundle}
else
    source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
fi

searchIntervalsArg=
if [ "${searchIntervals}" != "" ]; then
    searchIntervalsArg=$(echo "${searchIntervals}" | sed 's/ /\n/g' | awk '{print "-searchInterval " $0}' | xargs)
fi
echo "searchIntervalsArg=${searchIntervalsArg}"

java -Xmx2g \
    -cp $SV_CLASSPATH \
    org.broadinstitute.sv.queue.ComputeDiscoveryPartitions \
    -searchWindowSize 5000000 \
    -searchWindowPadding 10000 \
    -searchMinimumSize 100 \
    -searchMaximumSize 100000 \
    ${searchIntervalsArg} \
    -R $referenceFile \
    -O ${partitionsFile} \
    || exit 1

