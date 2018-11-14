#!/bin/bash

windowSize=$1
windowOverlap=$2
intervals="$3"
partitionsFile=$4
referenceBundle=$5

if [ -z "${referenceBundle}" ]; then
    echo "Usage: compute_partitions.sh <windowSize> <windowOverlap> <intervals> <partitionsFile> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1

genomeMaskArg="-genomeMaskFile ${svMaskFile} -genomeMaskFile ${lcMaskFile}"

intervalsArg=
if [ "${intervals}" != "" ]; then
    intervalsArg=$(echo "${intervals}" | sed 's/ /\n/g' | awk '{print "-L " $0}' | xargs)
fi
echo "intervalsArg=${intervalsArg}"

java -Xmx2g -cp $SV_CLASSPATH \
    org.broadinstitute.sv.discovery.SVDepthScanner \
    -R $referenceFile \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    ${parameterFileArg} \
    ${genomeMaskArg} \
    -tilingWindowSize ${windowSize} \
    -tilingWindowOverlap ${windowOverlap} \
    -maximumReferenceGapLength 1000000 \
    ${intervalsArg} \
    -O large.windows.vcf.gz \
    || exit 1

java -Xmx2g -cp $SV_CLASSPATH \
    org.broadinstitute.sv.discovery.SVDepthScanner \
    -R $referenceFile \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    ${parameterFileArg} \
    ${genomeMaskArg} \
    -tilingWindowSize 1000 \
    -tilingWindowOverlap 500 \
    -maximumReferenceGapLength 1000 \
    ${intervalsArg} \
    -O small.windows.vcf.gz \
    || exit 1

vcfextract -fields CHROM,POS,INFO:END small.windows.vcf.gz > small.windows.dat || exit 1

# When an interval chr:pos-end with pos>1 is passed to SVDepthScanner, it will start the first window at pos-1
vcfextract -fields CHROM,POS,INFO:END large.windows.vcf.gz \
| tail -n +2 \
| while read chrom pos end; do
    aligned=$(
        awk -v chrom=${chrom} -v pos=${pos} -v end=${end}  '{
            if ($1 != chrom) next
            if ($2 >= pos && !printed) {
                printf($2 + 1)
                printed = 1
            }
            if ($3 >= end) {
                printf("-%s", $3)
                exit
            }
        }' small.windows.dat
    )
    echo "${chrom}:${aligned}"
done \
| awk '{
    printf("P%04d\t%s\n", NR, $0)
}' > ${partitionsFile}
