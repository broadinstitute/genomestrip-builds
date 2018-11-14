#!/bin/bash

delOutput=$1
referenceBundle=$2

if [ "${referenceBundle}" == "" ]; then
    echo "Usage: <delOutput> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

tar -xvzf ${delOutput} del_output/genotyping --strip-components 1 || exit 1
rm ${delOutput}

vcfFile=genotyping/gs_dels.genotypes.vcf.gz
gcContentReport=genotyping/eval/GCContent.report.dat
clusterSepReport=genotyping/eval/ClusterSeparation.report.dat

auxFilePrefix=genotyping/gs_dels.genotypes

mkdir genotyping/eval

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVAnnotator \
    -R ${referenceFile} \
    -vcf ${vcfFile} \
    -auxFilePrefix ${auxFilePrefix} \
    -A GCContent \
    -A ClusterSeparation \
    -writeReport true \
    -writeSummary true \
    -reportDirectory genotyping/eval \
    || exit 1

vcfextract -fields ID,REF ${vcfFile} \
| tail -n +2 \
| joinfiles \
    - 1 \
    ${clusterSepReport} 1 \
| awk '$6 != "NA" && $6 >= 3' \
| cut -f 1-2 \
| joinfiles \
    - 1 \
    ${gcContentReport} 1 \
| awk -v OFS="\t" '{
    print $4, $5, $6, $2
}' \
> SelectedCalls.dat

