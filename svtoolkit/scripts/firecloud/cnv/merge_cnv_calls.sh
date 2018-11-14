#!/bin/bash

cnvCallsList=$1
referenceBundle=$2

if [ "${referenceBundle}" == "" ]; then
    echo "Usage: <cnvCallsList> <referenceBundle>" && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

mergedSitesVcf="gs_cnv.sites.vcf.gz"

xargs cat < ${cnvCallsList} \
| sort -u \
> SelectedSites.dat
echo "There are "$(cat SelectedSites.dat | wc -l)" unique sites"

(cat ${SV_DIR}/scripts/firecloud/templates/vcf_template.vcf;
joinfiles \
    SelectedSites.dat 2 \
    <(awk '{print $1 "\t" NR}' ${referenceFile}.fai) 1 \
| sort -k6n -k3n -k4n \
| awk -v OFS="\t" '{
    print $2, $3, $1, "N", "<CNV>", ".", ".", "END=" $4 ";SVTYPE=CNV"
}') \
| bgzip -c > ${mergedSitesVcf} || exit 1

echo "Completed successfully"

