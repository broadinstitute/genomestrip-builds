#/bin/bash

vcfFile=$1

echo "Fixing $vcfFile"

zcat ${vcfFile} \
| sed '/^#/! s/ /_/g' \
| bgzip -c > tmp.vcf.gz

mv tmp.vcf.gz ${vcfFile}
