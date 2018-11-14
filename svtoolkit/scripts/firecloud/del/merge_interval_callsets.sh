#!/bin/bash

batchList=$1
callsetList=$2
homologyFileList=$3
discoveryPairsFileList=$4
discoveryPairsIndexFileList=$5
genotypesPairsFileList=$6
genotypesPairsIndexFileList=$7
referenceBundle=$8

if [ "${referenceBundle}" == "" ]; then
    echo "Usage: " && exit 1
fi

source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

mkdir merged
mkdir final 

java -cp $SV_CLASSPATH -Xmx4g \
    org.broadinstitute.sv.apps.VCFMerge \
    -R ${referenceFile} \
    -vcf ${callsetList} \
    -includeInfoTag END \
    -includeInfoTag GSELENGTH \
    -includeInfoTag SVTYPE \
    -includeInfoTag CIPOS \
    -includeInfoTag CIEND \
    -includeInfoTag GSBKPT \
    -includeInfoTag GSNPAIRS \
    -includeInfoTag GSNSAMPLES \
    -O merged/gs_del.genotypes.vcf.gz \
    || exit 1

java -cp $SV_CLASSPATH -Xmx4g \
    org.broadinstitute.sv.apps.FilterRedundantSites \
    -R ${referenceFile} \
    -vcf merged/gs_del.genotypes.vcf.gz \
    -O final/gs_del.genotypes.vcf.gz \
    || exit 1

${SV_DIR}/scripts/firecloud/del/merge_homology_files.sh ${homologyFileList} final/gs_dels.homology.gz || exit 1

numbatches=$(cat ${batchList} | wc -l)
for idx in `seq 1 $numbatches`; do
    batch=$(head -${idx} ${batchlist} | tail -1)
    mkdir -p final/${batch}

    mv $(head -${idx} ${discoveryPairsFileList} | tail -1) final/${batch}
    mv $(head -${idx} ${discoveryPairsIndexFileList} | tail -1) final/${batch}
    mv $(head -${idx} ${genotypesPairsFileList} | tail -1) final/${batch}
    mv $(head -${idx} ${genotypesPairsIndexFileList} | tail -1) final/${batch}
done

tar -cvzf del_output.tar.gz final/

