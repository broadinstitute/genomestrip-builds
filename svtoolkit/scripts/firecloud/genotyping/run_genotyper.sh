#!/bin/bash

vcfFile=$1
mdPath=$2
bamFileList=$3
parameterArgs="$4"
referenceBundle=$5
useLcMask=$6

if [ -z "${useLcMask}" ]; then
    echo "Usage: run_genotyper.sh <vcfFile> <mdPath> <bamFileList> <parameterArgs> <referenceBundle> <useLcMask>" && exit 1
fi

if [[ $referenceBundle =~ .fasta$ ]]; then
    referenceFile=${referenceBundle}
    source ${SV_DIR}/scripts/firecloud/set_sv_dataset.sh ${referenceFile} || exit 1
else
    source ${SV_DIR}/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1
fi

SV_TMPDIR=tmp
mkdir -p ${SV_TMPDIR}
runDir=gtrun
mkdir -p ${runDir}/logs || exit 1

outputVcf="${runDir}/"$(echo ${vcfFile} | awk -F / '{ print $NF}' | sed 's/.sites.vcf.gz$/.vcf.gz/' | sed 's/.genotypes.vcf.gz$/.vcf.gz/' | sed 's/.vcf.gz$/.genotypes.vcf.gz/')

extraArgs=
bamFileArg=
if [ "${bamFileList}" == "NULL" ]; then
    extraArgs="-P genotyping.modules:depth"
else
    bamFileArg="-I ${bamFileList}"
fi

genomeMaskArg="-genomeMaskFile ${svMaskFile}"
if [ "${useLcMask}" == "true" ]; then
    genomeMaskArg="${genomeMaskArg} -genomeMaskFile ${lcMaskFile}"
fi

for attempt in {1..3}; do
    echo $(date +"[%b %d %H:%M:%S] Running genotyping for partition ${partitionName}, attempt=${attempt}")

java -cp ${SV_CLASSPATH} -Xmx4g \
    org.broadinstitute.sv.main.SVGenotyper \
    -T SVGenotyperWalker \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -tempDir ${SV_TMPDIR} \
    -R ${referenceFile} \
    -md ${mdPath} \
    ${genomeMaskArg} \
    ${bamFileArg} \
    -disableGATKTraversal true \
    -vcf ${vcfFile} \
    ${extraArgs} \
    ${parameterArgs} \
    -runDirectory ${runDir} \
    -O ${outputVcf}

    rc=$?
    echo $(date +"[%b %d %H:%M:%S] Genotyping completed, rc=${rc}")

    if [ "${rc}" -eq 0 ]; then
        break
    fi
done

if [ "${rc}" -ne 0 ]; then
    echo "Error running genotyping" && exit 1
fi

echo "Archiving output..."
tar -cvzf gtrun.tar.gz ${runDir}

