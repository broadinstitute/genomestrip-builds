#!/bin/bash

inputVcf=$1
outputVcf=$2
referenceFile=$3

if [ -z "${referenceFile}" ]; then
    echo "Usage: apply_default_filters.sh <inputVcf> <outputVcf> <referenceFile>" && exit 1
fi

java -Xmx4g -jar $GATK_JAR \
    -T VariantFiltration \
    -V ${inputVcf} \
    -o ${outputVcf} \
    -R $referenceFile \
    -filterName COVERAGE -filter "GSDEPTHCALLTHRESHOLD == \"NA\" || GSDEPTHCALLTHRESHOLD >= 1.0" \
    -filterName COHERENCE -filter "GSCOHPVALUE == \"NA\" || GSCOHPVALUE <= 0.01" \
    -filterName DEPTHPVAL -filter "GSDEPTHPVALUE == \"NA\" || GSDEPTHPVALUE >= 0.01" \
    -filterName DEPTH -filter "GSDEPTHRATIO == \"NA\" || GSDEPTHRATIO > 0.8 || (GSDEPTHRATIO > 0.63 && (GSMEMBPVALUE == \"NA\" || GSMEMBPVALUE >= 0.01))" \
    -filterName PAIRSPERSAMPLE -filter "GSNPAIRS <= 1.1 * GSNSAMPLES"

