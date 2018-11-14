#!/bin/bash

referenceFile=$1

ploidyMapFile=$(echo ${referenceFile} | sed 's/.fasta$/.ploidymap.txt/')
svMaskFile=$(echo ${referenceFile} | sed 's/.fasta$/.svmask.fasta/')
lcMaskFile=$(echo ${referenceFile} | sed 's/.fasta$/.lcmask.fasta/')

parameterFile=$(echo ${referenceFile} | sed 's/.fasta$/.gsparams.txt/')
parameterFileArg=
if [ -e "${parameterFile}" ]; then
    parameterFileArg="-configFile ${parameterFile}"
fi
