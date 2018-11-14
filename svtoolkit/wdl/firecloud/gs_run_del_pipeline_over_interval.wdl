task parseSampleInfo {
    File sampleInfo

    String docker
       
    command {
        cut -f4 ${sampleInfo} | tail -n +2 | sort -u -k1,1 > batch.list
        cut -f4,5 ${sampleInfo} | tail -n +2 | sort -u -k1,1 | cut -f2 > md_path.list

        tail -n +2 ${sampleInfo} | sort -t$'\t' -k4,4 > sample_info.dat
        $SV_DIR/scripts/firecloud/common/create_grouped_tsv.sh sample_info.dat 3 4 > batch_crams.tsv
    }

    output {
        Array[String] batchList = read_lines("batch.list")
        Array[String] mdPathList = read_lines("md_path.list")
        Array[Array[String]] batchCramFiles = read_tsv("batch_crams.tsv")
    }

    runtime {
        docker: "${docker}"
    }
}

task runBatchDiscovery {
    String searchInterval
    String mdPath
    Array[String] bamFileList
    String parameterArgs
    File referenceBundle
    File repeatTrackFile

    String docker
    Int memory
    Int diskSize
    Int numPreempt
   
    command {
        echo ${sep=' ' bamFileList} | sed 's/ /\n/g' > bam_file.list || exit 1

        $SV_DIR/scripts/firecloud/del/run_complete_discovery.sh "gs_dels" ${searchInterval} ${mdPath} bam_file.list "${parameterArgs}" ${referenceBundle} ${repeatTrackFile} || exit 1

        # Store the pairs bam file separately from the other output
        mv del_output/gs_dels.discovery.pairs.bam* .

        echo "Archiving output..."
        tar -cvzf del_output.tar.gz del_output
    }

    output {
        File delOutputArchive = "del_output.tar.gz"
        File discoveryPairsFile = "del_output/gs_dels.discovery.pairs.dat"
        File homologyFile = "del_output/gs_dels.discovery.homology.dat"
        File bamFile = "gs_dels.discovery.pairs.bam"
        File bamFileIndex = "gs_dels.discovery.pairs.bam.bai"

        File vcfFile = "gs_dels.sites.vcf.gz"
    }
    
    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task mergeDelCalls {
    Array[File] vcfFileList
    Array[File] homologyFileList
    File referenceBundle

    String docker
    Int diskSize
    Int numPreempt

    command {
        echo ${sep=' ' vcfFileList} | sed 's/ /\n/g' > vcf_file.list
        $SV_DIR/scripts/firecloud/del/merge_del_calls.sh vcf_file.list ${referenceBundle} || exit 1

        echo ${sep=' ' homologyFileList} | sed 's/ /\n/g' > homology_file.list
        $SV_DIR/scripts/firecloud/del/merge_homology_files.sh homology_file.list gs_dels.homology.gz || exit 1
    }

    output {
        File mergedSitesVcf = "gs_dels.sites.vcf.gz"
        File mergedHomologyFile = "gs_dels.homology.gz"
        File mergedHomologyIndexFile = "gs_dels.homology.gz.tbi"
    }

    runtime {
        docker: "${docker}"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task runBatchGenotyper {
    File vcfFile
    String mdPath
    File bamFile
    File bamFileIndex
    String parameterArgs
    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt
   
    Int adjustedDiskSize = diskSize + round(size(bamFile, "G"))

    command {
        $SV_DIR/scripts/firecloud/genotyping/run_genotyper.sh ${vcfFile} ${mdPath} ${bamFile} "${parameterArgs}" ${referenceBundle} false
    }

    output {
        File batchVcf = "gtrun/gs_dels.genotypes.vcf.gz"
        File genotypesPairsFile = "gtrun/gs_dels.genotypes.pairs.dat"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${adjustedDiskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task mergeBatches {
    Array[File] batchVcfList
    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt
    
    command {
        echo ${sep=' ' batchVcfList} | sed 's/ /\n/g' > vcf.list

        source $SV_DIR/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

        java -cp $SV_CLASSPATH -Xmx4g \
            org.broadinstitute.sv.apps.VCFMerge \
            -R $referenceFile \
            -vcf vcf.list \
            -includeInfoTag END \
            -includeInfoTag GSELENGTH \
            -includeInfoTag SVTYPE \
            -includeInfoTag SVLEN \
            -includeInfoTag CIPOS \
            -includeInfoTag CIEND \
            -includeInfoTag GSBKPT \
            -includeInfoTag GSNPAIRS \
            -includeInfoTag GSNSAMPLES \
            -O gs_del.genotypes.vcf.gz \
            || exit 1
    }

    output {
        File mergedGenotypesVcf = "gs_del.genotypes.vcf.gz"
        File mergedGenotypesVcfIndex = "gs_del.genotypes.vcf.gz.tbi"
    }
    
    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task filterRedundantSites {
    File vcfFile
    File vcfFileIndex
    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt

    Int adjustedDiskSize = diskSize + round(size(vcfFile, "G"))
    
    command {
        source $SV_DIR/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

        java -cp $SV_CLASSPATH -Xmx4g \
            org.broadinstitute.sv.apps.FilterRedundantSites \
            -R $referenceFile \
            -vcf ${vcfFile} \
            -O gs_del.genotypes.vcf.gz \
            || exit 1
    }

    output {
        File dedupedVcf = "gs_del.genotypes.vcf.gz"
    }
    
    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${adjustedDiskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task createFinalCallset {
    File vcfFile
    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/del/create_final_callset.sh ${vcfFile} ${referenceBundle}
    }

    output {
        File delCallset = "del_callset/final/del.genotypes.vcf.gz"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

workflow gs_run_del_pipeline_by_interval {
    File sampleInfo
    String searchInterval
    String parameterArgs
    File referenceBundle
    File repeatTrackFile

    String docker
    Int memory
    Int diskSize
    Int numPreempt

    call parseSampleInfo {
        input:        	
            sampleInfo = sampleInfo,
            docker = docker
    }

    scatter(batchIdx in range(length(parseSampleInfo.batchList))) {
        String batch = parseSampleInfo.batchList[batchIdx]
        call runBatchDiscovery {
            input:
                searchInterval = searchInterval,
                mdPath = parseSampleInfo.mdPathList[batchIdx],
                bamFileList = parseSampleInfo.batchCramFiles[batchIdx],
                parameterArgs = parameterArgs,
                referenceBundle = referenceBundle,
                repeatTrackFile = repeatTrackFile,
                docker = docker,
                memory = memory,
                diskSize = diskSize,
                numPreempt = numPreempt
        }
    }

    call mergeDelCalls {
        input:
            vcfFileList = runBatchDiscovery.vcfFile,
            homologyFileList = runBatchDiscovery.homologyFile,
            referenceBundle = referenceBundle,
            docker = docker,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    scatter(batchIdx in range(length(parseSampleInfo.batchList))) {
        call runBatchGenotyper {
            input: 
                vcfFile = mergeDelCalls.mergedSitesVcf,
                mdPath = parseSampleInfo.mdPathList[batchIdx],
                bamFile = runBatchDiscovery.bamFile[batchIdx],
                bamFileIndex = runBatchDiscovery.bamFileIndex[batchIdx],
                parameterArgs = parameterArgs,
                referenceBundle = referenceBundle,
                docker = docker,
                memory = memory,
                diskSize = diskSize,
                numPreempt = numPreempt
        }
    }

    call mergeBatches {
        input:
            batchVcfList = runBatchGenotyper.batchVcf,
            referenceBundle = referenceBundle,
            docker = docker,
            memory = memory,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    call filterRedundantSites {
        input:
            vcfFile = mergeBatches.mergedGenotypesVcf,
            vcfFileIndex = mergeBatches.mergedGenotypesVcfIndex,
            referenceBundle = referenceBundle,
            docker = docker,
            memory = memory,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    call createFinalCallset {
        input:
            vcfFile = filterRedundantSites.dedupedVcf,
            referenceBundle = referenceBundle,
            docker = docker,
            memory = memory,
            diskSize = diskSize,
            numPreempt = numPreempt
    } 

    output {
        File delCallset = createFinalCallset.delCallset
        Array[File] discoveryPairsFileList = runBatchDiscovery.discoveryPairsFile
        Array[File] genotypesPairsFileList = runBatchGenotyper.genotypesPairsFile
        File homologyFile = mergeDelCalls.mergedHomologyFile
        File homologyIndexFile = mergeDelCalls.mergedHomologyIndexFile
    }
}
