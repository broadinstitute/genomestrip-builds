import "https://api.firecloud.org/ga4gh/v1/tools/GenomeStrip:RunDeletionPipelineOverInterval/versions/2/plain-WDL/descriptor" as runDelDiscoveryWorkflow

task parseSampleInfo {
    File sampleInfo

    String docker
       
    command {
        cut -f4 ${sampleInfo} | tail -n +2 | sort -u -k1,1 > batch.list
    }

    output {
        Array[String] batchList = read_lines("batch.list")
    }

    runtime {
        docker: "${docker}"
    }
}

task computeDiscoveryPartitions {
    String searchInterval
    File referenceBundle

    String docker
    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/del/compute_partitions.sh ${searchInterval} partitions.dat ${referenceBundle} || exit 1
        cut -f 1 partitions.dat > partitions.list
        cut -f 5 partitions.dat > search_interval.list
    }

    output {
        File partitionFile = "partitions.dat"
        Array[String] partitionList = read_lines("partitions.list")
        Array[String] searchIntervalList = read_lines("search_interval.list")
    }
    
    runtime {
        docker: "${docker}"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task mergeBatchDataFiles {
    Array[File] batchDiscoveryPairsFileList
    Array[File] batchGenotypesPairsFileList

    String docker
    Int memory
    Int diskSize
    Int numPreempt
    
    command {
        echo ${sep=' ' batchDiscoveryPairsFileList} | sed 's/ /\n/g' > discovery_pairs_file.list
        $SV_DIR/scripts/firecloud/del/merge_pair_files.sh "discovery" discovery_pairs_file.list gs_dels.discovery.pairs.gz || exit 1

        echo ${sep=' ' batchGenotypesPairsFileList} | sed 's/ /\n/g' > genotypes_pairs_file.list
        $SV_DIR/scripts/firecloud/del/merge_pair_files.sh "genotypes" genotypes_pairs_file.list gs_dels.genotypes.pairs.gz || exit 1
    }

    output {
        File discoveryPairsFile = "gs_dels.discovery.pairs.gz"
        File discoveryPairsIndexFile = "gs_dels.discovery.pairs.gz.tbi"
        File genotypesPairsFile = "gs_dels.genotypes.pairs.gz"
        File genotypesPairsIndexFile = "gs_dels.genotypes.pairs.gz.tbi"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task mergeIntervalCallsets {
    Array[String] batchList
    Array[File] callsetList
    Array[File] homologyFileList
    Array[File] batchDiscoveryPairsFileList
    Array[File] batchDiscoveryPairsIndexFileList
    Array[File] batchGenotypesPairsFileList
    Array[File] batchGenotypesPairsIndexFileList

    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt
    
    command {
        echo ${sep=' ' batchList} | sed 's/ /\n/g' > batch.list
        echo ${sep=' ' callsetList} | sed 's/ /\n/g' > vcf.list
        echo ${sep=' ' homologyFileList} | sed 's/ /\n/g' > homology_file.list
        echo ${sep=' ' batchDiscoveryPairsFileList} | sed 's/ /\n/g' > discovery_file.list
        echo ${sep=' ' batchDiscoveryPairsIndexFileList} | sed 's/ /\n/g' > discovery_idx_file.list
        echo ${sep=' ' batchGenotypesPairsFileList} | sed 's/ /\n/g' > genotypes_file.list
        echo ${sep=' ' batchGenotypesPairsIndexFileList} | sed 's/ /\n/g' > genotypes_idx_file.list

        source $SV_DIR/scripts/firecloud/del/merge_interval_callsets.sh batch.list vcf.list homology_file.list discovery_file.list discovery_idx_file.list genotypes_file.list genotypes_idx_file.list ${referenceBundle}

    }

    output {
        File delCallset = "del_output.tar.gz"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

workflow gs_create_del_callset_by_interval {
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

    call computeDiscoveryPartitions {
        input:
            searchInterval = searchInterval,
            referenceBundle = referenceBundle,
            docker = docker,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    scatter(searchInterval in computeDiscoveryPartitions.searchIntervalList) {
        call runDelDiscoveryWorkflow.gs_run_del_pipeline_by_interval as runDelPipeline {
            input:
                sampleInfo = sampleInfo,
                searchInterval = searchInterval,
                parameterArgs = parameterArgs,
                referenceBundle = referenceBundle,
                repeatTrackFile = repeatTrackFile,
                docker = docker,
                memory = memory,
                diskSize = diskSize,
                numPreempt = numPreempt
        }
    }

    Array[Array[File]] discoveryPairsFileList = transpose(runDelPipeline.discoveryPairsFileList)
    Array[Array[File]] genotypesPairsFileList = transpose(runDelPipeline.genotypesPairsFileList)

    scatter(batchIdx in range(length(parseSampleInfo.batchList))) {
        call mergeBatchDataFiles {
            input:
                batchDiscoveryPairsFileList = discoveryPairsFileList[batchIdx],
                batchGenotypesPairsFileList = genotypesPairsFileList[batchIdx],
                docker = docker,
                memory = memory,
                diskSize = diskSize,
                numPreempt = numPreempt
        }
    }

    call mergeIntervalCallsets {
        input:
            batchList = parseSampleInfo.batchList,
            callsetList = runDelPipeline.delCallset,
            homologyFileList = runDelPipeline.homologyFile,
            batchDiscoveryPairsFileList = mergeBatchDataFiles.discoveryPairsFile,
            batchDiscoveryPairsIndexFileList = mergeBatchDataFiles.discoveryPairsIndexFile,
            batchGenotypesPairsFileList = mergeBatchDataFiles.genotypesPairsFile,
            batchGenotypesPairsIndexFileList = mergeBatchDataFiles.genotypesPairsIndexFile,
            referenceBundle = referenceBundle,
            docker = docker,
            memory = memory,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    output {
        File delCallset = mergeIntervalCallsets.delCallset
    }
}

