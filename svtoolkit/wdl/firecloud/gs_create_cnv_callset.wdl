task parseBatchInfo {
    File batchInfo

    command {
        cut -f1 ${batchInfo} > "batch.list"
        cut -f1,2 ${batchInfo} > "batch_md_path.map"
        cut -f4 ${batchInfo} > "batch_cnv_calls_path.list"
    }

    output {
        Array[String] batchList = read_lines("batch.list")
        Map[String, String] batchMdPathMap = read_map("batch_md_path.map")
        Array[String] batchCnvCallsList = read_lines("batch_cnv_calls_path.list")
    }

    runtime {
        docker: "skashin/genome-strip:latest"
    }
}

task mergeCnvCalls {

    Array[File] cnvCallsList
    File referenceBundle
 
    Int diskSize
    Int numPreempt

    command {
        echo ${sep=' ' cnvCallsList} | sed 's/ /\n/g' > cnv_calls.list
        $SV_DIR/scripts/firecloud/cnv/merge_cnv_calls.sh cnv_calls.list ${referenceBundle}
    }

    output {
        File mergedSitesVcf = "gs_cnv.sites.vcf.gz"
    }

    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task setupGenotyping {
    File vcfFile
    Int parallelRecords

    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/compute_vcf_partitions.sh ${vcfFile} ${parallelRecords} partitions.dat
        cut -f 1 partitions.dat > partitions.list
    }

    output {
        File partitionFile = "partitions.dat"
        Array[String] partitionList = read_lines("partitions.list")
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        preemptible: "${numPreempt}"
    }
}

task runParallelGenotyper {
    File vcfFile
    String batch
    File partitionFile
    String partitionName
    File mdPath
    File referenceBundle
    File credentialsKeyFile

    Int memory
    Int numPreempt
   
    Int diskSize = round(size(mdPath, "G")) + 30

    command {
        cat ${partitionFile} | awk -v partitionName=${partitionName} '$1 == partitionName' | cut -f 2 > partition.arg

        $SV_DIR/scripts/firecloud/genotyping/run_parallel_genotyper.sh ${vcfFile} ${partitionName} "$(cat partition.arg)" ${mdPath} NULL ${referenceBundle} true ${credentialsKeyFile}
    }

    output {
        String batchList = "${batch}"
        File partitionVcf = "gtrun/${partitionName}.genotypes.vcf.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task computeBatchVcfList {
    String batch
    Array[String] batchList
    Array[String] partitionVcfList

    Int numPreempt

    command {
        echo ${sep=' ' batchList} | sed 's/ /\n/g' > batch.list
        echo ${sep=' ' partitionVcfList} | sed 's/ /\n/g' > vcf.list
        paste batch.list vcf.list | awk -v batch=${batch} '$1 == batch' | cut -f 2 > batch_vcf.list
    }

    output {
        Array[String] batchVcfList = read_lines("batch_vcf.list")
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        preemptible: "${numPreempt}"
    }
}

task mergeBatchPartitions {
    Array[File] partitionVcfList
    File referenceBundle

    Int diskSize
    Int numPreempt

    # Use size(partitionVcfList), once it's available
    Int adjustedDiskSize = diskSize + 5
    
    command {
        echo ${sep=' ' partitionVcfList} | sed 's/ /\n/g' > vcf.list

        source $SV_DIR/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

        java -cp $SV_CLASSPATH -Xmx4g \
            org.broadinstitute.sv.apps.VCFMerge \
            -R $referenceFile \
            -vcf vcf.list \
            -includeInfoTag END \
            -includeInfoTag GSELENGTH \
            -includeInfoTag SVTYPE \
            -O gs_cnv.genotypes.vcf.gz \
            || exit 1
    }

    output {
        File batchVcf = "gs_cnv.genotypes.vcf.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${adjustedDiskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task mergeBatches {
    Array[File] batchVcfList
    File referenceBundle

    Int cpu
    Int memory
    Int diskSize
    
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
            -O gs_cnv.genotypes.vcf.gz \
            || exit 1
    }

    output {
        File mergedGenotypesVcf = "gs_cnv.genotypes.vcf.gz"
        File mergedGenotypesVcfIndex = "gs_cnv.genotypes.vcf.gz.tbi"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        cpu: "${cpu}"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
    }
}

task setupRedundancyFiltering {

    File vcfFile
    File referenceBundle

    Int windowSize

    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/common/create_genome_partitions.sh ${vcfFile} ${windowSize} intervals.list ${referenceBundle}
    }

    output {
        Array[String] intervalsList = read_lines("intervals.list")
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task filterRedundantSites {
    File vcfFile
    File vcfFileIndex
    String interval
    File referenceBundle

    Int diskSize
    Int numPreempt
    
    Int adjustedDiskSize = diskSize + round(size(vcfFile, "G"))

    command {
        source $SV_DIR/scripts/firecloud/gs_extract_reference.sh ${referenceBundle} || exit 1

        java -cp $SV_CLASSPATH -Xmx4g \
            org.broadinstitute.sv.apps.FilterRedundantSites \
            -R $referenceFile \
            -vcf ${vcfFile} \
            -L ${interval} \
            -O gs_cnv.genotypes.vcf.gz \
            || exit 1
    }

    output {
        File dedupedVcf = "gs_cnv.genotypes.vcf.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${adjustedDiskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task createFinalCallset {
    Array[File] vcfFileList
    File referenceBundle

    Int diskSize

    command {
        echo ${sep=' ' vcfFileList} | sed 's/ /\n/g' > vcf_files.list
        $SV_DIR/scripts/firecloud/cnv/create_final_callset.sh vcf_files.list ${referenceBundle}
    }

    output {
        File cnvCallset = "cnv_callset.tar.gz"
    }

    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "7GB"
        cpu: "2"
        disks: "local-disk ${diskSize} HDD"
        preemptible: 0
    }
}

workflow gs_create_cnv_callset_wf {
    File batchInfo
    Int genotypingParallelRecords
    File referenceBundle
    File credentialsKeyFile

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt

    call parseBatchInfo {
        input:
            batchInfo = batchInfo
    }

    call mergeCnvCalls {
        input:
            cnvCallsList = parseBatchInfo.batchCnvCallsList,
            referenceBundle = referenceBundle,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    call setupGenotyping {
        input:
            vcfFile = mergeCnvCalls.mergedSitesVcf,
            parallelRecords = genotypingParallelRecords,
            numPreempt = numPreempt
    }

    Array[Pair[String, String]] batchPartitionPairs = cross(parseBatchInfo.batchList, setupGenotyping.partitionList)
    scatter(pair in batchPartitionPairs) {
        String batch = pair.left
        call runParallelGenotyper {
            input: 
                vcfFile = mergeCnvCalls.mergedSitesVcf,
                batch = batch,
                partitionFile = setupGenotyping.partitionFile,
                partitionName = pair.right,
                mdPath = parseBatchInfo.batchMdPathMap[batch],
                referenceBundle = referenceBundle,
                credentialsKeyFile = credentialsKeyFile,
                memory = memory,
                numPreempt = numPreempt
        }
    }

    scatter(batch in parseBatchInfo.batchList) {
        call computeBatchVcfList {
            input:
                batch = batch,
                batchList = runParallelGenotyper.batchList,
                partitionVcfList = runParallelGenotyper.partitionVcf,
                numPreempt = numPreempt
        }
    }
    
    scatter(batchVcfList in computeBatchVcfList.batchVcfList) {
        call mergeBatchPartitions {
            input:
                partitionVcfList = batchVcfList,
                referenceBundle = referenceBundle,
                diskSize = diskSize,
                numPreempt = numPreempt
        }
    }

    call mergeBatches {
        input:
            batchVcfList = mergeBatchPartitions.batchVcf,
            referenceBundle = referenceBundle
    }

    call setupRedundancyFiltering {
        input:
            vcfFile = mergeCnvCalls.mergedSitesVcf,
            referenceBundle = referenceBundle,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    scatter(interval in setupRedundancyFiltering.intervalsList) {
        call filterRedundantSites {
            input:
                vcfFile = mergeBatches.mergedGenotypesVcf,
                vcfFileIndex = mergeBatches.mergedGenotypesVcfIndex,
                interval = interval,
                referenceBundle = referenceBundle,
                diskSize = diskSize,
                numPreempt = numPreempt
        }
    }

    call createFinalCallset {
        input:
            vcfFileList = filterRedundantSites.dedupedVcf,
            referenceBundle = referenceBundle
    } 

    output {
        File mergedVcf = mergeBatches.mergedGenotypesVcf
        File cnvCallset = createFinalCallset.cnvCallset
    }
}

