task computeOverlappingVcfs {
    File partitionFile
    Array[File] callSpansList
    File vcfListFile
    File vcfIdxListFile

    command {
        printf "P%04d\n" $(cut -f4 ${partitionFile} | sort -uk1n) > partitions.list
        $SV_DIR/scripts/firecloud/joint/compute_partition_intervals.sh ${partitionFile} intervals.dat || exit 1

        echo ${sep=' ' callSpansList} | sed 's/ /\n/g' | xargs -a - cat > spans.dat

        $SV_DIR/scripts/firecloud/joint/compute_overlapping_vcfs.sh ${partitionFile} spans.dat ${vcfListFile} overlapping_vcf.dat || exit 1
        $SV_DIR/scripts/firecloud/joint/compute_overlapping_vcfs.sh ${partitionFile} spans.dat ${vcfIdxListFile} overlapping_vcf_idx.dat || exit 1
    }

    output {
        Array[String] partitionList = read_lines("partitions.list")
        Array[Array[String]] intervalList = read_tsv("intervals.dat")
        Array[Array[String]] overlappingVcfList = read_tsv("overlapping_vcf.dat")
        Array[Array[String]] overlappingVcfIdxList = read_tsv("overlapping_vcf_idx.dat")
    }

    runtime {
        docker: "skashin/genome-strip:latest"
    }
}

task filterRedundantCalls {
    File sampleListFile
    String partition
    Array[String] intervalList
    Array[File] vcfList
    Array[File] vcfIdxList
    File referenceBundle

    Int memory
    Int numThreads
    Int diskSize
    Int numPreempt
    
    Int adjustedDiskSize = length(vcfList) * diskSize + 30

    String filteredVcfName = "${partition}.filtered.genotypes.vcf.gz"

    command {
        echo ${sep=' ' intervalList} | sed 's/ /\n/g' | sed -e 's/:/\t/' -e 's/-/\t/' > interval.list
        echo ${sep=' ' vcfList} | sed 's/ /\n/g' > vcf.list
        $SV_DIR/scripts/firecloud/joint/filter_redundant_calls.sh interval.list vcf.list ${sampleListFile} ${filteredVcfName} ${referenceBundle} || exit 1
    }

    output {
    	File filteredVcf = "${filteredVcfName}"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        cpu: "${numThreads}"
        memory: "${memory}GB"
        disks: "local-disk ${adjustedDiskSiz} HDD"
        preemptible: "${numPreempt}"
    }
}

task createVcfListFile {
    Array[String] vcfList

    command {
    }

    output {
        File vcfListFile = write_lines(vcfList)
    }

    runtime {
        docker: "skashin/genome-strip:latest"
    }
}

workflow gs_filter_redundant_calls_wf {
    File partitionFile
    File sampleListFile
    File callSpansFile
    File vcfList
    File vcfIdxList
    File referenceBundle

    Array[String] callSpansList = read_lines(callSpansFile)

    call computeOverlappingVcfs {
        input:
            partitionFile = partitionFile,
            callSpansList = callSpansList,
            vcfListFile = vcfList,
            vcfIdxListFile = vcfIdxList
    }

    scatter (idx in range(length(computeOverlappingVcfs.partitionList))) {
        call filterRedundantCalls {
            input:
                sampleListFile = sampleListFile,
                partition = computeOverlappingVcfs.partitionList[idx],
                intervalList = computeOverlappingVcfs.intervalList[idx],
                vcfList = computeOverlappingVcfs.overlappingVcfList[idx],
                vcfIdxList = computeOverlappingVcfs.overlappingVcfIdxList[idx],
                referenceBundle = referenceBundle
        }
    }

    call createVcfListFile {
        input:
            vcfList = filterRedundantCalls.filteredVcf
    }

    output {
        File filteredVcfListFile = createVcfListFile.vcfListFile
    }
}

