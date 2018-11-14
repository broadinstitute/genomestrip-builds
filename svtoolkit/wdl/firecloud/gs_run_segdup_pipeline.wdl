task runSegdup {
    File mdPath
    File segdupFile
    File sequenceMapFile
    File referenceBundle

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt

    Int adjustedDiskSize = diskSize + round(size(mdPath, "G")) + 1

    command {
        $SV_DIR/scripts/firecloud/segdup/run_segdup_pipeline.sh ${mdPath} ${segdupFile} ${sequenceMapFile} ${referenceBundle}
    }

    output {
        File segdupOutput = "segdup_output.tar.gz"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "${memory}GB"
        disks: "local-disk ${adjustedDiskSize} HDD"
        cpu: "${numThreads}"
        preemptible: "${numPreempt}"
    }
}

workflow gs_run_segdup_pipeline_wf {

    call runSegdup

    output {
        File segdupOutput = runSegdup.segdupOutput
    }
}

