import "run_multi_batch_genotyper.wdl" as runGenotyperWorkflow

task createCallset {
    String callType
    File genotypesVcf
    File callsBedFile
    File? sampleListFile
    File genderMapFile
    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/segpipe/create_segpipe_callset.sh ${callType} ${genotypesVcf} ${callsBedFile} "${sampleListFile}" ${genderMapFile} ${referenceBundle}
    }

    output {
        File callsetVcf = "segpipe_${callType}.genotypes.vcf.gz"
    }

    runtime {
        docker: "${docker}"
        memory: "7GB"
        disks: "local-disk 40 HDD"
        preemptible: "${numPreempt}"
    }
}

workflow gs_create_segpipe_callset {
    String callType
    File batchInfo
    File callsVcfFile
    File callsBedFile
    File? sampleListFile
    Boolean useLcMask
    File genderMapFile
    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt

    call runGenotyperWorkflow.gs_run_multi_batch_genotyper as runGenotyper {
        input:
            vcfFile = callsVcfFile,
            batchInfo = batchInfo,
            useLcMask = useLcMask,
            referenceBundle = referenceBundle,
            docker = docker,
            memory = memory,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    call createCallset {
        input:        	
            callType = callType,
            genotypesVcf = runGenotyper.genotypesVcf,
            callsBedFile = callsBedFile,
            sampleListFile = sampleListFile,
            genderMapFile = genderMapFile,
            referenceBundle = referenceBundle,
            docker = docker,
            memory = memory,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    output {
        File callsetVcf = createCallset.callsetVcf
    }
}

