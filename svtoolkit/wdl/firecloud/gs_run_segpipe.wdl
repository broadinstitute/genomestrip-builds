task createProfile {
    String chrom
    Int binSize
    String mdPath
    Boolean useLcMask
    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/build_profiles.sh ${chrom} ${binSize} ${mdPath} ${referenceBundle} ${useLcMask}
    }

    output {
        File profileFile = "profiles_${binSize}/profile_seq_${chrom}_${binSize}.dat.gz"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

task runSegpipe {
    File profileFile
    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt

    command {
        $SV_DIR/scripts/firecloud/segpipe/run_segpipe.sh ${profileFile} ${referenceBundle}
    }

    output {
        File profileFile = "profiles_${binSize}/profile_seq_${chrom}_${binSize}.dat.gz"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

workflow gs_run_segpipe {
    String chrom
    Int binSize
    String mdPath
    Boolean useLcMask
    File referenceBundle

    String docker
    Int memory
    Int diskSize
    Int numPreempt

    call createProfile {
        input:
            chrom = chrom,
            binSize = binSize,
            mdPath = mdPath,
            useLcMask = useLcMask,
            referenceBundle = referenceBundle,
            docker = docker,
            memory = memory,
            diskSize = diskSize,
            numPreempt = numPreempt
    }

    call runSegpipe {
        input:        	
            profileFile = crearteProfile.profileFile,
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

