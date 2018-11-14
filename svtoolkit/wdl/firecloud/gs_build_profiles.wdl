task gs_build_profiles {

    Int binSize
    Array[File] mdArchives
    File gsReferenceBundle

    Int memory
    Int diskSize
    Int numThreads
    Int numPreempt
 
    String mdArchivesList = "md_archives.list"

    command {
        echo ${sep=' ' mdArchives} | sed 's/ /\n/g' > ${mdArchivesList}
        ${SV_DIR}/scripts/firecloud/gs_build_profiles.sh ${mdArchivesList} ${binSize} ${gsReferenceBundle}
    }

    output {
        File profiles = "profiles_${binSize}.tar.gz"
    }

    runtime {
        docker: "skashin/genome-strip:latest"
        memory: "${memory}GB"
        disks: "local-disk ${diskSize} HDD"
        cpu: "${numThreads}"
        preemptible: "${numPreempt}"
    }

    meta {
        author: "Seva Kashin"
    }
}

workflow gs_build_profiles_workflow {
    call gs_build_profiles
}

