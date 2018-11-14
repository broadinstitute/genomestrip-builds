task genomestrip_preprocessing {

    File bam_file
    File bam_index
    File genomestrip_reference
    String sample_id

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
 
    command {
        ${SV_DIR}/scripts/firecloud/gs_preprocess.sh ${sample_id} ${bam_file} ${genomestrip_reference}
    }
    
    output {
        File genomestrip_metadata = "${sample_id}.tar.gz"
    }

    runtime {
        docker: "skashin/genome-strip:r1736"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
        author: "Seva Kashin"
    }
}

workflow genomestrip_preprocessing_workflow {
    call genomestrip_preprocessing
}
