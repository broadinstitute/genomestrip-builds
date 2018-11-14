task filterCalls {

    File cnvOutput
    String mdPath
    Int madRange
    File referenceBundle
    File credentialsKeyFile

    Int numPreempt

    Int diskSize = 2 * round(size(cnvOutput, "G")) + 30

    command {
        $SV_DIR/scripts/firecloud/cnv/filter_cnv_sites.sh ${cnvOutput} ${mdPath} ${madRange} ${referenceBundle} ${credentialsKeyFile}
    }

    output {
        File vpsReport = "eval/VariantsPerSample.report.dat"
        File selectedSites = "eval/SelectedSites.dat"
    }
    
    runtime {
        docker: "skashin/genome-strip:latest"
        disks: "local-disk ${diskSize} HDD"
        preemptible: "${numPreempt}"
    }
}

workflow gs_filter_cnv_sites_wf {
    call filterCalls
}

