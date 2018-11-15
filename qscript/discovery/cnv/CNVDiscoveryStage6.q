/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2013 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.sv.discovery.cnv

import java.io.File
import java.io.PrintWriter
import scala.io.Source

class CNVDiscoveryStage6 extends CNVDiscoveryStageBase {

    @Input(shortName="selectedSamplesList", required=true, doc="Filtered list of samples")
    var selectedSamplesList: File = null

    @Output(shortName="selectedSamplesMergedHeadersBam", fullName="selectedSamplesMergedHeadersBam", required=true,
            doc="BAM file containing that is created by merging the headers from the selected samples' BAM files")
    var selectedSamplesMergedHeadersBam: File = null

    // BAM creation
    def bamHeadersDirectory = {
        new File(runDirectory, "bam_headers")
    }

    // In Stage6 we create the headers-only bam file for the discovery samples
    // Then we rerun the CopyNumberClass & SiteFilters annotators on Stage4 per-chromosome vcf files using the discovery samples
    def script = {
        val mergedBam = addCommand(new CreateMergedBamHeaders())
        val mergedBamIndex = addCommand(new CreateBamIndex(mergedBam))
        val variantSiteReports = addCommand(new CreateVariantSiteReports(vcfFile, mergedBamIndex))
        val classifyVariants = addCommand(new ClassifySelectedVariants(List(variantSiteReports)))
        addCommand(new TouchSentinelFile(sentinelFile, classifyVariants))
    }

    class CreateMergedBamHeaders extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.ExtractBAMSubset"
        this.inputFile = bamInputs
        this.outputFile = selectedSamplesMergedHeadersBam
        commandArguments +=
            required(" -L ", "NONE") +
            required(" -sample ", selectedSamplesList)
    }

    class CreateVariantSiteReports(genotypesVcfFile: File, mergedBamIndex: File)
        extends SVAnnotator(genotypesVcfFile, null) {
        val evalDirectory = new File(runDirectory, "eval")
        val siteFiltersReportFile = new File(evalDir, "GenotypeSiteFilters.report.dat")

        this.dependsOnFile :+= mergedBamIndex
        this.commandResultFile = siteFiltersReportFile

        commandArguments +=
            required(" -A ", "CopyNumberClass") +
            required(" -A ", "SiteFilters") +
            required(" -sample ", selectedSamplesList) +
            required(" -writeReport ", "true") +
            required(" -reportFileMap ", "SiteFilters:" + siteFiltersReportFile.getPath)
    }
}
