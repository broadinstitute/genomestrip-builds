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

class CNVDiscoveryStage5 extends CNVDiscoveryStageBase {

    @Input(shortName="vpsReportsDirectory", required=true, doc="Directory where per-sequence VPS reports are located")
    var vpsReportsDirectory: File = null

    @Input(shortName="selectedSamplesList", required=true, doc="Output file for the selected samples")
    var selectedSamplesList: File = null

    @Argument(shortName="useAllSamplesInDiscovery", required=false, doc="Flag specifying whether to yse all the samples in discovery")
    var useAllSamplesInDiscovery: Boolean = false


    def script = {
        addCommand(new FilterSamples)
        val vpsReportPdf = addCommand(new MakeVpsPlot)
        addCommand(new TouchSentinelFile(sentinelFile, vpsReportPdf))
    }

    class FilterSamples extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.discovery.AnalyzeVPSReports"
        this.commandResultFile = selectedSamplesList
        val evalDir = new File(runDirectory, "eval")
        createDirectory(evalDir)
        commandArguments +=
            required(" -R ", referenceFile) +
            required(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
            required(" -vpsReportsDirectory ", vpsReportsDirectory) +
            required(" -outputDirectory ", selectedSamplesList.getParent()) +
            flag(" -useAllSamplesInDiscovery ", useAllSamplesInDiscovery)
    }

    class MakeVpsPlot extends SimpleCommand {
        val vpsReport = new File(selectedSamplesList.getParent(), "VariantsPerSample.report.dat")
        val vpsReportPdf = new File(selectedSamplesList.getParent(), "VariantsPerSample.report.pdf")

        this.dependsOnFile :+= selectedSamplesList
        this.commandResultFile = vpsReportPdf

        val plotVPSScriptName = "discovery/plot_vps.R"
        val plotVPSScriptPath = org.broadinstitute.sv.util.RUtilities.findRScript(plotVPSScriptName)
        commandArguments =
            required("Rscript") +
            required(plotVPSScriptPath) +
            required(vpsReport) +
            required(vpsReportPdf)
    }
}
