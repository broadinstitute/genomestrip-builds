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

import org.broadinstitute.sv.qscript.SVQScript
import java.io.PrintWriter
import scala.io.Source
import org.broadinstitute.sv.util.io.ErrorCheckingPrintWriter

abstract class CNVDiscoveryStageBase extends SVQScript {

    @Input(shortName="vcf", required=false, doc="Input VCF file for this stage")
    var vcfFile: File = null

    @Argument(shortName="sequenceName", required=false, doc="Sequence name for this stage")
    var sequenceName: String = null

    @Argument(fullName="genotypeFilterFile", shortName="genotypeFilterFile", doc="File(s) conatining per-sample regions to be excluded from genotyping", required=false)
    var genotypeFilterFileList: List[File] = Nil

    @Argument(shortName="keepTempFiles", required=false, doc="Categories of temporary files to keep")
    var keepTempFiles: List[String] = Nil

    @Output(shortName="sentinelFile", fullName="sentinelFile", required=true, doc="Sentinel file to be produced by this stage")
    var sentinelFile: File = null

    def shouldDeleteTempFiles(category: String) : Boolean =
        !(keepTempFiles.contains(category) || keepTempFiles.contains("ALL"))

    class ClassifySelectedVariants(dependsOn: List[File]) extends InProcessFunction
            with LogDirectory with UnifiedInputOutput {
        val evalDirectory = new File(runDirectory, "eval")
        val copyNumberReportFile = new File(evalDirectory, "CopyNumberClass.report.dat")
        val siteFiltersReportFile = new File(evalDirectory, "GenotypeSiteFilters.report.dat")
        val selectedVariantsList = new File(evalDirectory, "SelectedVariants.list")

        this.dependsOnFile :::= dependsOn
        this.commandResultFile = selectedVariantsList

        def run = {
            val listWriter = new ErrorCheckingPrintWriter(selectedVariantsList)
            var src = Source.fromFile(siteFiltersReportFile)
            val passingVariants = src
                .getLines.filter(line => !line.startsWith("#"))
                .map(_.split("\t", -1))
                .filter(fields => fields(1).equals("PASS"))
                .map(fields => fields(0))
                .toSet
            src.close

            src = Source.fromFile(copyNumberReportFile)
            val lines = src.getLines.drop(1)
            lines.foreach(line => {
                val fields = line.split("\t")
                if (fields(5).toInt > 0 && passingVariants.contains(fields(0))) {
                    listWriter.println(fields(0))
                }
            })
            src.close
            listWriter.close
        }
    }
}
