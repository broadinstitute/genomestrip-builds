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

class CNVDiscoveryStage9 extends CNVDiscoveryStageBase {

    @Argument(shortName="adjacentMergedSitesVcf", required=true, doc="Adjacent merged sites output vcf file")
    var adjacentMergedSitesVcf: File = null

    def script = {
        addCommand(new MergeAdjacentSites)
        addCommand(new TouchSentinelFile(sentinelFile, adjacentMergedSitesVcf))
    }

    class MergeAdjacentSites extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.discovery.MergeDepthScannerSites"
        this.outputFile = adjacentMergedSitesVcf
        val evalDirectory = new File(runDirectory, "eval")
        createDirectory(evalDirectory)

        commandArguments +=
            required(" -R ", referenceFile) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            required(" -vcf ", vcfFile) +
            optional(" -reportFile ", new File(evalDirectory, "MergeAdjacentSites.report.dat")) +
            optional(" -outputVariants ", "ALL") +
            optional(" -adjacentSiteDistance ", "1000000") +
            optional(" -duplicateScoreThresholdMax ", "0")
    }
}
