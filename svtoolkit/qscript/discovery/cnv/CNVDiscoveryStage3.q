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
import collection.JavaConversions

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._


class CNVDiscoveryStage3 extends CNVDiscoveryStageBase {

    @Argument(shortName="mergedVcfFile", required=true, doc="Merged Windows VCF file")
    var mergedVcfFile: File = null

    @Argument(shortName="duplicateScoreThresholdMin", required=false, doc="Minimum threshold for all matches for merging (default none)")
    var duplicateScoreThresholdMin: java.lang.Double = null

    @Argument(shortName="duplicateScoreThresholdMax", required=false, doc="Minimum threshold for best match for merging (default none)")
    var duplicateScoreThresholdMax: java.lang.Double = null

    @Argument(shortName="correlationThresholdMin", required=false, doc="Minimum positive correlation (r2) for all matches for merging (default none)")
    var correlationThresholdMin: java.lang.Double = null

    @Argument(shortName="correlationThresholdMax", required=false, doc="Minimum positive correlation (r2) for best match for merging (default none)")
    var correlationThresholdMax: java.lang.Double = null


    def script = {
        addCommand(new MergeDepthScannerSites(vcfFile))
        addCommand(new TouchSentinelFile(sentinelFile, mergedVcfFile))
    }

    class MergeDepthScannerSites(vcfFile: File) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.discovery.MergeDepthScannerSites"
        this.dependsOnFile :+= vcfFile
        this.outputFile = mergedVcfFile
        commandArguments +=
            optional(" -debug ", debug) +
            required(" -vcf ", vcfFile) +
            required(" -R ", referenceFile) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            required(" -reportFile ", new File(runDirectory, "MergeScannerSites.report.dat")) +
            required(" -outputVariants ", "ALL") +
            optional(" -duplicateScoreThresholdMin ", duplicateScoreThresholdMin) +
            optional(" -duplicateScoreThresholdMax ", duplicateScoreThresholdMax) +
            optional(" -correlationThresholdMin ", correlationThresholdMin) +
            optional(" -correlationThresholdMax ", correlationThresholdMax)
    }
}
