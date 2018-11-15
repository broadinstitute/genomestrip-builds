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

// Initial windows scan
class CNVDiscoveryStage1 extends CNVDiscoveryStageBase {

    @Argument(shortName="intervalList", required=false, doc="Interval(s) or .list for which intervals to process (default is whole genome)")
    var intervalList: List[String] = null

    @Argument(shortName="tilingWindowSize", required=true, doc="Tiling Window Size")
    var tilingWindowSize: Int = _

    @Argument(shortName="tilingWindowOverlap", required=true, doc="Tiling Window Overlap")
    var tilingWindowOverlap: Int = _

    @Argument(shortName="maximumReferenceGapLength", required=false, doc="Max reference gap")
    var maximumReferenceGapLength: java.lang.Integer = null

    @Argument(shortName="testMaximumWindowCount", required=false, doc="For testing, max # windows to generate")
    var testMaximumWindowCount: java.lang.Integer = null

    @Argument(shortName="scannedWindowsVcfFile", required=true, doc="Scanned Windows output file")
    var scannedWindowsVcfFile: File = null


    def script = {
        addCommand(new WindowScanner)
        addCommand(new TouchSentinelFile(sentinelFile, scannedWindowsVcfFile))
    }

    class WindowScanner extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.discovery.SVDepthScanner"
        this.outputFile = scannedWindowsVcfFile
        commandArguments +=
            required(" -R ", referenceFile) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            repeat(" -L ", intervalList) +
            required(" -tilingWindowSize ", tilingWindowSize) +
            required(" -tilingWindowOverlap ", tilingWindowOverlap) +
            optional(" -testMaximumWindowCount ", testMaximumWindowCount) +
            optional(" -maximumReferenceGapLength ", maximumReferenceGapLength)
    }
}
