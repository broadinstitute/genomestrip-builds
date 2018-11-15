/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2010 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
import java.io.File
import java.io.PrintWriter
import org.broadinstitute.sv.qscript.SVQScript
import scala.io.Source

class SegDupDiscoveryStage1 extends SVQScript {

    @Argument(shortName="sequenceName", required=true, doc="Sequence name")
    var sequenceName: String = null

    @Input(shortName="sequenceMapFile", required=true, doc="Sequence map file")
    var sequenceMapFile: String = null

    @Input(shortName="segdupFile", required=true, doc="Segdup file")
    var segdupFile: File = null

    @Argument(shortName="siteVcfFile", required=true, doc="Sites VCF file")
    var siteVcfFile: File = null

    @Output(shortName="sentinelFile", required=true, doc="Sentinel fFile to be produced by this stage")
    var sentinelFile: File = null

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        addCommand(new GenerateSegdupSites)
        add(new TouchSentinelFile(sentinelFile, scannedWindowsVcfFile, new File(jobLogDirectory, "TouchSentinelFile")))
    }

    class GenerateSegdupSites extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.playground.discovery.GenerateSegdupSites"
        this.outputFile = siteVcfFile
        commandArguments +=
            required(" -R ", referenceFile) +
            required(" -sequenceMapFile ", sequenceMapFile) +
            required(" -segdupFile ", segdupFile)
    }
}
