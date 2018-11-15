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
import java.io.PrintWriter
import org.broadinstitute.sv.qscript.SVQScript
import scala.io.Source

class SegDupDiscoveryPipeline extends SVQScript {

    @Input(shortName="segdupFile", required=true, doc="File of annotated segmental duplications from UCSC browser")
    var segdupFile: File = null

    @Input(shortName="sequenceMapFile", required=false, doc="Sequence map file")
    var sequenceMapFile: File = null

    var siteVcfFile: File = null
    var genotypesVcfFile: File = null

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        val vcfsDirectory = new File(runDirectory, "vcfs")
        createDirectory(vcfsDirectory)
        siteVcfFile = new File(vcfsDirectory, "segdup_scan.sites.vcf")
        genotypesVcfFile = new File(runDirectory, "segdup_scan.genotypes.vcf")

        addCommand(new GenerateSegdupSites)
        addCommand(new SVGenotyperScript)
    }

    class GenerateSegdupSites extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.discovery.GenerateSegdupSites"
        this.outputFile = siteVcfFile
        commandArguments +=
            required(" -R ", referenceFile) +
            required(" -segdupFile ", segdupFile) +
            optional(" -sequenceMapFile ", sequenceMapFile)
    }

    class SVGenotyperScript extends QueueScriptCommand {
        this.dependsOnFile +:= siteVcfFile
        this.commandResultFile = genotypesVcfFile
        override def scriptName = "SVGenotyper2.q"
        override def commandLine =
            super.commandLine +
            optional(" -configFile ", configFile) +
            required(" -P ", "genotyping.modules:depth") +
            repeat(" -P ", parameterList) +
            optional(" -R ", referenceFile) +
            optional(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
            repeat(" -genderMapFile ", getGenderMapFileList) +
            repeat(" -md ", metaDataLocationList) +
            flag(" -disableGATKTraversal ", disableGATKTraversal) +
            required(" -vcf ", siteVcfFile) +
            required(" -parallelRecords ", "100") +
            repeat(" -I ", bamInputs) +
            required(" -O ", genotypesVcfFile) +
            required(" -runDirectory ", runDirectory)
    }

    // TBD: Replace this functionality
    class ClassifySelectedVariants(copyNumberReportFile: File) extends SimpleCommand {
        this.dependsOnFile +:= copyNumberReportFile
        this.outputFile = new File(jobLogDirectory, "ClassifySelectedVariants")

        val evalDirectory = new File(runDirectory, "eval")
        commandArguments =
            required("scripts/classify_selected_variants.sh") +
            required(evalDirectory) +
            required(genotypesVcfFile) +
            required(copyNumberReportFile)
    }
}
