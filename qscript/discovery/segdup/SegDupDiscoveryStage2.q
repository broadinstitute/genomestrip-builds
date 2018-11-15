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

class SegDupDiscoveryStage2 extends SVQScript {

    @Argument(shortName="sequenceName", required=true, doc="Sequence name")
    var sequenceName: String = null

    @Input(shortName="vcf", required=true, doc="Input VCF file")
    var vcfFile: File = null

    @Output(shortName="sentinelFile", required=true, doc="Sentinel fFile to be produced by this stage")
    var sentinelFile: File = null

    var genotypesVcfFile: File = null

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {

        genotypesVcfFile = swapExt(vcfFile, "sites.vcf", "genotypes.vcf")
        addCommand(new SVGenotyperScript)
        val copyNumberReportFile = addCommand(new CopyNumberClassAnnotator(genotypesVcfFile))
        val classifyVariants = addCommand(new ClassifySelectedVariants(genotypesVcfFile, copyNumberReportFile))
        //val gzipVcfFile = addCommand(new GzipFile(genotypesVcfFile))
        val deleteAuxFiles = addCommand(new DeleteAuxiliaryFiles(classifyVariants))
        add(new TouchSentinelFile(sentinelFile, deleteAuxFiles, new File(jobLogDirectory, "TouchSentinelFile")))
    }

    class SVGenotyperScript extends QueueScriptCommand {
        this.dependsOnFile +:= vcfFile
        this.commandResultFile = genotypesVcfFile
        override def scriptName = "SVGenotyper2.q"
        override def commandLine =
            super.commandLine +
            repeat(" -configFile ", parameterFiles) +
            required(" -P ", "depth.readCountCacheIgnoreGenomeMask:true") +
            required(" -P ", "genotyping.modules:depth") +
            repeat(" -P ", parameterList) +
            optional(" -R ", referenceFile) +
            repeat(" -genomeMaskFile ", genomeMaskFiles) +
            optional(" -ploidyMapFile ", ploidyMapFile) +
            optional(" -genderMapFile ", getGenderMapFileList) +
            optional(" -md ", metaDataLocation) +
            flag(" -disableGATKTraversal ", disableGATKTraversal) +
            required(" -vcf ", vcfFile) +
            required(" -parallelRecords ", "100") +
            repeat(" -I ", bamInputs) +
            required(" -O ", genotypesVcfFile) +
            required(" -runDirectory ", runDirectory)
    }

    class CopyNumberClassAnnotator(genotypesVcfFile: File)
        extends SVAnnotator(genotypesVcfFile, null) {
        val evalDirectory = new File(runDirectory, "eval")
        val copyNumberReportFile = new File(evalDirectory, "CopyNumberClass.report.dat")
        this.commandResultFile = copyNumberReportFile

        commandArguments +=
            required(" -A ", "CopyNumberClass") +
            required(" -writeReport ", "true") +
            required(" -reportFile ", copyNumberReportFile)
    }

    class ClassifySelectedVariants(genotypesVcfFile:File, copyNumberReportFile: File) extends SimpleCommand {
        this.dependsOnFile +:= copyNumberReportFile
        this.outputFile = new File(jobLogDirectory, "ClassifySelectedVariants")

        val evalDirectory = new File(runDirectory, "eval")
        commandArguments =
            required("scripts/classify_selected_variants.sh") +
            required(evalDirectory) +
            required(genotypesVcfFile) +
            required(copyNumberReportFile)
    }

    class GzipFile(fileToCOmpress: File) extends SimpleCommand {
        this.dependsOnFile +:= fileToCOmpress
        this.commandResultFile = new File(fileToCOmpress.getPath + ".gz")

        override def description = "GzipFile " + fileToCOmpress.getPath
        commandArguments =
            "gzip " + fileToCOmpress
    }

    class DeleteAuxiliaryFiles(depends: File) extends SimpleCommand {
        this.dependsOnFile +:= depends
        this.outputFile = new File(jobLogDirectory, "DeleteAuxiliaryFiles")

        override def description = "DeleteAuxiliaryFiles seq_" + sequenceName
        commandArguments =
            "find " + runDirectory +
            "  -name '*' -regex '.*/P[0-9]*.*' " +
            "  -exec rm -f {} \\;"
    }
}
