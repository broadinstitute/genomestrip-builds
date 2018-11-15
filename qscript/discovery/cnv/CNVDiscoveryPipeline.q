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

import java.io.PrintWriter
import org.broadinstitute.sv.commandline.ArgumentException;
import org.broadinstitute.sv.util.GenomeInterval

import org.broadinstitute.sv.qscript.SVQScript
import scala.io.Source

class CNVDiscoveryPipeline extends SVQScript {

    @Argument(shortName="tilingWindowSize", required=true, doc="Tiling Window Size")
    var tilingWindowSize: java.lang.Integer = null

    @Argument(shortName="tilingWindowOverlap", required=true, doc="Tiling Window Overlap")
    var tilingWindowOverlap: java.lang.Integer = null

    @Argument(shortName="maximumReferenceGapLength", required=true, doc="Max reference gap")
    var maximumReferenceGapLength: java.lang.Integer = null

    @Argument(shortName="boundaryPrecision", required=true, doc="Boundary refinement precision")
    var boundaryPrecision: java.lang.Integer = null

    @Argument(shortName="minimumRefinedLength", required=true, doc="Minimum refined length")
    var minimumRefinedLength: java.lang.Integer = null

    @Argument(shortName="maximumRefinedLength", required=false, doc="Maximum length of refined interval (effective length)")
    var maximumRefinedLength: java.lang.Integer = null

    @Argument(shortName="useAllSamplesInDiscovery", required=false, doc="Flag specifying whether to yse all the samples in discovery")
    var useAllSamplesInDiscovery: Boolean = false

    @Argument(shortName="produceAuxiliaryFiles", required=false, doc="Flag specifying whether to produce auxiliary files and keep them")
    var produceAuxiliaryFiles: Boolean = false

    @Argument(fullName="intervalList", shortName="intervalList", required=false, doc="Interval(s) or .list for which intervals to process (default is whole genome)")
    var intervalList: List[String] = null

    @Argument(fullName="genotypeFilterFile", shortName="genotypeFilterFile", doc="File(s) conatining per-sample regions to be excluded from genotyping", required=false)
    var genotypeFilterFileList: List[File] = Nil

    @Argument(shortName="lastStage", required=false, doc="For testing, last stage to be executed")
    var lastStage: Int = 12

    @Argument(shortName="genotypingParallelRecords", required=false, doc="The number of parallel records to use in genotyping")
    var genotypingParallelRecords: java.lang.Integer = 1000

    @Argument(shortName="maxConcurrentStageJobs", required=false, doc="The maximum number of parallel jobs allowed to run in each stage (default no limit)")
    var maxConcurrentStageJobs: java.lang.Integer = null

    @Argument(fullName="outputVcfName", shortName="outputVcfName", required=false, doc="Name of final vcf output file")
    var outputVcfName: String = "gs_cnv.genotypes.vcf.gz"

    @Argument(shortName="keepTempFiles", required=false, doc="Categories of temporary files to keep")
    var keepTempFiles: List[String] = Nil

    @Argument(shortName="testMaximumWindowCount", required=false, doc="For testing, max # windows to generate")
    var testMaximumWindowCount: java.lang.Integer = null

    lazy val sequenceIntervalMap: Map[String, List[String]] = {
        var map: Map[String, List[String]] = Map()
        val refSequenceSet = computeReferenceSequenceNames().toSet

        if (intervalList == null) {
            val rmdIntervalFile = getReferenceMetadata.getMetadataFile(null, "interval.list")
            if (rmdIntervalFile != null) {
                intervalList = List(rmdIntervalFile.getPath)
            }
        }
        if (intervalList == null) {
            map = refSequenceSet.map(x => x -> List(x)).toMap
        } else {
            intervalList = expandListFiles(intervalList)
            intervalList.foreach(interval => {
                val genomeInterval = GenomeInterval.parse(interval)
                val sequenceName = genomeInterval.getSequenceName
                if (!refSequenceSet.contains(sequenceName)) {
                    throw new RuntimeException(genomeInterval.getSequenceName + " is an invalid sequence")
                }
                map += sequenceName -> (map.getOrElse(sequenceName, List[String]()) ++ List(interval))
            })
        }
        map
    }

    lazy val getGenomeMaskFiles = {
        if (getReferenceMetadata.genomeMasks.size == 1) {
            getReferenceMetadata.genomeMasks :+ getReferenceMetadata.lcGenomeMask
        } else {
            getReferenceMetadata.genomeMasks
        }
    }

    var sentinelFileDirectory: File = null

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        if (runDirectory == null) {
            throw new ArgumentException("Argument with name '-runDirectory' is missing")
        }

        sentinelFileDirectory = new File(runDirectory, "cnv_sentinel_files")
        createDirectory(sentinelFileDirectory)

        mergeBamHeaders(bamHeadersMergedBam)

        createDirectory(genderGenotypeFilterDir)
        addCommand(new CreateGenderGenotypeFilterFile)

        createStages(1, math.min(4, lastStage))
        if (lastStage >= 5) {
            addCommand(new StageScript_5)
        }
        if (lastStage >= 6) {
            createStages(6, math.min(11, lastStage))
        }
        if (lastStage >= 12) {
            addCommand(new StageScript_12)
        }
    }

    def getStageDirectory(stageNumber: Int): File = {
        new File(runDirectory, "cnv_stage" + stageNumber)
    }

    def getSeqDirectory(stageNumber: Int, sequenceName: String): File = {
        val stageDirectory = getStageDirectory(stageNumber)
        new File(stageDirectory, "seq_" + sequenceName)
    }

    def sentinelFile(stageNumber: Int, sequenceName: String) = {
        val fileName = sequenceName match {
            case seq: Any => "stage_" + stageNumber + "_seq_" + seq + ".sent"
            case _        => "stage_" + stageNumber + ".sent"
        }
        new File(sentinelFileDirectory, fileName)
    }

    // The following functions define output files from the pipeline stages that are used be other stages

    // BAM creation
    def bamHeadersDirectory = {
        new File(runDirectory, "bam_headers")
    }

    def bamHeadersMergedBam = {
        new File(bamHeadersDirectory, "merged_headers.bam")
    }

    def bamHeadersMergedBamIndex = {
        new File(bamHeadersDirectory, "merged_headers.bam.bai")
    }

    def genderGenotypeFilterDir = {
        new File(runDirectory, "gender_gt_filters")
    }

    def genderGenotypeFilterFile = {
        new File(genderGenotypeFilterDir, "gender_gt_filter.txt")
    }

    def genderGenotypeFilterDescriptionFile = {
        new File(genderGenotypeFilterDir, "gender_gt_filter_descr.txt")
    }

    // Stage_1 annotated VCF file of the initial scanned windows
    def scannedWindowsVcfFile(sequenceName: String): File = {
        new File(getSeqDirectory(1, sequenceName), "seq_" + sequenceName + ".sites.vcf.gz")
    }

    // Stage_2 Genotyped VCF file of the initial scanned windows
    def genotypedWindowsVcfFile(sequenceName: String): File = {
        new File(getSeqDirectory(2, sequenceName), "seq_" + sequenceName + ".genotypes.vcf.gz")
    }

    // Stage_3 Merged scanned windows VCF file
    def mergedSitesVcfFile(sequenceName: String): File = {
        new File(getSeqDirectory(3, sequenceName), "seq_" + sequenceName + ".merged.sites.vcf.gz")
    }

    // Stage_4 Merged windows genotypes VCF file
    def mergedGenotypesVcfFile(sequenceName: String): File = {
        new File(getSeqDirectory(4, sequenceName), "seq_" + sequenceName + ".merged.genotypes.vcf.gz")
    }

    // Stage_5 Discovery samples list file
    def selectedSamplesList = {
        new File(getStageDirectory(5), "eval/DiscoverySamples.list")
    }

    // Stage_5 VPS report created by merging VPS reports for all the sequences
    def mergedVPSReportFile = {
        new File(getStageDirectory(5), "eval/VariantsPerSample.report.dat")
    }

    // Stage_6 Selected sites list file
    def selectedSiteListFile(sequenceName: String): File = {
        new File(getSeqDirectory(6, sequenceName), "eval/SelectedVariants.list")
    }

    // Stage_6 Based on a sub-list of samples, produce a new headers-only BAM file
    def selectedSamplesMergedHeadersBam(sequenceName: String): File = {
        new File(getSeqDirectory(6, sequenceName), "seq_" + sequenceName + ".merged_headers.bam")
    }

    // Stage_7 Brig merged site Vcf file
    def brigSiteVcfFile(sequenceName: String): File = {
        new File(getSeqDirectory(7, sequenceName), "seq_" + sequenceName + ".brig.sites.vcf.gz")
    }

    // Stage_8 Brig merged genotypes Vcf file
    def brigGenotypesVcfFile(sequenceName: String): File = {
        new File(getSeqDirectory(8, sequenceName), "seq_" + sequenceName + ".brig.genotypes.vcf.gz")
    }

    // Stage_9 Adjacent merged sites vcf file
    def adjacentMergedSitesVcfFile(sequenceName: String): File = {
        new File(getSeqDirectory(9, sequenceName), "seq_" + sequenceName + ".adjacent_merged.sites.vcf.gz")
    }

    // Stage_10 Adjacent merged sites genotypes vcf file
    def adjacentMergedGenotypesVcfFile(sequenceName: String): File = {
        new File(getSeqDirectory(10, sequenceName), "seq_" + sequenceName + ".adjacent_merged.genotypes.vcf.gz")
    }

    // Stage_11 Quality & length filtered vcf file
    def qualLengthFilteredVcfFile(sequenceName: String): File = {
        new File(getSeqDirectory(11, sequenceName), "seq_" + sequenceName + ".filtered.genotypes.vcf.gz")
    }

    // Stage_12 Merged vcf file
    def mergedVcfFile = {
        new File(getStageDirectory(12), outputVcfName)
    }

    // Creates a range of stages from firstStage to lastStage
    def createStages(firstStage: Int, lastStage: Int) = {
        sequenceIntervalMap.keys.foreach (sequenceName => {
            var stageSentinel: File = null

            (firstStage to lastStage).foreach (stageNumber => {
                val seqDirectory = getSeqDirectory(stageNumber, sequenceName)
                createDirectory(seqDirectory)

                // Class.forName("CNVDiscoveryPipeline#StageScript_" + stageNumber).newInstance.asInstanceOf[StageScript]
                var script: StageScript = null
                stageNumber match {
                    case 1 =>
                        script = new StageScript_1(sequenceName)
                    case 2 =>
                        script = new StageScript_2(sequenceName)
                    case 3 =>
                        script = new StageScript_3(sequenceName)
                    case 4 =>
                        script = new StageScript_4(sequenceName)
                    case 6 =>
                        script = new StageScript_6(sequenceName)
                    case 7 =>
                        script = new StageScript_7(sequenceName)
                    case 8 =>
                        script = new StageScript_8(sequenceName)
                    case 9 =>
                        script = new StageScript_9(sequenceName)
                    case 10 =>
                        script = new StageScript_10(sequenceName)
                    case 11 =>
                        script = new StageScript_11(sequenceName)
                    case _ =>
                        throw new RuntimeException("Stage " + stageNumber + " is not defined")
                }
                stageSentinel = addScript(script, stageSentinel)
            })
        })
    }

    def addScript(script: StageScript, previousStageSentinel: File): File = {
        if (previousStageSentinel != null) {
            script.dependsOnFile :+= previousStageSentinel
        }
        addCommand(script)
    }

    abstract class StageScript(sequenceName: String) extends QueueScriptCommand with UnifiedInputOutput {
        val stageNumber = getClass.getName.split("_")(1).toInt
        createDirectory(getStageDirectory(stageNumber))

        this.commandResultFile = sentinelFile(stageNumber, sequenceName)

        override def scriptName = "discovery/cnv/CNVDiscoveryStage" + stageNumber + ".q"
        includeScripts :+= "discovery/cnv/CNVDiscoveryStageBase.q"
        includeScripts :+= "discovery/cnv/CNVDiscoveryGenotyper.q"
        def workingDirectory = getSeqDirectory(stageNumber, sequenceName)
        override def jobLogDir = new File(workingDirectory, "logs")

        // BobH: Not sure if this memory limit is too high, but will try this for now.
        // With many windows, the Queue jobs were being killed for exceeding the LSF memory limit.
        this.lsfMemLimit = Some(12)

        // BobH: Not sure if these should go here, but I removed them from the base class QueueScriptCommand.
        // Not all Queue scripts will take all of these arguments, so the base class should only process arguments that we know are common to all Queue scripts.
        // It would be nice if we had a mechanism where we could introspect and say "if command/script X takes parameter Y, then pass it along...".
        // Note: We now pass -disableJobReport by default.  The job reports cause no end of problems (seems to be a performance/synchronization bottleneck, can create permissions problems, etc.).
        val commonStageArguments =
            flag(" --disableJobReport ", true) +
            optional(" -configFile ", configFile) +
            repeat(" -P ", parameterList) +
            optional(" -R ", referenceFile) +
            optional(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
            repeat(" -genomeMaskFile ", getGenomeMaskFiles) +
            optional(" -copyNumberMaskFile ", getReferenceMetadata.copyNumberMask) +
            optional(" -readDepthMaskFile ", getReferenceMetadata.readDepthMask) +
            optional(" -genderMaskBedFile ", getReferenceMetadata.genderMaskBed) +
            optional(" -vdjBedFile ", getReferenceMetadata.vdjRegionsBed) +
            repeat(" -genderMapFile ", getGenderMapFileList) +
            repeat(" -md ", metaDataLocationList) +
            repeat(" -excludeReadGroup ", excludeReadGroup) +
            repeat(" -genotypeFilterFile ", genotypeFilterFileList) +
            repeat(" -keepTempFiles ", keepTempFiles) +
            flag(" -disableGATKTraversal ", disableGATKTraversal) +
            optional(" -maxConcurrentRun ", maxConcurrentStageJobs)

        override def commandLine =
            super.commandLine +
            optional(" -sequenceName ", sequenceName) +
            required(" -runDirectory ", workingDirectory) +
            required(" -sentinelFile ", sentinelFile(stageNumber, sequenceName)) +
            commonStageArguments
    }

    class StageScript_1(sequenceName: String) extends StageScript(sequenceName) {
        this.dependsOnFile :+= genderGenotypeFilterFile

        override def commandLine =
            super.commandLine +
            required(" -I ", bamHeadersMergedBam) +
            repeat(" -intervalList ", sequenceIntervalMap(sequenceName)) +
            required(" -scannedWindowsVcfFile ", scannedWindowsVcfFile(sequenceName)) +
            required(" -tilingWindowSize ", tilingWindowSize) +
            required(" -tilingWindowOverlap ", tilingWindowOverlap) +
            optional(" -maximumReferenceGapLength ", maximumReferenceGapLength) +
            optional(" -testMaximumWindowCount ", testMaximumWindowCount)
    }

    class StageScript_2(sequenceName: String) extends StageScript(sequenceName) {
        override def commandLine =
            super.commandLine +
            required(" -I ", bamHeadersMergedBam) +
            required(" -vcf ", scannedWindowsVcfFile(sequenceName)) +
            required(" -genderGenotypeFilterFile ", genderGenotypeFilterFile) +
            required(" -filterDescriptionFile ", genderGenotypeFilterDescriptionFile) +
            optional(" -genotypingParallelRecords ", genotypingParallelRecords)
    }

    class StageScript_3(sequenceName: String) extends StageScript(sequenceName) {
        override def commandLine =
            super.commandLine +
            required(" -I ", bamHeadersMergedBam) +
            required(" -vcf ", genotypedWindowsVcfFile(sequenceName)) +
            required(" -mergedVcfFile ", mergedSitesVcfFile(sequenceName)) +
            required(" -duplicateScoreThresholdMax ", 0)
    }

    class StageScript_4(sequenceName: String) extends StageScript(sequenceName) {
        override def commandLine =
            super.commandLine +
            required(" -I ", bamHeadersMergedBam) +
            required(" -vcf ", mergedSitesVcfFile(sequenceName)) +
            required(" -genderGenotypeFilterFile ", genderGenotypeFilterFile) +
            required(" -filterDescriptionFile ", genderGenotypeFilterDescriptionFile) +
            optional(" -genotypingParallelRecords ", genotypingParallelRecords)
    }

    // Merging of results from ALL the sequences in Stage_4
    class StageScript_5 extends StageScript(null) {
        this.commandResultFile = selectedSamplesList
        sequenceIntervalMap.keys.foreach (sequenceName => {
            this.dependsOnFile :+= sentinelFile(4, sequenceName)
        })

        override def workingDirectory = getStageDirectory(5)
        override def commandLine =
            super.commandLine +
            required(" -I ", bamHeadersMergedBam) +
            required(" -vpsReportsDirectory ", getStageDirectory(4)) +
            required(" -selectedSamplesList ", selectedSamplesList) +
            flag(" -useAllSamplesInDiscovery ", useAllSamplesInDiscovery)
    }

    class StageScript_6(sequenceName: String) extends StageScript(sequenceName) {
        this.dependsOnFile :+= selectedSamplesList

        // See if there is a sequence-specific SelectedSamplesList. If not, use the generic one
        var seqSelectedSamplesList = new File(getStageDirectory(5), "eval/seq_" + sequenceName + "_DiscoverySamples.list")
         if (!seqSelectedSamplesList.exists) {
            seqSelectedSamplesList = selectedSamplesList
         }

        override def commandLine =
            super.commandLine +
            required(" -I ", bamHeadersMergedBam) +
            required(" -vcf ", mergedGenotypesVcfFile(sequenceName)) +
            required(" -selectedSamplesList ", seqSelectedSamplesList) +
            required(" -selectedSamplesMergedHeadersBam ", selectedSamplesMergedHeadersBam(sequenceName))
    }

    class StageScript_7(sequenceName: String) extends StageScript(sequenceName) {
        override def commandLine =
            super.commandLine +
            required(" -I ", selectedSamplesMergedHeadersBam(sequenceName)) +
            required(" -vcf ", mergedGenotypesVcfFile(sequenceName)) +
            required(" -siteListFile ", selectedSiteListFile(sequenceName)) +
            required(" -boundaryPrecision ", boundaryPrecision) +
            optional(" -maximumReferenceGapLength ", maximumReferenceGapLength) +
            required(" -minimumRefinedLength ", minimumRefinedLength) +
            required(" -brigSiteVcfFile ", brigSiteVcfFile(sequenceName)) +
            flag(" -produceAuxiliaryFiles ", produceAuxiliaryFiles)
    }

    class StageScript_8(sequenceName: String) extends StageScript(sequenceName) {
        override def commandLine =
            super.commandLine +
            required(" -I ", selectedSamplesMergedHeadersBam(sequenceName)) +
            required(" -vcf ", brigSiteVcfFile(sequenceName)) +
            required(" -genderGenotypeFilterFile ", genderGenotypeFilterFile) +
            required(" -filterDescriptionFile ", genderGenotypeFilterDescriptionFile) +
            optional(" -genotypingParallelRecords ", genotypingParallelRecords) +
            optional(" -duplicateScoreThreshold ", duplicateScoreThreshold)
    }

    class StageScript_9(sequenceName: String) extends StageScript(sequenceName) {
        override def commandLine =
            super.commandLine +
            required(" -I ", selectedSamplesMergedHeadersBam(sequenceName)) +
            required(" -vcf ", brigGenotypesVcfFile(sequenceName)) +
            required(" -adjacentMergedSitesVcf ", adjacentMergedSitesVcfFile(sequenceName))
    }

    class StageScript_10(sequenceName: String) extends StageScript(sequenceName) {
        override def commandLine =
            super.commandLine +
            required(" -I ", bamHeadersMergedBam) +
            required(" -vcf ", adjacentMergedSitesVcfFile(sequenceName)) +
            required(" -genderGenotypeFilterFile ", genderGenotypeFilterFile) +
            required(" -filterDescriptionFile ", genderGenotypeFilterDescriptionFile) +
            optional(" -genotypingParallelRecords ", genotypingParallelRecords)
    }

    class StageScript_11(sequenceName: String) extends StageScript(sequenceName) {
        override def commandLine =
            super.commandLine +
            required(" -I ", selectedSamplesMergedHeadersBam(sequenceName)) +
            required(" -vcf ", adjacentMergedGenotypesVcfFile(sequenceName)) +
            required(" -filteredVcf ", qualLengthFilteredVcfFile(sequenceName))
    }

    class StageScript_12 extends StageScript(null) {
        var vcfFileList: List[File] = Nil
        sequenceIntervalMap.keys.foreach (sequenceName => {
            this.dependsOnFile :+= sentinelFile(11, sequenceName)
            vcfFileList +:= qualLengthFilteredVcfFile(sequenceName)
        })
        override def workingDirectory = getStageDirectory(12)
        override def commandLine =
            super.commandLine +
            repeat(" -filteredVcf ", vcfFileList) +
            required(" -mergedVcfFile ", mergedVcfFile) +
            required(" -outputDirectory ", runDirectory) +
            required(" -genderGenotypeFilterFile ", genderGenotypeFilterFile) +
            required(" -filterDescriptionFile ", genderGenotypeFilterDescriptionFile)
    }

    class CreateGenderGenotypeFilterFile extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.CreateGenderGenotypeFilterFile"
        this.dependsOnFile :+= bamHeadersMergedBam
        this.outputFile = genderGenotypeFilterFile
        commandArguments +=
            required(" -R ", referenceFile) +
            required(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
            repeat(" -genderMapFile ", getGenderMapFileList) +
            required(" -filterDescriptionFile ", genderGenotypeFilterDescriptionFile)
    }
}
