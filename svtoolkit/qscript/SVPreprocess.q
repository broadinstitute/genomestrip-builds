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
package org.broadinstitute.sv.preprocess

import scala.io.Source
import scala.collection.JavaConverters._
import org.broadinstitute.sv.qscript.SVQScript
import org.broadinstitute.sv.util.GenomeInterval
import org.broadinstitute.sv.util.bed.BedFileLine
import org.broadinstitute.sv.util.bed.BedFileReader
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature

//@QscriptDocStart
@DocumentedGATKFeature(groupName = "Queue Scripts")
class SVPreprocess extends SVQScript {

    @Argument(shortName="useMultiStep", required=false, doc="Use more parallel jobs to compute the metadata")
    var useMultiStep: Boolean = false

    @Argument(shortName="reduceInsertSizeDistributions", required=false, doc="Reduces memory footprint by creating reduced representations of insert size distributions")
    var reduceInsertSizeDistributionsArg: String = null
    var reduceInsertSizeDistributions: Boolean = true

    @Argument(shortName="computeGCProfiles", required=false, doc="Compute GC bias profiles needed to do GC-bias normalization (requires copy number mask)")
    var computeGCProfilesArg: String = null
    var computeGCProfiles: Boolean = true

    @Argument(shortName="computeReadCounts", required=false, doc="Pre-compute read counts and create on-disk cache (experimental, not for general use yet)")
    var computeReadCountsArg: String = null
    var computeReadCounts: Boolean = true

    @Argument(shortName="deleteIntermediateDirs", required=false, doc="Delete the intermediate directories created during preprocessing")
    var deleteIntermediateDirectoriesArg: String = null
    var deleteIntermediateDirectories: Boolean = true

    @Argument(shortName="profileBinSize", required=false, doc="Size of profile bins to use for depth profiles")
    var profileBinSize: java.lang.Integer = 100000;

    @Argument(shortName="maximumReferenceGapLength", required=false, doc="Max reference gap in depth profiles")
    var maximumReferenceGapLength: java.lang.Integer = 10000

//@QscriptDocEnd

    var isdHistFiles: List[File] = Nil
    var depthFiles: List[File] = Nil
    var countFiles: List[File] = Nil
    var spanFiles: List[File] = Nil
    var gcProfileFiles: List[File] = Nil
    var gcRefProfile: File = null

    var mergeOutputFiles: List[File] = Nil

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        if (bamFilesAreDisjointArg != null) {
            bamFilesAreDisjoint = parseBoolean(bamFilesAreDisjointArg)
        }
        if (reduceInsertSizeDistributionsArg != null) {
            reduceInsertSizeDistributions = parseBoolean(reduceInsertSizeDistributionsArg)
        }
        if (computeGCProfilesArg != null) {
            computeGCProfiles = parseBoolean(computeGCProfilesArg)
        }
        if (computeReadCountsArg != null) {
            computeReadCounts = parseBoolean(computeReadCountsArg)
        }
        if (deleteIntermediateDirectoriesArg != null) {
            deleteIntermediateDirectories = parseBoolean(deleteIntermediateDirectoriesArg)
        }
        val mdVersionFile = addCommand(new CreateMetaDataDirectory())
        addCommand(new StampSVToolkitVesion(mdVersionFile))

        mergeBamHeaders(headersBam)

        addCommand(new ComputeGenomeSizes)
        if (computeGCProfiles) {
            gcRefProfile = addCommand(new ComputeGCReferenceProfile)
        }

        var chimerismFiles: List[File] = Nil
        for (bamLocation <- bamLocations) {
            val computeHistorgamCmd = new ComputeInsertSizeHistograms(bamLocation)
            isdHistFiles :+= addCommand(computeHistorgamCmd)
            chimerismFiles :+= computeHistorgamCmd.chimerismFile
        }
        if (reduceInsertSizeDistributions) {
            var isdDistFile:File = null
            if (bamFilesAreDisjoint) {
                var isdDistFiles: List[File] = Nil
                for (histFile <- isdHistFiles) {
                    isdDistFiles :+= addCommand(new ReduceInsertSizeHistograms(histFile))
                }
                isdDistFile = addCommand(new MergeInsertSizeDistributions(isdDistFiles))
            } else {
                val isdHistFile = addCommand(new MergeInsertSizeHistograms(isdHistFiles))
                isdDistFile = addCommand(new ReduceInsertSizeHistograms(isdHistFile))
            }
            val isdStatsFile = addCommand(new ComputeInsertStatistics(isdDistFile))
        } else {
            val isdHistFile = addCommand(new MergeInsertSizeHistograms(isdHistFiles))
            val isdStatsFile = addCommand(new ComputeInsertStatistics(isdHistFile))
        }
        mergeOutputFiles :+= addCommand(new MergeChimerismFiles(chimerismFiles))
        
        if (useMultiStep) {
            createMultiStepScript
        } else {
            createTwoStepScript
        }

        // Delete intermediate directories only after all the other tasks have finished
        val leafOutputFiles = createProfilesCallGenders()

        if (deleteIntermediateDirectories) {
            addCommand(new DeleteIntermediateDirectories(leafOutputFiles))
        }
    }

    def createMultiStepScript = {
        for (bamLocation <- bamLocations) {
            depthFiles :+= addCommand(new ComputeReadDepthCoverage(bamLocation))
        }
        mergeOutputFiles :+= addCommand(new MergeReadDepthCoverage(depthFiles))

        for (bamLocation <- bamLocations) {
            spanFiles :+= addCommand(new ComputeReadSpanCoverage(bamLocation))
        }
        mergeOutputFiles :+= addCommand(new MergeReadSpanCoverage(spanFiles))

        if (computeGCProfiles) {
            for (bamLocation <- bamLocations) {
                gcProfileFiles :+= addCommand(new ComputeGCProfiles(bamLocation, gcRefProfile))
            }
            mergeOutputFiles :+= addCommand(new MergeGCProfiles(gcProfileFiles))
        }

        if (computeReadCounts) {
            for (bamLocation <- bamLocations) {
                val countFile = addCommand(new ComputeReadCounts(bamLocation))
                addCommand(new IndexReadCountFile(countFile))
                countFiles :+= countFile
            }
            mergeOutputFiles :+= mergeReadCounts(countFiles)
        }
    }

    def createTwoStepScript = {
        for (bamLocation <- bamLocations) {
            val computeMetadataCommand = new ComputeMetadata(bamLocation, computeGCProfiles, computeReadCounts, gcRefProfile)
            addCommand(computeMetadataCommand)

            depthFiles :+= computeMetadataCommand.depthFile
            spanFiles :+= computeMetadataCommand.spanFile

            if (computeGCProfiles) {
                computeMetadataCommand.dependsOnFile :+= gcRefProfile
                gcProfileFiles :+= computeMetadataCommand.gcProfileFile
            }
            if (computeReadCounts) {
                countFiles :+= computeMetadataCommand.readCountFile
                addCommand(new IndexReadCountFile(computeMetadataCommand.readCountFile))
            }
        }

        mergeOutputFiles :+= addCommand(new MergeReadDepthCoverage(depthFiles))
        mergeOutputFiles :+= addCommand(new MergeReadSpanCoverage(spanFiles))

        if (computeGCProfiles) {
            mergeOutputFiles :+= addCommand(new MergeGCProfiles(gcProfileFiles))
        }
        if (computeReadCounts) {
            mergeOutputFiles :+= mergeReadCounts(countFiles)
        }
    }

    def headersBam : File = {
        new File(metaDataLocation, "headers.bam")
    }

    def headersBamIndex : File = {
        new File(metaDataLocation, "headers.bam.bai")
    }

    def mergeReadCounts(inputFiles: List[File], cacheName: String = "rccache") : File = {
        // Skip doing a by-locus merge if we are only processing a single input file.
        if (inputFiles.isEmpty) {
            return null
        } else if (inputFiles.length == 1) {
            val mergedCountFile = new File(metaDataLocation, cacheName + ".bin")
            addCommand(new CopyFile(inputFiles(0), mergedCountFile))
            val mergedCountFileIdx = addCommand(new IndexReadCountFile(mergedCountFile))
            return mergedCountFileIdx
        } else {
            mergeReadCountsByLocus(inputFiles, cacheName)
        }
    }

    def createProfilesCallGenders(): List[File] = {
        var profilesDirName: String = "profiles_"
        if (profileBinSize < 1000) {
            profilesDirName += profileBinSize
        } else {
            profilesDirName += (profileBinSize / 1000) + "Kb"
        }
        val profilesDir = new File(metaDataLocation, profilesDirName)
        createDirectory(profilesDir)
        var profileIndexFiles: List[File] = Nil
        for ((sequenceName, intervalList) <- profileIntervalMap) {
            var profileFile = addCommand(new ComputeDepthProfile(profilesDir, sequenceName, intervalList))
            profileIndexFiles :+= addCommand(new IndexDepthProfile(profileFile))
        }

        var leafOutputFiles: List[File] = Nil
        leafOutputFiles :+= addCommand(new ComputeSampleReadDepth(profilesDir, profileIndexFiles))

        val sampleGenderReport = addCommand(new CallSampleGender)
        val plotPdf = new File(metaDataLocation, "chrY_vs_chrX.pdf")
        leafOutputFiles :+= addCommand(new PlotChrYvsChrX(sampleGenderReport, plotPdf))
        return leafOutputFiles
    }

    def profileIntervalMap = {
        // Currently, we use one of three methods, in priority order.
        // If -L was supplied, generate one profile per chromosome (for every one mentioned in -L).
        // If we have reference metadata interval list (as a bed file), use the name field to build the interval map.
        // Otherwise, fall back to the historical method of one profile per reference sequence.
        if (genomeIntervalList != null) {
            val intervalList = expandListFiles(genomeIntervalList).map { interval => GenomeInterval.parse(interval) }
                                   .map { genomeInterval => genomeInterval.getSequenceName() -> genomeInterval }
            groupIntervalsByKeyPreservingOrder(intervalList)
        } else if (getReferenceMetadata.profileIntervalBed != null) {
            val intervalList = (new BedFileReader(getReferenceMetadata.profileIntervalBed) : java.lang.Iterable[BedFileLine]).asScala.toList
                                   .map { line => line.getField(3) -> line.getInterval }
            groupIntervalsByKeyPreservingOrder(intervalList)
        } else {
            val intervalList = computeReferenceSequenceNames.map { seq => seq -> GenomeInterval.parse(seq) }
            groupIntervalsByKeyPreservingOrder(intervalList)
        }
    }

    def groupIntervalsByKeyPreservingOrder(inputList: List[(String, GenomeInterval)]) = {
        val intervalMap = inputList.groupBy(_._1).map { case (k,v) => (k, v.map(_._2)) }
        val keyList = inputList.map(_._1)
        scala.collection.mutable.LinkedHashMap(keyList.map { k => k -> intervalMap(k) }: _*)
    }

    class StampSVToolkitVesion(mdVersionFile: File) extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
            with FileInputOutput {
        override def javaOpts = super.javaOpts + optional(" -cp ", classpath)
        this.javaMainClass = "org.broadinstitute.sv.apps.SVToolkitInfo"

        val svVersionFile = new File(metaDataLocation, "svtoolkit.version.dat")

        this.dependsOnFile :+= mdVersionFile
        this.outputFile = svVersionFile
        this.jobOutputFile = svVersionFile
    }

    class CallSampleGender() extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.CallSampleGender"
        this.dependsOnFile = headersBam :: headersBamIndex :: mergeOutputFiles
        this.outputFile = new File(metaDataLocation, "sample_gender.report.txt")

        commandArguments +=
            required(" -I ", headersBam) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            required(" -R ", referenceFile) +
            required(" -md ", metaDataLocation) +
            required(" -genderBedFile ", getReferenceMetadata.genderMaskBed) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -L ", genomeIntervalList) +
            optional(" -ploidyMapFile ", getReferenceMetadata.ploidyMap)
    }

    class MergeChimerismFiles(chimerismFiles: List[File]) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeChimerismFiles"
        this.dependsOnFile ++= isdHistFiles
        this.inputFile = chimerismFiles
        this.outputFile = new File(metaDataLocation, "chimerism.dat")
    }

    class ComputeDepthProfile(profilesDir: File, sequenceName: String, intervalList: List[GenomeInterval]) extends JavaCommand with BAMInputOutput {
        this.javaMainClass = "org.broadinstitute.sv.apps.ComputeDepthProfiles"
        this.dependsOnFile = headersBam :: headersBamIndex :: mergeOutputFiles
        this.outputFile = new File(profilesDir, "profile_seq_%s_%d.dat.gz".format(sequenceName, profileBinSize))

        commandArguments +=
            required(" -I ", headersBam) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            required(" -R ", referenceFile) +
            repeat(" -L ", intervalList) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -md ", metaDataLocationList) +
            required(" -profileBinSize ", profileBinSize) +
            optional(" -maximumReferenceGapLength ", maximumReferenceGapLength)
    }

    // Generate a tabix index for the profile file
    class IndexDepthProfile(profileFile: File)  extends SimpleCommand {
        val indexFile = new File(profileFile.getPath() + ".tbi")
        this.inputFile :+= profileFile
        this.outputFile = indexFile

        commandArguments = "tabix -f -S 1 -s 2 -b 3 -e 4 %s".format(profileFile)
    }

    class ComputeSampleReadDepth(profilesDir: File, profileIndexFiles: List[File]) extends JavaCommand with BAMInputOutput {
        this.javaMainClass = "org.broadinstitute.sv.apps.ComputeSampleReadDepth"
        this.dependsOnFile = profileIndexFiles
        this.outputFile = new File(profilesDir, "rd.dat")

        commandArguments +=
            required(" -profileDirectory ", profilesDir) +
            required(" -profileBinSize ", profileBinSize) +
            repeat(" -sequence ", profileIntervalMap.keySet)
    }

    class PlotChrYvsChrX(sampleGenderReport: File, plotPdf: File) extends SimpleCommand {
        this.dependsOnFile :+= sampleGenderReport
        this.commandResultFile = plotPdf

        val plotReadDepthScriptName = "metadata/plot_chr_vs_chr_readdepth.R"
        val plotReadDepthScriptPath = org.broadinstitute.sv.util.RUtilities.findRScript(plotReadDepthScriptName)
        commandArguments =
            required("Rscript") +
            required(plotReadDepthScriptPath) +
            required(sampleGenderReport) +
            required(plotPdf) +
            required("seq_Y vs. seq_X Read Depth") +
            required("DOSAGE_X") +
            required("DOSAGE_Y")
    }

    class DeleteIntermediateDirectories(leafOutputFiles: List[File]) extends SimpleCommand {
        this.dependsOnFile = leafOutputFiles
        this.outputFile = new File(metaDataLocation, "delete_intermediate_dirs.out")

        val intermediateDirectories = List(isdDir, depthDir, spanDir, gcProfDir, rcCacheDir, rcCacheMergeDir("rccache"))
        val builder = StringBuilder.newBuilder
        intermediateDirectories.foreach { dir =>
          builder
            .append("rm -rf ")
            .append(dir.getPath)
            .append("; ")
        }
        commandArguments = builder.toString
    }
}
