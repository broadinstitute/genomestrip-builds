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
package org.broadinstitute.sv.qscript

import java.io.File
import java.io.FileOutputStream
import java.net.URL
import collection.JavaConverters._
import collection.Map
import scala.io.Source.fromFile
import org.broadinstitute.sv.queue.ComputeDiscoveryPartitions
import org.broadinstitute.sv.queue.ComputeGenomeLocusPartitions
import org.broadinstitute.sv.queue.ComputeVCFPartitions
import org.broadinstitute.sv.util.GenomeInterval
import org.broadinstitute.sv.util.fasta.IndexedFastaFile
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.function.QFunction
import org.broadinstitute.gatk.utils.jna.lsf.v7_0_6.LibBat
import org.broadinstitute.gatk.utils.jna.lsf.v7_0_6.LibBat.submit
import org.broadinstitute.gatk.utils.jna.lsf.v7_0_6.LibLsf
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature

//@QscriptDocStart
/**
 * Queue utilities for writing structural variation discovery/genotyping pipepines.
 */
@DocumentedGATKFeature(groupName = "Queue Scripts")
abstract class SVQScript extends QScript {

    @Input(shortName="gatk", doc="The path to the GenomeAnalysisTK.jar file")
    var gatkJar: File = null

    @Argument(shortName="cp", required=false, doc="The java classpath")
    var classpath: String = null

    @Argument(shortName="variantType", required=false, doc="Variant discovery")
    var variantType: String = null

    @Input(shortName="R", doc="The reference file for the bam files")
    var referenceFile: File = null

    @Input(shortName="altMaskFile", required=false, doc="The genome mask file(s) for the alternate alleles")
    var altMaskFiles: List[File] = Nil

    @Input(shortName="readDepthMaskFile", required=false, doc="Genome mask file specifying regions over which the sequencing depth should be estimated")
    var readDepthMaskFile: File = null

    @Argument(shortName="I", required=false, doc="One or more bam files or a bam file list (with extension .list)")
    var bamInputs: List[String] = Nil

    @Input(fullName="altAlleleAlignments", shortName="altAlignments", required=false, doc="Alternate allele alignments")
    var altAlignmentFiles: List[File] = Nil

    @Argument(fullName="debug", shortName="debug", required=false, doc="Enable verbose output for debugging")
    var debug: String = null

    @Argument(fullName="verbose", shortName="verbose", required=false, doc="Enable verbose output")
    var verbose: String = null

    @Argument(fullName="metaDataLocation", shortName="md", required=false, doc="Path to meta-data directory")
    var metaDataLocationList: List[String] = Nil

    @Argument(fullName="referenceMetaDataLocation", shortName="rmd", required=false, doc="Path to the directory containing the reference companinon files, such as ploidy map, gcmask file, etc.")
    var referenceMetaDataLocation: File = null

    @Argument(shortName="configFile", required=false, doc="Configuration property file")
    var configFile: File = null

    @Argument(fullName="parameter", shortName="P", required=false, doc="Override individual parameters from the configuration property file")
    var parameterList: List[String] = Nil

    @Argument(shortName="depthMinimumMappingQuality", required=false, doc="Minimum mapping quality for read depth analysis")
    var depthMinimumMappingQuality: java.lang.Integer = 10

    @Argument(shortName="depthMaximumInsertSizeRadius", required=false, doc="Maximum insert size radius for read depth analysis")
    var depthMaximumInsertSizeRadius: java.lang.Double = 10.0

    @Argument(shortName="excludeReadGroup", required=false, doc="Read group(s) to exclude or .list files of read groups")
    var excludeReadGroup: List[String] = null

    @Argument(shortName="sample", required=false, doc="Sample(s) to process or .list files of samples")
    var sampleList: List[String] = Nil

    @Argument(shortName="minimumInsertSizeMapFile", required=false, doc="File containing a manual overrides of minimum insert size used to select read pairs")
    var minimumInsertSizeMapFile: File = null

    @Argument(shortName="bamFilesAreDisjoint", required=false, doc="Reduces memory footprint if each input BAM file contains separate samples/libraries/read groups")
    var bamFilesAreDisjointArg: String = null
    var bamFilesAreDisjoint: Boolean = false

    @Argument(shortName="altAlleleFlankLength", required=false, doc="Length of flanking sequence to use for alt allele mapping")
    var altAlleleFlankLength: java.lang.Integer = null

    @Argument(shortName="alignUnmappedMates", required=false, doc="Whether to align unmapped mates to alternate alleles (default true)")
    var alignUnmappedMates: Boolean = true

    @Argument(fullName="jobLogDirectory", shortName="jobLogDir", required=false, doc="Directory for queue output files")
    var jobLogDirectory: File = null

    @Argument(shortName="runDirectory", required=false, doc="Directory for auxilliary output files")
    var runDirectory: File = null

    @Argument(fullName="discoveryWindowSize", shortName="windowSize", required=false, doc="Window size for SV discovery")
    var discoveryWindowSize: java.lang.Integer = null;

    @Argument(fullName="discoveryWindowPadding", shortName="windowPadding", required=false, doc="Window padding for SV discovery")
    var discoveryWindowPadding: java.lang.Integer = null;

    @Argument(fullName="discoveryMinimumSize", shortName="minimumSize", required=false, doc="Minimum event size for SV discovery")
    var discoveryMinimumSize: java.lang.Integer = null;

    @Argument(fullName="discoveryMaximumSize", shortName="maximumSize", required=false, doc="Maximum event size for SV discovery")
    var discoveryMaximumSize: java.lang.Integer = null;

    @Argument(fullName="genomeInterval", shortName="L", required=false, doc="Interval(s) (or file of intervals) used in preprocessing and discovery")
    var genomeIntervalList: List[String] = null;

    @Argument(shortName="parallelJobs", required=false, doc="Number of parallel jobs to run")
    var parallelJobs: java.lang.Integer = null;

    @Argument(shortName="parallelRecords", required=false, doc="Number of records to process in each parallel job")
    var parallelRecords: java.lang.Integer = null;

    @Argument(shortName="lsfResource", required=false, doc="Resource request rusage string to supply to LSF")
    var lsfResource: String = null;

    @Argument(shortName="lsfMemLimit", required=false, doc="Override default LSF MEMLIMIT setting (in GB)")
    var lsfDefaultMemLimit: java.lang.Integer = null;

    @Argument(shortName="vcfCachingDistance", required=false, doc="Caching distance for VCF records to preserve sort order (default 10K)")
    var vcfCachingDistance: java.lang.Integer = null;

    @Argument(shortName="exome", required=false, doc="True to enable exome genotyping mode (experimental) (default is false)")
    var exomeMode: String = null;

    @Argument(shortName="exomeDepthMatrixFile", required=false, doc="Input data matrix for exome genotyping (experimental)")
    var exomeDepthMatrixFile: String = null;

    @Argument(shortName="duplicateOverlapThreshold", required=false, doc="Overlap threshold for duplicate site filtering (default 0.5)")
    var duplicateOverlapThreshold: java.lang.Double = 0.5;

    @Argument(shortName="duplicateScoreThreshold", required=false, doc="Concordance LOD score threshold (minimum per-sample LOD score at which events are considered non-redundant, default zero)")
    var duplicateScoreThreshold: java.lang.Double = 0.0;

    @Argument(shortName="nonVariantScoreThreshold", required=false, doc="Score threshold for non-variant site filtering (phred scale, default 13.0)")
    var nonVariantScoreThreshold: java.lang.Double = 13.0;

    @Argument(shortName="suppressVCFCommandLines", required=false, doc="Flag indicating whether to suppress the command line header in the output VCF file")
    var suppressVCFCommandLines: Boolean = false;

    @Argument(shortName="disableGATKTraversal", required=false, doc="Disable GATK normal walker traversal")
    var disableGATKTraversal: Boolean = true;

    @Argument(fullName="inputFileIndexCache", shortName="IC", required=false, doc="Location where the input BAM files' indices are cached")
    var inputFileIndexCache: String = null;

    @Argument(shortName="gatkJobRunner", required=false, doc="Job runner to be used when submitting scripts from within this script (Lsf706/GridEngine/Drmaa)")
    var gatkJobRunner: String = "Lsf706";

    @Argument(shortName="jobWrapperScript", required=false, doc="Wrapper script to be invoked by Queue. Can be used to setup necessary environment variables required by the wrapped job to run")
    var jobWrapperScript: String = null;

    @Argument(shortName="genderMapFile", required=false, doc="File(s) containing the declared gender for each sample")
    private var genderMapFile: List[File] = Nil

    // These are reference-dependent metadata files.
    // If the value is not supplied, it will be defaulted from the reference location
    // Therefore these should be accessed through ReferenceMetadata object
    @Input(shortName="genomeMaskFile", required=false, doc="The genome mask file(s) for the reference file")
    private var genomeMaskFiles: List[File] = Nil

    @Input(shortName="copyNumberMaskFile", required=false, doc="Genome mask file for regions of potentially polymorphic copy number")
    private var copyNumberMaskFile: File = null

    @Input(shortName="genderMaskBedFile", required=false, doc="Bed file defining the chrX and chrY regions that should be used in gender calling")
    private var genderMaskBedFile: File = null

    @Input(shortName="vdjBedFile", required=false, doc="Bed file defining the vdj regions")
    private var vdjBedFile: File = null

    @Input(shortName="ploidyMapFile", required=false, doc="Map file defining the gender-specific ploidy of each region of the reference sequence")
    private var ploidyMapFile: File = null

    @Input(shortName="profileIntervalBedFile", required=false, doc="Bed file defining the default intervals for depth profiling")
    private var profileIntervalBedFile: File = null


//@QscriptDocEnd

    var bamLocationList: List[String] = Nil

    var referenceMetadata: ReferenceMetadata = null

    def metaDataLocation = {
        if (metaDataLocationList.size > 1) {
            throw new IllegalArgumentException("Only one metadataLocation must be specified")
        }
        metaDataLocationList(0)
    }

    // List conataining locations of all BAM files provided as input
    def bamLocations = {
        if (bamLocationList == Nil) {
            bamLocationList = bamInputs flatMap parseBAMsInput
        }
        bamLocationList
    }

    def getReferenceMetadata = {
        if (referenceMetadata == null) {
            referenceMetadata = new ReferenceMetadata(referenceFile, referenceMetaDataLocation, "")
        }
        referenceMetadata
    }

    def addCommand(function: QFunction with UnifiedInputOutput) = {
        super.add(function)
        if (function.commandResultFile != null) function.commandResultFile else function.outputFile
    }

    def required(prefix: String, param: Any, suffix: String = "") = {
        "%s%s%s".format(prefix, param, suffix)
    }

    def flag(prefix: String, param: Boolean) = {
        if (param) prefix else ""
    }

    def createDirectory(file: File) = {
        file.mkdirs()
    }

    def parseBoolean(value: String): Boolean = {
        value.toLowerCase() match {
            case "true" => true
            case "false" => false
            case _ => throw new IllegalArgumentException("Invalid boolean value: " + value);
        }
    }

    def parseBAMsInput(bamsIn: String): List[String] = bamsIn match {
        case bamPath if bamPath.endsWith("bam")  => return List(bamsIn)
        case bamPath if bamPath.endsWith("cram")  => return List(bamsIn)
        case bamPath if bamPath.endsWith("list") => return (for (line <- fromFile(new File(bamsIn)).getLines) yield line).toList
        case _ => throw new RuntimeException("Unexpected BAM input type: " + bamsIn + "; only permitted extensions are .bam .cram and .list")
    }

    def parameterFiles: List[File] = {
        List(configFile, getReferenceMetadata.parameterFile).filter(_ != null)
    }

    def expandListFiles(args: List[String]) : List[String] = {
        return args flatMap expandListFile
    }

    def expandListFile(arg: String) : List[String] = {
        if (arg.endsWith(".list") && new File(arg).exists) {
            return fromFile(new File(arg)).getLines.toList
        } else {
            return List(arg)
        }
    }

    def expandFileListFiles(args: List[File]) : List[File] = {
        return expandListFiles(args.map { file => file.getPath() }) map { path => new File(path) }
    }

    def computeReferenceSequenceNames() : List[String] = {
        new IndexedFastaFile(referenceFile).getSequenceNames().asScala.toList
    }

    def computeVCFPartitions(file: File): Map[String, String]  = {
        computeVCFPartitions(file, parallelJobs, parallelRecords)
    }

    def computeVCFPartitions(file: File, parallelJobs: java.lang.Integer, parallelRecords: java.lang.Integer): Map[String, String] = {
        val alg = new ComputeVCFPartitions()
        alg.setParallelJobs(parallelJobs)
        alg.setParallelRecords(parallelRecords)
        alg.computePartitions(file).asScala
    }

    def computeLocusPartitions(locusSize: java.lang.Integer) = {
        val intervalList = if (genomeIntervalList == null) null else genomeIntervalList.asJava
        val alg = new ComputeGenomeLocusPartitions()
        alg.setReferenceFile(referenceFile)
        alg.setLocusSize(locusSize)
        alg.setGenomeIntervalList(intervalList)
        alg.computePartitions().asScala
    }

    def computeDiscoveryPartitions() = {
        val alg = new ComputeDiscoveryPartitions()
        alg.setReferenceFile(referenceFile)
        alg.setSearchWindowSize(discoveryWindowSize)
        alg.setSearchWindowPadding(discoveryWindowPadding)
        alg.setSearchMinimumSize(discoveryMinimumSize)
        alg.setSearchMaximumSize(discoveryMaximumSize)
        alg.setSearchInterval(discoverySearchLocus)
        alg.computePartitions().asScala
    }

    def computeDiscoverySearchWindow(locus: String, padding: java.lang.Integer) = {
        if (locus == null || padding == null) {
            locus
        } else {
            var interval = GenomeInterval.parse(locus)
            var start = if (interval.getStart() == 0) 0 else math.max(1, interval.getStart() - padding.intValue)
            var end = if (interval.getEnd() == 0) 0 else interval.getEnd() + padding.intValue
            new GenomeInterval(interval.getSequenceName(), start, end).toString()
        }
    }

    def discoverySearchLocus: String = {
        var searchLocus: String = null
        if (genomeIntervalList != null && genomeIntervalList.size > 0) {
            if (genomeIntervalList.size == 1) {
                searchLocus = genomeIntervalList(0)
            } else {
                throw new IllegalArgumentException("Discovery can be run on only one locus interval");
            }
        }
        searchLocus
    }

    trait UnifiedInputOutput {
        @Input(fullName="dependsOnFile", shortName="dependsOnFile", doc="Auxilliary dependencies, not passed as arguments", required=false)
        var dependsOnFile: List[File] = Nil
        @Output(fullName="outputFile", shortName="O", doc="Output file", required=false)
        var outputFile: File = _
        @Output(fullName="commandResultFile", shortName="commandResultFile", doc="Final command result for dependency tracking (default outputFile)", required=false)
        var commandResultFile: File = _
    }

    trait FileInputOutput extends UnifiedInputOutput {
        @Input(fullName="inputFile", shortName="I", doc="Input file", required=false)
        var inputFile: List[File] = Nil
    }

    trait BAMInputOutput extends UnifiedInputOutput {
        @Argument(fullName="inputLocation", shortName="inputLocation", doc="Input location, either a local file or a URL", required=false)
        var inputLocation: List[String] = Nil
    }

    // The Queue commandLine method causes the command line arguments to be bound very late, well after the command
    // object has been constructed.  This is useful sometimes, but more often it is useful to bind the command arguments
    // at construction time.  For example, to give a different run directory to a pipeline running one part of a
    // larger process.  We introduce a commandArguments string that is bound early, at command construction time,
    // to allow us to have the best of both worlds.
    // Unfortunately, we have to merge this with commandLine on every implementing class, as it cannot be done once in the trait
    // because the method on the base class is abstract.
    trait EarlyBoundCommandArgs extends QFunction {
        var commandArguments: String = ""
    }

    trait LogDirectory extends QFunction {
        override def freezeFieldValues = {
            super.freezeFieldValues
            this.jobOutputFile = new File(jobName + ".out")
            if (jobLogDirectory != null) {
                this.jobOutputFile = new File(jobLogDirectory, this.jobOutputFile.getPath())
            }
        }
        // Disable the new Queue behavior of tracking job output files as job output artifacts
        override def statusPaths = {
            outputs
        }
    }

    trait DefaultMemoryLimit extends CommandLineFunction {
        override def freezeFieldValues = {
            if (this.memoryLimit.isEmpty && this.qSettings.memoryLimit.isEmpty) {
                this.memoryLimit = Some(3)
            }
            super.freezeFieldValues
        }
    }

    trait LSFOptions extends CommandLineFunction {
        var lsfMemLimit : Option[Int] = None
        this.updateJobRun = {
            case lsfJob: submit => {
                if (lsfResource != null) {
                    if (lsfResource.indexOf('[') < 0) {
                        throw new IllegalArgumentException("Malformed LSF resource string: " + lsfResource)
                    }
                    lsfJob.resReq = mergeLsfResourceStrings(lsfJob.resReq, lsfResource)
                    lsfJob.options |= LibBat.SUB_RES_REQ
                }
                if (!this.lsfMemLimit.isEmpty) {
                    lsfJob.rLimits(LibLsf.LSF_RLIMIT_RSS) = this.lsfMemLimit * 1000
                } else if (lsfDefaultMemLimit != null) {
                    lsfJob.rLimits(LibLsf.LSF_RLIMIT_RSS) = lsfDefaultMemLimit * 1000
                }
            }
        }

        // This code is not really correct, but will work in simple cases.
        // Some known defects:
        // - The code is intolerant of spaces in the resource specifications.
        // - The rules for merging multiple resoruce types should be type specific,
        //   for example, select should use &&, not comma.
        def mergeLsfResourceStrings(res1: String, res2: String) : String = {
            if (res1 == null) {
                res2
            } else if (res2 == null) {
                res1
            } else {
                val map1: Map[String, List[String]] = parseLsfResource(res1)
                val map2: Map[String, List[String]] = parseLsfResource(res2)
                val mergedMap: Map[String, List[String]] = mergeLsfResources(map1, map2)
                formatLsfResource(mergedMap)
            }
        }

        def mergeLsfResources(map1: Map[String, List[String]], map2: Map[String, List[String]]) : Map[String, List[String]] = {
            val result = collection.mutable.Map(map1.toSeq: _*)
            for ((key, value2) <- map2) {
                result.get(key) match {
                    case Some(value1) => {
                        result.put(key, value1 ++ value2)
                    }
                    case None => {
                        result.put(key, value2)
                    }
                }
            }
            return result
        }

        def parseLsfResource(res: String) : Map[String, List[String]] = {
            val map = new collection.mutable.HashMap[String, List[String]]
            for (fieldx <- res.split("\\s+")) {
                val ind1 = fieldx.indexOf('[')
                val ind2 = fieldx.length-1
                if (ind1 > 0 && fieldx.charAt(ind2) == ']') {
                    val key = fieldx.substring(0,ind1)
                    val values : List[String] = fieldx.substring(ind1+1,ind2).split(",\\s*").toList
                    map.put(key, values)
                }
            }
            return map
        }

        def formatLsfResource(m: Map[String, List[String]]) : String = {
            return ((m map { case (key,list) => key + "[" + list.mkString(",") + "]" }).mkString(" "))
        }
    }

    trait DependsOnMetaData {
        @Input(fullName="mdDependency", shortName="mdDependency", doc="Metadata version dependency", required=false)
        var mdDependency: File = new File(metaDataLocation, "mdversion.txt")
    }

    trait DependsOnISDFile {
        @Input(fullName="isdDependency", shortName="isdDependency", doc="Insert size distribution dependency", required=false)
        var isdDependency: File = new File(metaDataLocation, "isd.stats.dat")
    }

    def baseName(path: String): String = {
        new File(path).getName()
    }

    // The default QScript version removes the path, which is not desirable.
    override def swapExt(file: File, oldExtension: String, newExtension: String) : File = {
        new File(file.getPath.stripSuffix("." + oldExtension) + "." + newExtension)
    }

    def stripExt(file: File, extension: String) : File = {
        new File(file.getPath.stripSuffix("." + extension))
    }

    def swapExt(fileName: String, oldExtension: String, newExtension: String) : String = {
        fileName.stripSuffix("." + oldExtension) + "." + newExtension
    }

    def getGenderMapFileList = {
        if (genderMapFile == Nil) {
            genderMapFile = metaDataLocationList
                .map(md => new File(md, "sample_gender.report.txt"))
                .filter(file => file.exists)
        }
        genderMapFile
    }

    def injectJobWrapperScript = {
        if (jobWrapperScript != null) jobWrapperScript + " " else ""
    }

    class ReferenceMetadata (referenceFile: File,
                             referenceMetaDataLocation: File,
                             basename: String) {

        object Extensions {
            val ReferenceFasta = "fasta"
            val GenomeMask = "svmask.fasta"
            val LCGenomeMask = "lcmask.fasta"
            val CopyNumberMask = "gcmask.fasta"
            val ReadDepthMask = "rdmask.bed"
            val GenderMaskBed = "gendermask.bed"
            val PloidyMap = "ploidymap.txt"
            val VDJRegionsBed = "vdjregions.bed"
            val ProfileIntervalBed = "profiles.bed"
            val GSParams = "gsparams.txt"
        }

        lazy val rmdLocation = {
            if (referenceMetaDataLocation == null) {
                new File(referenceFile.getParent)
            } else {
                referenceMetaDataLocation
            }
        }

        def getMetadataFile(file: File, extension: String) = {
            if (file == null) {
                val mdFile = new File(rmdLocation, referenceFile.getName.stripSuffix("." + Extensions.ReferenceFasta) + "." + extension)
                if (mdFile.exists) mdFile else null
            } else {
                file
            }
        }

        def getMetadataFileList(fileList: List[File], extension: String) = {
            if (fileList == Nil) {
                val file = getMetadataFile(null, extension)
                if (file != null) List(file) else Nil
            } else {
                fileList
            }
        }

        lazy val genomeMasks = getMetadataFileList(genomeMaskFiles, Extensions.GenomeMask)
        lazy val lcGenomeMask = getMetadataFile(null, Extensions.LCGenomeMask)
        lazy val copyNumberMask = getMetadataFile(copyNumberMaskFile, Extensions.CopyNumberMask)
        lazy val readDepthMask = getMetadataFile(readDepthMaskFile, Extensions.ReadDepthMask)
        lazy val genderMaskBed = getMetadataFile(genderMaskBedFile, Extensions.GenderMaskBed)
        lazy val ploidyMap = getMetadataFile(ploidyMapFile, Extensions.PloidyMap)
        lazy val vdjRegionsBed = getMetadataFile(vdjBedFile, Extensions.VDJRegionsBed)
        lazy val profileIntervalBed = getMetadataFile(profileIntervalBedFile, Extensions.ProfileIntervalBed)
        lazy val parameterFile = getMetadataFile(null, Extensions.GSParams)
    }

    class TouchSentinelFile(val fileToTouch: File, val dependsOn: List[File]) extends InProcessFunction
           with LogDirectory with UnifiedInputOutput {

        def this(fileToTouch: File, dependsOn: File) = this(fileToTouch, List(dependsOn))

        this.dependsOnFile ++= dependsOn
        this.outputFile = fileToTouch

        def run = {
            if (!fileToTouch.exists()) {
                new FileOutputStream(fileToTouch).close();
            }
            fileToTouch.setLastModified(System.currentTimeMillis);
        }
    }

    abstract class SimpleCommand extends org.broadinstitute.gatk.queue.function.CommandLineFunction
        with FileInputOutput with LogDirectory with DefaultMemoryLimit with LSFOptions with EarlyBoundCommandArgs {
        override def commandLine =
	    injectJobWrapperScript +
	    commandArguments
    }

    abstract class JavaCommand extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
        with FileInputOutput with LogDirectory with DefaultMemoryLimit with LSFOptions with EarlyBoundCommandArgs {
        override def javaOpts = super.javaOpts + optional(" -cp ", classpath)
        override def commandLine =
            injectJobWrapperScript +
            super.commandLine +
            repeat(" -I ", inputFile) +
            optional(" -O ", outputFile) +
            commandArguments
    }

    abstract class PicardCommand extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
        with FileInputOutput with LogDirectory with DefaultMemoryLimit with LSFOptions with EarlyBoundCommandArgs {
        override def javaOpts = super.javaOpts + optional(" -cp ", classpath)
        override def commandLine =
            super.commandLine +
            repeat(" I=", inputFile) +
            optional(" O=", outputFile) +
            commandArguments
    }

    abstract class GATKWalkerCommand extends org.broadinstitute.gatk.queue.extensions.gatk.CommandLineGATK
        with FileInputOutput with LogDirectory with DefaultMemoryLimit with LSFOptions with EarlyBoundCommandArgs {
        this.jarFile = gatkJar
        override def commandLine =
	    injectJobWrapperScript +
	    super.commandLine +
            commandArguments
    }

    abstract class SVWalkerCommand extends org.broadinstitute.gatk.queue.extensions.gatk.CommandLineGATK
        with BAMInputOutput with LogDirectory with DefaultMemoryLimit with LSFOptions with EarlyBoundCommandArgs {
        override def javaExecutable = "org.broadinstitute.sv.main.SVCommandLine"
        this.jarFile = gatkJar
        this.reference_sequence = referenceFile
        override def javaOpts = super.javaOpts + optional(" -cp ", classpath)
        override def commandLine =
            injectJobWrapperScript +
            super.commandLine +
            repeat(" -I ", inputLocation) +
            optional(" -O ", outputFile) +
            optional(" -IC ", inputFileIndexCache) +
            optional(" -disableGATKTraversal ", disableGATKTraversal) +
            commandArguments
        commandArguments += repeat(" -md ", metaDataLocationList)
    }

    abstract class QueueScriptCommand(submitToFarm: Boolean = true) extends JavaCommand {
        def scriptName: String
        def jobLogDir: File = jobLogDirectory
        var includeScripts: List[String] = Nil
        this.javaMainClass = "org.broadinstitute.gatk.queue.QCommandLine"
        this.commandResultFile = new File(runDirectory, baseName(scriptName) + ".out")

        val svHome = System.getenv().get("SV_DIR")
        val qscriptHome = new File(svHome, "qscript")

        def submissionArgs =
            if (submitToFarm) {
                required(" -jobRunner ", gatkJobRunner) +
                required(" -gatkJobRunner ", gatkJobRunner) +
                optional(" -resMemReqParam ", qSettings.residentRequestParameter) +
                repeat(" -jobNative ", qSettings.jobNativeArgs) +
                optional(" -jobQueue ", jobQueue) +
                optional(" -jobWrapperScript ", jobWrapperScript)
            } else {
                ""
            }

        def resourceArgs =
            optional(" -memLimit ", qSettings.memoryLimit)

        override def commandLine =
            super.commandLine +
            optional(" -cp ", classpath) +
            required(" -S ", new File(qscriptHome, scriptName)) +
            repeat(" -S ", includeScripts map (script => new File(qscriptHome, script))) +
            required(" -S ", new File(qscriptHome, "SVQScript.q")) +
            required(" -gatk ", gatkJar) +
//            " --disableJobReport " +
            optional( " -jobLogDir ", jobLogDir) +
            resourceArgs +
            submissionArgs +
            " -run "
    }

    // Utility Commands

    class WriteFileList(inputFiles: List[File], outFile: File) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.queue.WriteFileList"
        this.inputFile = inputFiles
        this.outputFile = outFile
    }

    class CopyFile(inFile: File, outFile: File) extends SimpleCommand {
        commandArguments = "cp %s %s".format(inFile, outFile)
        this.inputFile :+= inFile
        this.outputFile = outFile
    }


    // SV Commands

    class CreateMetaDataDirectory() extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.CreateMetaDataDirectory"
        val mdVersionFile = new File(metaDataLocation, "mdversion.txt")
        this.commandResultFile = mdVersionFile
        commandArguments +=
            required(" -md ", metaDataLocation) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList)
    }

    class ComputeGenomeSizes extends JavaCommand with DependsOnMetaData {
        this.javaMainClass = "org.broadinstitute.sv.apps.ComputeGenomeSizes"
        val genomeSizesFile = new File(metaDataLocation, "genome_sizes.txt")
        this.outputFile = genomeSizesFile
        commandArguments +=
            required(" -R ", referenceFile) +
            required(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -genomeInterval ", genomeIntervalList)
    }

    def isdDir = new File(metaDataLocation, "isd")
    def depthDir = new File(metaDataLocation + "/depth")
    def spanDir = new File(metaDataLocation + "/spans")
    def gcProfDir = new File(metaDataLocation + "/gcprofile")
    def rcCacheDir = new File(metaDataLocation + "/rccache")
    def rcCacheMergeDir(cacheName: String) = new File(metaDataLocation, cacheName + ".merge")

    class ComputeMetadata(bamLocation: String, computeGCProfiles: Boolean, computeReadCounts: Boolean, gcRefProfile: File) extends SVWalkerCommand with DependsOnISDFile {
        @Output(fullName="depthFile", shortName="depthFile", doc="Output file", required=false)
        var depthFile: File = _
        @Output(fullName="spanFile", shortName="spanFile", doc="Output file", required=false)
        var spanFile: File = _
        @Output(fullName="gcProfileFile", shortName="gcProfileFile", doc="Output file", required=false)
        var gcProfileFile: File = _
        @Output(fullName="readCountFile", shortName="readCountFile", doc="Output file", required=false)
        var readCountFile: File = _

        this.analysis_type = "ComputeMetadataWalker"
        this.inputLocation :+= bamLocation

        // Read-depth
        createDirectory(depthDir)
        depthFile = new File(depthDir, swapExt(baseName(bamLocation), "bam", "depth.txt"))

        // Spans
        createDirectory(spanDir)
        spanFile = new File(spanDir, swapExt(baseName(bamLocation), "bam", "spans.txt"))

        // GC Profile
        if (computeGCProfiles) {
            createDirectory(gcProfDir)
            gcProfileFile = new File(gcProfDir, baseName(bamLocation) + ".gcprof.zip")
        }

        // Read counts cache
        if (computeReadCounts) {
            createDirectory(rcCacheDir)
            readCountFile = new File(rcCacheDir, swapExt(baseName(bamLocation), "bam", "rc.bin"))
        }

        commandArguments +=
                required(" -computeGCProfiles", String.valueOf(computeGCProfiles)) +
                required(" -computeReadCounts", String.valueOf(computeReadCounts)) +
                required(" -depthFile ", depthFile) +
                required(" -spanFile ", spanFile) +
                optional(" -gcProfileFile ", gcProfileFile) +
                optional(" -readCountFile ", readCountFile) +
                repeat(" -configFile ", parameterFiles) +
                repeat(" -P ", parameterList) +
                required(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
                repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
                optional(" -copyNumberMaskFile ", getReferenceMetadata.copyNumberMask) +
                optional(" -readDepthMaskFile ", getReferenceMetadata.readDepthMask) +
                repeat(" -genomeInterval ", genomeIntervalList) +
                optional(" -minMapQ ", depthMinimumMappingQuality) +
                required(" -maxInsertSizeStandardDeviations ", 3) +
                optional(" -insertSizeRadius ", depthMaximumInsertSizeRadius) +
                repeat(" -excludeReadGroup ", excludeReadGroup) +
                required(" -referenceProfile ", gcRefProfile)
    }

    class ComputeInsertSizeHistograms(bamLocation: String, createHistogramFile: Boolean = true) extends SVWalkerCommand with DependsOnMetaData {
        this.analysis_type = "ComputeInsertSizeHistogramsWalker"
        this.inputLocation :+= bamLocation
        createDirectory(isdDir)

        def histFile = new File(isdDir, swapExt(baseName(bamLocation), "bam", "hist.bin"))
        def chimerismFile = new File(isdDir, swapExt(baseName(bamLocation), "bam", "chimer.dat"))

        this.outputFile = if (createHistogramFile) histFile else chimerismFile;

        commandArguments +=
            repeat(" -genomeInterval ", genomeIntervalList) +
            repeat(" -excludeReadGroup ", excludeReadGroup) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            required(" -chimerismFile ", chimerismFile) +
            optional(" -createHistogramFile ", createHistogramFile) +
            " -createEmpty "
    }

    class MergeInsertSizeHistograms(inputFiles: List[File], disjoint: Boolean = bamFilesAreDisjoint) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeInsertSizeHistograms"
        this.inputFile = inputFiles
        this.outputFile = new File(metaDataLocation, "isd.hist.bin")
        commandArguments += flag(" -disjoint", disjoint)
    }

    class ReduceInsertSizeHistograms(histFile: File) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.ReduceInsertSizeHistograms"
        this.inputFile :+= histFile
        this.outputFile = swapExt(histFile, "hist.bin", "dist.bin")
    }

    class MergeInsertSizeDistributions(inputFiles: List[File]) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeInsertSizeDistributions"
        this.inputFile = inputFiles
        this.outputFile = new File(metaDataLocation, "isd.dist.bin")
    }

    class ComputeInsertStatistics(inFile: File) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.ComputeInsertStatistics"
        this.inputFile :+= inFile
        this.outputFile = new File(metaDataLocation, "isd.stats.dat")
    }

    class ComputeReadDepthCoverage(bamLocation: String) extends SVWalkerCommand with DependsOnISDFile {
        this.analysis_type = "ComputeReadDepthCoverageWalker"
        this.inputLocation :+= bamLocation
        createDirectory(depthDir)
        this.outputFile = new File(depthDir, swapExt(baseName(bamLocation), "bam", "depth.txt"))
        commandArguments +=
            required(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -genomeInterval ", genomeIntervalList) +
            optional(" -readDepthMaskFile ", getReferenceMetadata.readDepthMask) +
            optional(" -minMapQ ", depthMinimumMappingQuality) +
            optional(" -insertSizeRadius ", depthMaximumInsertSizeRadius) +
            repeat(" -excludeReadGroup ", excludeReadGroup)
    }

    class MergeReadDepthCoverage(inputFiles: List[File]) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeReadDepthCoverage"
        this.inputFile = inputFiles
        this.outputFile = new File(metaDataLocation, "depth.dat")
    }

    class ComputeReadCounts(bamLocation: String) extends SVWalkerCommand with DependsOnISDFile {
        this.analysis_type = "ComputeReadCountsWalker"
        this.inputLocation :+= bamLocation
        createDirectory(rcCacheDir)
        this.outputFile = new File(rcCacheDir, swapExt(baseName(bamLocation), "bam", "rc.bin"))
        commandArguments +=
            repeat(" -genomeInterval ", genomeIntervalList) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            optional(" -insertSizeRadius ", depthMaximumInsertSizeRadius) +
            repeat(" -excludeReadGroup ", excludeReadGroup)
    }

    class MergeReadCountsBase(inputFiles: List[File], outFile: File, partitionArgs: Array[String] = null) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeReadCounts"
        this.inputFile = inputFiles
        this.outputFile = outFile
        // This can be removed if we implement indexing on the fly, but for now we need to generate a dependency on the index files.
        this.dependsOnFile = expandFileListFiles(inputFiles).map { rcFile => swapExt(rcFile, "bin", "bin.idx") }
        var partitionArgString = ""
        if (partitionArgs != null) {
            partitionArgs.foreach(arg => partitionArgString += (" " + arg))
        }
        commandArguments +=
            required(" -R ", referenceFile) +
            partitionArgString
    }

    class MergeReadCounts(inputFiles: List[File], outFile: File = new File(metaDataLocation, "rccache.bin")) extends MergeReadCountsBase(inputFiles, outFile) {
        // commandArguments += required(" -debug ", "true");
    }

    class WriteMergeReadCountsInputList(inputFiles: List[File], outFile: File) extends WriteFileList(inputFiles, outFile) {
        // Necessary to track dependencies on rc index files as well as rc data files
        this.dependsOnFile = inputFiles.map { rcFile => swapExt(rcFile, "bin", "bin.idx") }
    }

    def mergeReadCountsByLocus(inputFiles: List[File], cacheName: String = "rccache") : File = {
        val rcListFile = new File(metaDataLocation, cacheName + ".list")
        addCommand(new WriteMergeReadCountsInputList(inputFiles, rcListFile))
        createDirectory(rcCacheMergeDir(cacheName))
        var mergeInputList: List[File] = Nil
        mergeInputList :+= rcListFile
        var mergedPartitionFiles: List[File] = Nil
        val partitions = computeLocusPartitions(10000000)
        for ((partitionName, partitionArgs) <- partitions) {
            val mergedFile = new File(rcCacheMergeDir(cacheName), partitionName + ".rccache.bin")
            val mergeCountsCmd = new MergeReadCountsBase(mergeInputList, mergedFile, partitionArgs)
            mergeCountsCmd.jobArrayName = "mergeReadCounts"
            mergedPartitionFiles :+= addCommand(mergeCountsCmd)
            val indexCountsCmd = new IndexReadCountFile(mergedFile)
            indexCountsCmd.jobArrayName = "indexReadCounts"
            addCommand(indexCountsCmd)
        }
        val mergedCountFile = new File(metaDataLocation, cacheName + ".bin")
        addCommand(new MergeReadCounts(mergedPartitionFiles, mergedCountFile))
        val mergedCountFileIdx = addCommand(new IndexReadCountFile(mergedCountFile))
        return mergedCountFileIdx
    }

    class IndexReadCountFile(rcFile: File) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.IndexReadCountFile"
        this.inputFile :+= rcFile
        this.outputFile = swapExt(rcFile, "bin", "bin.idx")
    }

    class ComputeReadSpanCoverage(bamLocation: String) extends SVWalkerCommand with DependsOnISDFile {
        this.analysis_type = "ComputeReadSpanCoverageWalker"
        this.inputLocation :+= bamLocation
        createDirectory(spanDir)
        this.outputFile = new File(spanDir, swapExt(baseName(bamLocation), "bam", "spans.txt"))
        commandArguments +=
            repeat(" -genomeInterval ", genomeIntervalList) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            required(" -maxInsertSizeStandardDeviations ", 3) +
            repeat(" -excludeReadGroup ", excludeReadGroup)
    }

    class MergeReadSpanCoverage(inputFiles: List[File]) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeReadSpanCoverage"
        this.inputFile = inputFiles
        this.outputFile = new File(metaDataLocation, "spans.dat")
    }

    class ComputeGCReferenceProfile extends JavaCommand with DependsOnMetaData {
        this.javaMainClass = "org.broadinstitute.sv.apps.ComputeGCProfiles"
        createDirectory(gcProfDir)
        this.outputFile = new File(gcProfDir, "reference.gcprof.zip")
        commandArguments +=
            required(" -R ", referenceFile) +
            required(" -md ", metaDataLocation) +
            required(" -writeReferenceProfile ", "true") +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            optional(" -copyNumberMaskFile ", getReferenceMetadata.copyNumberMask) +
            repeat(" -genomeInterval ", genomeIntervalList) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList)
    }

    class ComputeGCProfiles(bamLocation: String, refProfile: File) extends SVWalkerCommand with DependsOnISDFile {
        this.analysis_type = "ComputeGCProfileWalker"
        this.dependsOnFile :+= refProfile
        this.inputLocation :+= bamLocation
        this.outputFile = new File(gcProfDir, baseName(bamLocation) + ".gcprof.zip")
        commandArguments +=
            required(" -referenceProfile ", refProfile) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -genomeInterval ", genomeIntervalList) +
            optional(" -copyNumberMaskFile ", getReferenceMetadata.copyNumberMask) +
            optional(" -insertSizeRadius ", depthMaximumInsertSizeRadius) +
            repeat(" -excludeReadGroup ", excludeReadGroup)
    }

    class MergeGCProfiles(inputFiles: List[File]) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeGCProfiles"
        this.inputFile = inputFiles
        this.outputFile = new File(metaDataLocation + "/gcprofiles.zip")
    }

    class MergeTextOutput(inputFiles: List[File], outFile: File, headerLines: java.lang.Integer = null) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeTextOutput"
        this.inputFile = inputFiles
        this.outputFile = outFile
        optional(" -headerLines ", headerLines)
    }

    class GenerateAltAlleleFasta(vcfFile: File) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.GenerateAltAlleleFasta"
        val altAlignDir = new File(runDirectory, "altalign")
        createDirectory(altAlignDir)
        val fastaFile = new File(altAlignDir, swapExt(vcfFile, "vcf", "alt.fasta").getName())
        this.inputFile :+= vcfFile
        this.outputFile = fastaFile
        commandArguments +=
            required(" -R ", referenceFile) +
            optional(" -flankLength ", altAlleleFlankLength)
    }

    class GenerateAltAlleleDictionary(altReference: File) extends PicardCommand {
        this.javaMainClass = "picard.sam.CreateSequenceDictionary"
        val dictFile = swapExt(altReference, "fasta", "dict")
        this.dependsOnFile :+= altReference
        this.outputFile = dictFile
        commandArguments +=
            required(" R=", altReference) +
            " TRUNCATE_NAMES_AT_WHITESPACE=TRUE"
    }

    class IndexFastaFile(fastaFile: File) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.IndexFastaFile"
        val indexFile = new File(fastaFile.getPath + ".fai")
        this.inputFile :+= fastaFile
        this.outputFile = indexFile
    }

    class CreateBwaIndex(fastaFile: File) extends SimpleCommand {
        val indexFile = new File(fastaFile.getPath() + ".bwt")
        commandArguments = "bwa index %s".format(fastaFile)
        this.inputFile :+= fastaFile
        this.outputFile = indexFile
    }

    class CreateBamIndex(bamFile: File) extends SimpleCommand {
        val indexFile = new File(bamFile.getPath() + ".bai")
        commandArguments = "samtools index %s".format(bamFile)
        this.inputFile :+= bamFile
        this.outputFile = indexFile
    }

    class MergeSamFiles(inFiles: List[File], outFile: File) extends PicardCommand {
        this.javaMainClass = "picard.sam.MergeSamFiles"
        this.inputFile = inFiles
        this.outputFile = outFile
        this.memoryLimit = Some(8)
        commandArguments += " VALIDATION_STRINGENCY=SILENT"
    }

    class SVAltAligner(bamLocation: String, altReference: File, outFile: File = null) extends SVWalkerCommand {
        override def javaExecutable = "org.broadinstitute.gatk.engine.CommandLineGATK"
        this.analysis_type = "SVAltAlignerWalker"
        this.inputLocation :+= bamLocation
        this.outputFile = outFile
        if (outFile == null) {
            val altAlignDir = new File(runDirectory, "altalign")
            createDirectory(altAlignDir)
            this.outputFile = new File(altAlignDir, swapExt(baseName(bamLocation), "bam", "alt.bam"))
        }
        val dictFile = swapExt(altReference, "fasta", "dict")
        val bwaIndexFile = new File(altReference.getPath() + ".bwt")
        this.dependsOnFile :+= dictFile
        this.dependsOnFile :+= bwaIndexFile
        // The SVAltAligner memory requirements are unusual.
        // We use about 4g of non-heap memory and a relatively small java heap.
        this.memoryLimit = Some(5)
        this.javaMemoryLimit = Some(1)
        commandArguments +=
            required(" -U ", "ALLOW_UNINDEXED_BAM") +
            required(" -altReference ", altReference) +
            optional(" -alignUnmappedMates ", alignUnmappedMates)
    }

    class SVDiscoveryBase(outFile: File) extends SVWalkerCommand {
        override def javaExecutable = "org.broadinstitute.sv.main.SVDiscovery"
        this.analysis_type = "SVDiscoveryWalker"
        this.inputLocation = bamInputs
        this.outputFile = outFile
        this.memoryLimit = Some(6)
        this.javaMemoryLimit = Some(3)
        this.lsfMemLimit = Some(12)
        override def javaOpts = super.javaOpts
        commandArguments +=
            optional(" -debug ", debug) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            optional(" -runDirectory ", runDirectory) +
            repeat(" -genderMapFile ", getGenderMapFileList) +
            repeat(" -excludeReadGroup ", excludeReadGroup) +
            optional(" -variantType", variantType) +
            optional(" -minimumInsertSizeMapFile ", minimumInsertSizeMapFile) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            optional(" -vcfCachingDistance ", vcfCachingDistance)
    }

    class SVDiscovery(outFile: File, runFilePrefix: String = null, storeReadPairFile: Boolean = true) extends SVDiscoveryBase(outFile) {
        val searchLocus = discoverySearchLocus
        val searchWindow = computeDiscoverySearchWindow(searchLocus, discoveryWindowPadding)
        commandArguments +=
            optional(" -L ", searchWindow) +
            optional(" -runFilePrefix ", runFilePrefix) +
            optional(" -searchLocus ", searchLocus) +
            optional(" -searchWindow ", searchWindow) +
            optional(" -searchMinimumSize ", discoveryMinimumSize) +
            optional(" -searchMaximumSize ", discoveryMaximumSize) +
            optional(" -storeReadPairFile ", storeReadPairFile)
    }

    class SVParallelDiscovery(partitionName: String, partitionArgs: Array[String], storeReadPairFile: Boolean = true)
        extends SVDiscoveryBase(new File(runDirectory, partitionName + ".discovery.vcf.gz")) {
        var partitionArgString = ""
        partitionArgs.foreach(arg => partitionArgString += (" " + arg))
        commandArguments +=
            required(" -partitionName ", partitionName) +
            required(" -runFilePrefix ", partitionName) +
            optional(" -storeReadPairFile ", storeReadPairFile) +
            partitionArgString
    }

    class MergeDiscoveryOutput(outFile: File, partitionFiles: List[File]) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeDiscoveryOutput"
        this.outputFile = outFile
        if (outFile == null) {
            this.outputFile = new File(runDirectory, "genstrip.discovery.vcf.gz")
        }
        this.dependsOnFile ++= partitionFiles
        commandArguments +=
            required(" -R ", referenceFile) +
            optional(" -runDirectory ", runDirectory)
    }

    class SVVariantFiltration(vcfFile: File, outFile: File) extends GATKWalkerCommand {
        this.analysis_type = "VariantFiltration"
        this.inputFile :+= vcfFile
        this.outputFile = outFile
        commandArguments +=
            flag(" -no_cmdline_in_header ", suppressVCFCommandLines) +
            required(" -V ", vcfFile) +
            required(" -o ", outFile) +
            required(" -R ", referenceFile)
    }

    class SVDiscoveryDefaultFilter(vcfFile: File, outFile: File) extends SVVariantFiltration(vcfFile, outFile) {
        this.dependsOnFile :+= vcfFile
        commandArguments +=
            " -filterName COVERAGE -filter \"GSDEPTHCALLTHRESHOLD == \\\"NA\\\" || GSDEPTHCALLTHRESHOLD >= 1.0\"" +
            " -filterName COHERENCE -filter \"GSCOHPVALUE == \\\"NA\\\" || GSCOHPVALUE <= 0.01\"" +
            " -filterName DEPTHPVAL -filter \"GSDEPTHPVALUE == \\\"NA\\\" || GSDEPTHPVALUE >= 0.01\"" +
            " -filterName DEPTH -filter \"GSDEPTHRATIO == \\\"NA\\\" || GSDEPTHRATIO > 0.8 || (GSDEPTHRATIO > 0.63 && (GSMEMBPVALUE == \\\"NA\\\" || GSMEMBPVALUE >= 0.01))\"" +
            " -filterName PAIRSPERSAMPLE -filter \"GSNPAIRS <= 1.1 * GSNSAMPLES\""
    }

    class SVDiscoveryPassFilter(vcfFile: File, outFile: File) extends GATKWalkerCommand {
        this.analysis_type = "SelectVariants"
        this.dependsOnFile :+= vcfFile
        this.outputFile = outFile
        commandArguments +=
            required(" -R ", referenceFile) +
            required(" -V ", vcfFile) +
            required(" -o ", outputFile) +
            " -ef "
    }

    class SVGenotyper(vcfFile: File, outFile: File, runFilePrefix: String = null, genotypeFilterFileList: List[File] = Nil) extends SVWalkerCommand {
        override def javaExecutable = "org.broadinstitute.sv.main.SVGenotyper"
        this.analysis_type = "SVGenotyperWalker"
        this.inputLocation = bamInputs
        this.outputFile = outFile
        this.memoryLimit = Some(6)
        this.javaMemoryLimit = Some(3)

        // Note: By default, SVGenotyper no longer utilizese "standard" GATK traversal for scalability.
        // The standard traversal is still supported, however, and would use a VCF rod like this:
        //   required(" -BTI ", "input") +
        //   required(" -B:input,VCF ", vcfFile)

        commandArguments +=
            optional(" -debug ", debug) +
            optional(" -verbose ", verbose) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            optional(" -runDirectory ", runDirectory) +
            optional(" -runFilePrefix ", runFilePrefix) +
            repeat(" -genotypeFilterFile ", genotypeFilterFileList) +
            repeat(" -genderMapFile ", getGenderMapFileList) +
            optional(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -altMaskFile ", altMaskFiles) +
            repeat(" -altAlignments ", altAlignmentFiles) +
            required(" -vcf ", vcfFile) +
            repeat(" -sample ", sampleList) +
            repeat(" -excludeReadGroup ", excludeReadGroup) +
            optional(" -exome ", exomeMode) +
            optional(" -exomeDepthMatrixFile ", exomeDepthMatrixFile)
    }

    class SVParallelGenotyper(vcfFile: File, partitionName: String, partitionArg: String, genotypeFilterFileList: List[File] = Nil)
        extends SVGenotyper(vcfFile, new File(runDirectory, partitionName + ".genotypes.vcf.gz"), null, genotypeFilterFileList) {
        commandArguments +=
            required(" -partitionName ", partitionName) +
            required(" -partition ", partitionArg)
    }

    class MergeGenotyperOutput(vcfFile: File, outFile: File, partitionFiles: List[File], runFilePrefix: String = null) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeGenotyperOutput"
        this.outputFile = outFile
        if (outFile == null) {
            var outFileName = vcfFile.getName.stripSuffix(".gz").stripSuffix(".discovery.vcf").stripSuffix(".sites.vcf").stripSuffix(".vcf") + ".genotypes.vcf.gz"
            this.outputFile = new File(runDirectory, outFileName)
        }
        this.dependsOnFile ++= partitionFiles
        commandArguments +=
            required(" -R ", referenceFile) +
            optional(" -runDirectory ", runDirectory) +
            optional(" -runFilePrefix ", runFilePrefix)
    }

    class SVGenotyperDefaultFilter(vcfFile: File, outFile: File) extends SVVariantFiltration(vcfFile, outFile) {
        commandArguments +=
            " -filterName ALIGNLENGTH -filter \"GSELENGTH < 200\"" +
            " -filterName CLUSTERSEP -filter \"GSCLUSTERSEP == \\\"NA\\\" || GSCLUSTERSEP <= 2.0\"" +
            " -filterName GTDEPTH -filter \"GSM1 == \\\"NA\\\" || GSM1 <= 0.5 || GSM1 >= 2.0\"" +
            " -filterName NONVARIANT -filter \"GSNONVARSCORE != \\\"NA\\\" && GSNONVARSCORE >= " + nonVariantScoreThreshold + "\"" +
            " -filterName DUPLICATE -filter \"GSDUPLICATESCORE != \\\"NA\\\" && GSDUPLICATEOVERLAP >= " + duplicateOverlapThreshold + " && GSDUPLICATESCORE >= " + duplicateScoreThreshold + "\"" +
            " -filterName INBREEDINGCOEFF -filter \"GLINBREEDINGCOEFF != \\\"NA\\\" && GLINBREEDINGCOEFF < -0.15\""
    }

    class SVAnnotator(vcfFile: File, outFile: File) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.main.SVAnnotator"
        this.outputFile = outFile
        this.dependsOnFile :+= vcfFile
        val evalDir = new File(runDirectory, "eval")
        createDirectory(evalDir)
        val auxFilePrefix = vcfFile.getPath.stripSuffix(".gz").stripSuffix(".vcf").stripSuffix(".unfiltered")
        commandArguments +=
            required(" -vcf ", vcfFile) +
            required(" -R ", referenceFile) +
            repeat(" -md ", metaDataLocationList) +
            required(" -auxFilePrefix ", auxFilePrefix) +
            required(" -reportDirectory ", evalDir) +
            repeat(" -genderMapFile ", getGenderMapFileList) +
            optional(" -ploidyMapFile ", getReferenceMetadata.ploidyMap)
    }

    val genotyperDefaultAnnotatorList = List(
        "ClusterSeparation",
        "GCContent",
        "GenotypeLikelihoodStats",
        "NonVariant",
        "Redundancy",
        "CopyNumberClass",
        "CNQuality"
    )
    class GenotyperDefaultFilterAnnotations(vcfFile: File, outFile: File, addAnnotators: List[String] = Nil, skipAnnotators: List[String] = Nil)
        extends SVAnnotator(vcfFile, outFile) {
        val annotatorList = List.concat(genotyperDefaultAnnotatorList, addAnnotators).distinct diff skipAnnotators
        commandArguments +=
            annotatorList.map(e => "-A " + e).mkString(" ") +
            required(" -filterVariants ", "false") +
            required(" -writeReport ", "true") +
            required(" -writeSummary ", "true") +
            required(" -comparisonFile ", vcfFile) +
            required(" -duplicateOverlapThreshold ", duplicateOverlapThreshold) +
            required(" -duplicateScoreThreshold ", duplicateScoreThreshold)

        if (getReferenceMetadata.vdjRegionsBed != null && getReferenceMetadata.vdjRegionsBed.exists) {
            commandArguments +=
                required(" -A ", "BedFileOverlap") +
                required(" -bedFileTag ", "VDJ") +
                required(" -bedFile ", getReferenceMetadata.vdjRegionsBed)
        }
    }

    class GenotyperDefaultQCAnnotations(vcfFile: File, outFile: File)
        extends SVAnnotator(vcfFile, outFile) {
        val siteFiltersReportFile = new File(evalDir, "GenotypeSiteFilters.report.dat")
        val siteFiltersSummaryFile = new File(evalDir, "GenotypeSiteFilters.summary.dat")
        this.commandResultFile = siteFiltersSummaryFile
        commandArguments +=
            required(" -A ", "AlleleFrequency") +
            required(" -A ", "SiteFilters") +
            required(" -A ", "VariantsPerSample") +
            required(" -reportFileMap ", "SiteFilters:" + siteFiltersReportFile.getPath) +
            required(" -summaryFileMap ", "SiteFilters:" + siteFiltersSummaryFile.getPath) +
            required(" -writeReport ", "true") +
            required(" -writeSummary ", "true")
    }

    class IntensityRankSumAnnotator(vcfFile: File, samplesFileName: String, arrayIntensityFile: File)
        extends SVAnnotator(vcfFile, null) {
        val reportFile = new File(evalDir, "IntensityRankSum.report.dat")
        this.commandResultFile = reportFile
        commandArguments +=
            required(" -A ", "IntensityRankSum") +
            required(" -arrayIntensityFile ", arrayIntensityFile) +
            required(" -sample", samplesFileName) +
            required(" -irsSampleTag ", "GSSAMPLES") +
            required(" -writeReport ", "true") +
            required(" -writeSummary ", "true")
    }

    class CopyNumberClassAnnotator(vcfFile: File)
        extends SVAnnotator(vcfFile, null) {
        val reportFile = new File(evalDir, "CopyNumberClass.report.dat")
        this.commandResultFile = reportFile
        commandArguments +=
            required(" -A ", "CopyNumberClass") +
            required(" -writeReport ", "true") +
            required(" -writeSummary ", "true")
    }

    class BedFileAnnotator(vcfFile: File, outFile: File, bedFileTag: String, bedFile: String)
        extends SVAnnotator(vcfFile, outFile) {
        val reportFile = new File(evalDir, "BedFileOverlap.report.dat")
        this.commandResultFile = reportFile
        commandArguments +=
            required(" -A ", "BedFileOverlap") +
            required(" -bedFileTag ", bedFileTag) +
            required(" -bedFile ", bedFile) +
            required(" -writeReport ", "true") +
            required(" -writeSummary ", "false")
    }

    class CreateMergedBamHeaders(bamFile: File) 
        extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.ExtractBAMSubset"
        this.inputFile = bamInputs
        this.outputFile = bamFile
        this.memoryLimit = Some(6)
        this.javaMemoryLimit = Some(3)
        val readGroupMappingFiles =
            metaDataLocationList
                .map(md => new File(md, "read_groups.dat"))
                .filter(file => file.exists)
        commandArguments +=
            optional(" -R ", referenceFile) +
            optional(" -IC ", inputFileIndexCache) +
            repeat(" -sample ", sampleList) +
            required(" -L ", "NONE") +
            repeat(" -readGroupMappingFile ", readGroupMappingFiles)
    }

    def mergeBamHeaders(mergedBamFile: File) = {
        addCommand(new CreateMergedBamHeaders(mergedBamFile))
        addCommand(new CreateBamIndex(mergedBamFile))
    }
}
