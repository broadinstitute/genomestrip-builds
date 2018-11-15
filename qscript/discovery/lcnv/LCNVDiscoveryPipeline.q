/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2016 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */

import java.io.File
import java.io.FileInputStream
import java.util.zip.GZIPInputStream
import scala.io.Source

import collection.JavaConverters._
import scala.collection.JavaConversions._

import org.broadinstitute.sv.commandline.CommandLineParser
import org.broadinstitute.sv.metadata.gender.{Gender, GenderMap}
import org.broadinstitute.sv.util.RUtilities
import org.broadinstitute.sv.util.bed.BedFileLine
import org.broadinstitute.sv.util.bed.BedFileReader

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.sv.qscript.SVQScript


class LCNVDiscoveryPipeline extends SVQScript {

    @Argument(shortName="profilesDir", required=true, doc="Profiles directory(ies)")
    var profilesDirList: List[File] = Nil

    @Argument(shortName="maxDepth", required=true, doc="Max depth of search")
    var maxDepth: java.lang.Integer = null

    @Argument(shortName="backgroundSample", required=false, doc="Samples to be used as background for the score computation")
    var backgroundSampleList: List[String] = Nil

    @Argument(shortName="perSampleScan", required=false, doc="If true, runs one scan job for each chrom/sample pair, otherwise runs one job for each chromosome (which scans all the samples internally)")
    var perSampleScanArg: String = null
    var perSampleScan: Boolean = true

    val RSCRIPT_NAME = "lcnv/lcnv_scan.R"
    var rScriptPath: String = null

    var binSize: Int = _

    def script = {
        rScriptPath = RUtilities.findRScript(RSCRIPT_NAME)

        validateProfilesDirectories()

        parseBinSize()
        if (maxDepth == null) {
            maxDepth = if (binSize >= 100000) 0 else 50
        }

        if (perSampleScanArg != null) {
            perSampleScan = parseBoolean(perSampleScanArg)
        }

        if (sampleList == Nil) {
            createSampleList()
        }
        println("sampleList.size: " + sampleList.size)

        var sentinelFiles: List[File] = Nil
        sequenceList.foreach(chrom => {
            sentinelFiles +:= createScanJobs(chrom)
        })

        addCommand(new MergeResults(sentinelFiles))
    }

    def validateProfilesDirectories() {
        profilesDirList.foreach(profilesDir =>
            if (!profilesDir.exists) {
                throw new RuntimeException("Profiles directory " + profilesDir + " does not exist")
            }
        )
    }

    def parseBinSize() {
        val profilesDir = profilesDirList.head
        val profiles = profilesDir.listFiles.filter(_.isFile).toList.filter { file =>
            file.getName.endsWith(".dat.gz") && file.length > 0
        }
        if (profiles.size == 0) {
            throw new RuntimeException(profilesDir + " does not contain any files with the extension dat.gz")
        }

        val gzStream = new GZIPInputStream(new FileInputStream(profiles.head))
        val src = Source.fromInputStream(gzStream)
        binSize = src.getLines().drop(1).take(1).map(line => line.split("\t", -1)(4)).toList.head.toInt
        gzStream.close()
        src.close()
    }

    def createSampleList() {
         profilesDirList.foreach(profilesDir => {
            val profiles = profilesDir.listFiles.filter(_.isFile).toList.filter { file =>
                file.getName.endsWith(".dat.gz") && file.length > 0
            }
            if (profiles.size == 0) {
                throw new RuntimeException(profilesDir + " does not contain any files with the extension dat.gz")
            }

            val gzStream = new GZIPInputStream(new FileInputStream(profiles.head))
            val src = Source.fromInputStream(gzStream)

            val fields = src.getLines().next.split("\t")
            val profileSamples = fields.slice(6, fields.size)
            gzStream.close()
            src.close()

            if (sampleList.intersect(profileSamples).size > 0) {
                throw new RuntimeException("Sample sets in the profile directories are not disjoint")
            }
            sampleList ++= profileSamples
        })
    }


    def sequenceList = {
        if (genomeIntervalList != null) {
            expandListFiles(genomeIntervalList)
        } else if (getReferenceMetadata.profileIntervalBed != null) {
            parseBedFileIntervalMap(getReferenceMetadata.profileIntervalBed)
        } else {
            computeReferenceSequenceNames
        }
    }

    def parseBedFileIntervalMap(bedFile: File) = {
        (new BedFileReader(bedFile) : java.lang.Iterable[BedFileLine]).asScala.toList
            .map { line => line.getField(0) }
    }

    def sentinelFileDirectory = {
        new File(runDirectory, "lcnv_sentinel_files")
    }

    def sentinelFile(chrom: String) = {
        new File(sentinelFileDirectory, "seq_" + chrom + ".sent")
    }

    def createScanJobs(chrom: String): File = {
        var scanFiles: List[File] = Nil
        if (perSampleScan) {
            val chromDirectory = new File(runDirectory, "seq_" + chrom)
            createDirectory(chromDirectory)
            val samples = CommandLineParser.parseStringList(sampleList).asScala.toList
            samples.foreach(sample => {
                // Strip '/' from the sample name
                val sampleName = ("/".r).replaceAllIn(sample, "")
                val scanFile = new File(chromDirectory, chrom + "_" + sampleName + ".dat")
                scanFiles +:= addCommand(new ScanCommand(chrom, List(sample), scanFile))
            })
        } else {
            val scanFile = new File(runDirectory, "seq_" + chrom + ".dat")
            scanFiles +:= addCommand(new ScanCommand(chrom, sampleList, scanFile))
        }

        addCommand(new TouchSentinelFile(sentinelFile(chrom), scanFiles))
    }

    class ScanCommand(chrom: String, targetSampleList: List[String], scanFile: File) extends CommandLineFunction
        with FileInputOutput {

        this.outputFile = scanFile

        val jobLogDir = new File(scanFile.getParent, "/logs")
        this.jobOutputFile = new File(jobLogDir, scanFile.getName().stripSuffix(".dat") + ".out")
        this.memoryLimit = Some(4)
        this.jobArrayName = "LCNV"

        val profileFileList = profilesDirList.map(profilesDir => new File(profilesDir, "profile_seq_" + chrom + "_" + binSize + ".dat.gz"))

        override def commandLine =
            required(jobWrapperScript) +
            required("Rscript") +
            required(rScriptPath) +
            repeat("--profileFile", profileFileList) +
            repeat("--targetSample", targetSampleList) +
            repeat("--backgroundSample", backgroundSampleList) +
            required("--maxDepth", maxDepth) +
            required("--ploidyMapFile", getReferenceMetadata.ploidyMap) +
            repeat("--genderMapFile", getGenderMapFileList) +
            required("--outputFile", scanFile)
    }

    class MergeResults(sentinelFiles: List[File]) extends CommandLineFunction 
        with FileInputOutput {
        val callsFile = runDirectory.getAbsolutePath + "/calls.dat"
        this.dependsOnFile ++= sentinelFiles
        this.outputFile = callsFile
        this.jobOutputFile = new File(jobLogDirectory, "MergeResults.out")

        override def commandLine =
            " find " + runDirectory + "/seq_* " +
            " -name '*.dat' | head -1 | xargs -d '\\n' head -1 " +
            " > " + callsFile + " || exit 1; " +
            " find " + runDirectory + "/seq_* " +
            " -name '*.dat' -exec tail -n +2 {} \\; " +
            " | sed 's/^X\\t/23\\t/g' " +
            " | sort -k1n -k2n -k3n " +
            " | sed 's/^23\\t/X\\t/g' " +
            " >> " + callsFile + " || exit 1; "
    }
}
