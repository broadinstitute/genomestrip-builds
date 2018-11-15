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
import collection.JavaConversions

import org.broadinstitute.sv.util.FileUtilities
import org.broadinstitute.sv.util.io.ErrorCheckingPrintWriter

class CNVDiscoveryStage7 extends CNVDiscoveryStageBase {

    @Input(shortName="siteListFile", required=true, doc="List of sites to process")
    var siteListFile: File = null

    @Argument(shortName="boundaryPrecision", required=false, doc="Target CNV boundary refinement precision (in effective length, default 200)")
    var boundaryPrecision: java.lang.Integer = null

    @Argument(shortName="maximumReferenceGapLength", required=false, doc="Do not generate event boundaries containing gaps larger than this (base pairs)")
    var maximumReferenceGapLength: java.lang.Integer = null

    @Argument(shortName="minimumRefinedLength", required=false, doc="Minimum length of refined interval (effective length)")
    var minimumRefinedLength: java.lang.Integer = null

    @Argument(shortName="maximumRefinedLength", required=false, doc="Maximum length of refined interval (effective length)")
    var maximumRefinedLength: java.lang.Integer = null

    @Argument(shortName="minimumRefinedLengthFraction", required=false, doc="Minimum relative length of refined interval (effective length)")
    var minimumRefinedLengthFraction: java.lang.Double = null

    @Argument(shortName="maximumRefinedLengthFraction", required=false, doc="Maximum relative length of refined interval (effective length)")
    var maximumRefinedLengthFraction: java.lang.Double = null

    @Argument(shortName="brigSiteVcfFile", required=true, doc="Brig sites vcf file")
    var brigSiteVcfFile: File = null

    @Argument(shortName="produceAuxiliaryFiles", required=false, doc="Flag specifying whether to produce auxiliary files and keep them")
    var produceAuxiliaryFiles: Boolean = false

    val numSitesPerJob = 10
    var partialBrigVcfFileName: String = null
    var brigVcfFileList: File = null
    var brigOutputFiles: List[File] = Nil

    def script = {
        partialBrigVcfFileName = vcfFile.getName.stripSuffix(".gz").stripSuffix(".vcf").stripSuffix(".sites").stripSuffix(".genotypes") + ".brig.vcf"

        brigVcfFileList = new File(runDirectory, "brig.vcf.file.list")
        addCommand(new CreatePartialVcfFileList)

        val siteList = FileUtilities.parseStringList(siteListFile)
        val numPartitions = Math.ceil(1.0 * siteList.size / numSitesPerJob).toInt
        (1 to numPartitions).foreach(partitionNumber => {
            val partSiteListFile = new File(runDirectory, generatePartitionId(partitionNumber) + ".sites.list")
            val partitionDir = new File(runDirectory, generatePartitionId(partitionNumber))
            val brigOutputFile = getPartialBrigVcfFile(partitionDir)

            val boundaryRefiner = new RefineCNVBoundaries(vcfFile, partSiteListFile, brigOutputFile)
            boundaryRefiner.jobArrayName = "BoundaryRefiner"
            addCommand(boundaryRefiner)
            brigOutputFiles :+= brigOutputFile
        })

        val mergeVcfFiles = addCommand(new MergeBrigVcfFiles)
        val tempFileSentinel = deleteBrigTempFiles(mergeVcfFiles)
        addCommand(new TouchSentinelFile(sentinelFile, tempFileSentinel))
    }

    class CreatePartialVcfFileList() extends InProcessFunction
        with LogDirectory with UnifiedInputOutput {
        this.commandResultFile = brigVcfFileList

        def run = {
            var numSites = 0
            var partitionNumber = 0
            var partitionDir: File = null
            var partSiteListWriter: PrintWriter = null
            val brigVcfFileWriter = new ErrorCheckingPrintWriter(brigVcfFileList)
            val siteList = FileUtilities.parseStringList(siteListFile)
            for (site <- JavaConversions.iterableAsScalaIterable(siteList)) {
                numSites += 1
                if (partSiteListWriter == null) {
                    partitionNumber += 1
                    val partSiteListFile = new File(runDirectory, generatePartitionId(partitionNumber) + ".sites.list")
                    partSiteListWriter = new ErrorCheckingPrintWriter(partSiteListFile)
                    partitionDir = createPartitionDir(partitionNumber)
                }
                partSiteListWriter.println(site)
                if (numSites % numSitesPerJob == 0 || numSites == siteList.size) {
                    val brigOutputFile = getPartialBrigVcfFile(partitionDir)
                    brigVcfFileWriter.println(brigOutputFile)
                    partSiteListWriter.close
                    partSiteListWriter = null
                }
            }
            brigVcfFileWriter.close
        }
    }

    def createPartitionDir(partitionNumber: Int): File = {
        val partitionDir = new File(runDirectory, generatePartitionId(partitionNumber))
        createDirectory(partitionDir)
        partitionDir
    }

    def generatePartitionId(n: Int): String = {
        var number = n.toString
        while (number.length() < 4) {
            number = "0" + number
        }
        "P" + number
    }

    def getPartialBrigVcfFile(partitionDir: File) = {
        new File(partitionDir, partialBrigVcfFileName)
    }

    // Notes:
    // -auxFilePrefix is used only for generating intermediate files for debugging (and plotting?)
    // -tempDir should not be used if we are using the internal java EM
    class RefineCNVBoundaries(vcfFile: File, siteListFile: File, outputVcf: File) extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.genotyping.RefineCNVBoundaries"

        var auxFilePrefix: String = null
        if (produceAuxiliaryFiles) {
            auxFilePrefix = outputVcf.getPath.stripSuffix(".brig.vcf")
        }

        this.dependsOnFile :+= brigVcfFileList
        this.inputFile = bamInputs
        this.outputFile = outputVcf

        commandArguments +=
            required(" -R ", referenceFile) +
            repeat(" -md ", metaDataLocationList) +
            repeat(" -configFile ", parameterFiles) +
            optional(" -P ", "depth.readCountCacheIgnoreGenomeMask:true") +
            repeat(" -P ", parameterList) +
            repeat(" -genotypeFilterFile ", genotypeFilterFileList) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -genderMapFile ", getGenderMapFileList) +
            optional(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
            repeat(" -excludeReadGroup ", excludeReadGroup) +
            required(" -vcf ", vcfFile) +
            required(" -site ", siteListFile) +
            optional(" -auxFilePrefix ", auxFilePrefix) +
            optional(" -boundaryPrecision ", boundaryPrecision) +
            optional(" -minimumRefinedLength ", minimumRefinedLength) +
            optional(" -maximumRefinedLength ", maximumRefinedLength) +
            optional(" -minimumRefinedLengthFraction ", minimumRefinedLengthFraction) +
            optional(" -maximumRefinedLengthFraction ", maximumRefinedLengthFraction) +
            optional(" -maximumReferenceGapLength ", maximumReferenceGapLength)
    }

    class MergeBrigVcfFiles extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.discovery.MergeBrigVcfFiles"
        this.dependsOnFile = brigOutputFiles
        this.commandResultFile = brigSiteVcfFile

        commandArguments +=
            required(" -R ", referenceFile) +
            required(" -vcfFile ", brigVcfFileList) +
            required(" -mergedVcfFile ", brigSiteVcfFile)
    }

    def deleteBrigTempFiles(dependsOn: File) : File = {
        var result = dependsOn
        if (shouldDeleteTempFiles("BRIGPARTITIONS")) {
            result = addCommand(new DeleteBrigPartitions(result))
        }
        if (shouldDeleteTempFiles("BRIGLISTS")) {
            result = addCommand(new DeleteBrigListFiles(result))
        }
        result
    }

    class DeleteBrigPartitions(dependsOn: File) extends SimpleCommand {
        this.dependsOnFile +:= dependsOn
        this.outputFile = new File(runDirectory, "DeleteBrigPartitions.sent")
        commandArguments =
            "find " + runDirectory +
            " -maxdepth 1 -type d -name 'P[0-9]*' " +
            " -exec rm -r {} \\;"
    }

    class DeleteBrigListFiles(dependsOn: File) extends SimpleCommand {
        this.dependsOnFile +:= dependsOn
        this.outputFile = new File(runDirectory, "DeleteBrigListFiles.sent")
        commandArguments =
            "find " + runDirectory +
            " -maxdepth 1 -type f -name 'P[0-9]*.sites.list' " +
            " -exec rm {} \\;"
    }
}
