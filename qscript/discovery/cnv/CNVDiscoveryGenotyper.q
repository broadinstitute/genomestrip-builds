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

abstract class CNVDiscoveryGenotyper extends CNVDiscoveryStageBase {

    @Input(shortName="genderGenotypeFilterFile", required=true, doc="File specifying what genotypes are to be filtered out")
    var genderGenotypeFilterFile: File = null

    @Input(shortName="filterDescriptionFile", required=true, doc="File specifying the genotype filter descriptions")
    var filterDescriptionFile: File = null

    @Argument(shortName="genotypingParallelRecords", required=false, doc="The number of parallel records to use in genotyping")
    var genotypingParallelRecords: java.lang.Integer = 1000

    def isFinalStage = false
    def applyDuplicateFilter = true

    def script = {

        parameterList :+= "genotyping.modules:depth"
        parameterList :+= "depth.readCountCacheIgnoreGenomeMask:true"

        val genotypesVcfFileName = swapExt(vcfFile.getName, "sites.vcf.gz", "genotypes.vcf.gz")
        val genotypesVcfBaseName = genotypesVcfFileName.stripSuffix(".gz").stripSuffix(".vcf")
        val genotypesVcfFile = new File(runDirectory, genotypesVcfFileName)

        val runFilePrefix = genotypesVcfBaseName.stripSuffix(".genotypes")
        val unfilteredVcfFile = new File(runDirectory, genotypesVcfBaseName + ".unfiltered.vcf.gz")
        val baseAnnotatedVcfFile = new File(runDirectory, genotypesVcfBaseName + ".base.annotated.vcf.gz")
        val annotatedVcfFile = new File(runDirectory, genotypesVcfBaseName + ".annotated.vcf.gz")
        val partitions = computeVCFPartitions(vcfFile, null, genotypingParallelRecords)
        if (partitions.isEmpty) {
            addCommand(new SVGenotyper(vcfFile, unfilteredVcfFile, runFilePrefix, genotypeFilterFileList))
        } else {
            var gtPartFiles: List[File] = Nil
            for ((partitionName, partitionArg) <- partitions) {
                val pGenotyper = new SVParallelGenotyper(vcfFile, partitionName, partitionArg, genotypeFilterFileList)
                pGenotyper.jobArrayName = "Genotypers"
                gtPartFiles :+= addCommand(pGenotyper)
            }
            addCommand(new MergeGenotyperOutput(vcfFile, unfilteredVcfFile, gtPartFiles, runFilePrefix))
        }

        addCommand(new CNVBaseFilterAnnotations(unfilteredVcfFile, baseAnnotatedVcfFile))
        addCommand(new CNVDefaultFilterAnnotations(baseAnnotatedVcfFile, annotatedVcfFile))
        addCommand(new CNVDefaultFilter(annotatedVcfFile, genotypesVcfFile))
        val siteFiltersSummaryFile = addCommand(new CNVDefaultQCAnnotations(genotypesVcfFile, null))
        val classifyVariants = addCommand(new ClassifySelectedVariants(List(annotatedVcfFile, siteFiltersSummaryFile)))
        val tempFileSentinel = deleteGenotypingTempFiles(genotypesVcfFile, classifyVariants)
        addCommand(new TouchSentinelFile(sentinelFile, tempFileSentinel))
    }

    val cnvBaseAnnotatorList = List(
        "ClusterSeparation",
        "GenotypeFilter"
    )

    class CNVBaseFilterAnnotations(vcfFile: File, outFile: File)
        extends SVAnnotator(vcfFile, outFile) {

        var annotatorList = cnvBaseAnnotatorList
        if (isFinalStage) {
            annotatorList +:= "CNQuality"
        }

        commandArguments +=
            annotatorList.map(e => "-A " + e).mkString(" ") +
            required(" -genotypeFilterFile ", genderGenotypeFilterFile) +
            required(" -filterDescriptionFile ", filterDescriptionFile) +
            required(" -filterVariants ", "false") +
            required(" -writeReport ", "true") +
            required(" -writeSummary ", "true")
    }

    val cnvDefaultAnnotatorList = List(
        "GCContent",
        "GenotypeLikelihoodStats",
        "NonVariant",
        "Redundancy",
        "CopyNumberClass"
    )

    class CNVDefaultFilterAnnotations(vcfFile: File, outFile: File)
        extends SVAnnotator(vcfFile, outFile) {

        commandArguments +=
            cnvDefaultAnnotatorList.map(e => "-A " + e).mkString(" ") +
            required(" -filterVariants ", "false") +
            required(" -writeReport ", "true") +
            required(" -writeSummary ", "true") +
            required(" -comparisonFile ", vcfFile) +
            required(" -duplicateOverlapThreshold ", duplicateOverlapThreshold) +
            required(" -duplicateScoreThreshold ", duplicateScoreThreshold)

        if (isFinalStage && getReferenceMetadata.vdjRegionsBed != null && getReferenceMetadata.vdjRegionsBed.exists) {
            commandArguments +=
                required(" -A ", "BedFileOverlap") +
                required(" -bedFileTag ", "VDJ") +
                required(" -bedFile ", getReferenceMetadata.vdjRegionsBed)
        }
    }

    class CNVDefaultFilter(vcfFile: File, outFile: File) extends SVVariantFiltration(vcfFile, outFile) {
        commandArguments +=
            " -filterName CLUSTERSEP -filter \"GSCLUSTERSEP == \\\"NA\\\" || GSCLUSTERSEP <= 2.0\"" +
            " -filterName GTDEPTH -filter \"GSM1 == \\\"NA\\\" || GSM1 <= 0.5 || GSM1 >= 2.0\"" +
            " -filterName NONVARIANT -filter \"GSNONVARSCORE != \\\"NA\\\" && GSNONVARSCORE >= " + nonVariantScoreThreshold + "\""
        if (applyDuplicateFilter) {
            commandArguments +=
                " -filterName DUPLICATE -filter \"GSDUPLICATESCORE != \\\"NA\\\" && GSDUPLICATEOVERLAP >= " + duplicateOverlapThreshold + " && GSDUPLICATESCORE >= " + duplicateScoreThreshold + "\""
        }
    }

    class CNVDefaultQCAnnotations(vcfFile: File, outFile: File)
        extends SVAnnotator(vcfFile, outFile) {
        val siteFiltersReportFile = new File(evalDir, "GenotypeSiteFilters.report.dat")
        val siteFiltersSummaryFile = new File(evalDir, "GenotypeSiteFilters.summary.dat")
        this.commandResultFile = siteFiltersSummaryFile
        commandArguments +=
            required(" -A ", "SiteFilters") +
            required(" -A ", "VariantsPerSample") +
            required(" -reportFileMap ", "SiteFilters:" + siteFiltersReportFile.getPath) +
            required(" -summaryFileMap ", "SiteFilters:" + siteFiltersSummaryFile.getPath) +
            required(" -writeReport ", "true") +
            required(" -writeSummary ", "true")
    }

    def deleteGenotypingTempFiles(vcfFile: File, dependsOn: File) : File = {
        val fileNameRoot = stripExt(vcfFile, "genotypes.vcf.gz").getName()
        var fileNamePatterns: List[String] = Nil
        if (shouldDeleteTempFiles("GTPARTITIONS")) {
            fileNamePatterns :+= " -name 'P[0-9]*.genotypes.*' "
        }
        if (shouldDeleteTempFiles("GTFILTERING")) {
            fileNamePatterns :+= " -name '%s.genotypes.unfiltered.*' ".format(fileNameRoot)
            fileNamePatterns :+= " -name '%s.genotypes.base.annotated.*' ".format(fileNameRoot)
            fileNamePatterns :+= " -name '%s.genotypes.annotated.*' ".format(fileNameRoot)
        }
        if (!isFinalStage && shouldDeleteTempFiles("GTAUXFILES")) {
            fileNamePatterns :+= " -name '%s.genotypes.*.dat' ".format(fileNameRoot)
        }
        if (fileNamePatterns != Nil) {
            val patterns = "\\( " + fileNamePatterns.mkString(" -o ") + " \\)"
            addCommand(new DeleteTempFiles(patterns, dependsOn))
        } else {
            dependsOn
        }
    }

    class DeleteTempFiles(patterns: String, dependsOn: File) extends SimpleCommand {
        this.dependsOnFile +:= dependsOn
        this.outputFile = new File(runDirectory, "DeleteTempFiles.sent")
        commandArguments =
            "find " + runDirectory +
            " -maxdepth 1 -type f " + patterns +
            " -exec rm {} \\;"
    }
}
