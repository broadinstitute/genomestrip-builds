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
import collection.JavaConversions._
import scala.io.Source

import org.broadinstitute.sv.util.CommandRunner
import org.broadinstitute.sv.util.io.ErrorCheckingPrintWriter

class CNVDiscoveryStage12 extends CNVDiscoveryGenotyper {
    // Note: We inherit from CNVDiscoveryGenotyper to get access to the CNVDefaultQCAnnotations inner class
    // We do not re-genotype in stage 12.

    @Input(shortName="filteredVcf", required=true, doc="The list of per-sequence vcf files to merge")
    var vcfFileList: List[File] = Nil

    @Output(shortName="mergedVcfFile", fullName="mergedVcfFile", required=true, doc="The output vcf file.")
    var mergedVcfFile: File = null

    @Argument(shortName="outputDirectory", required=true, doc="Where to write the final pipeline results.")
    var outputDirectory: File = null

    /**
     * In script, you create and add functions to the pipeline.
     */
    override def script = {
        val mergedFile = addCommand(new MergeVcfFiles)
        val annotations = addCommand(new CNVFinalQCAnnotations(mergedFile, null))
        addCommand(new CreatePipelineResults(annotations))
    }

    class MergeVcfFiles extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.MergeVCFFiles"
        this.outputFile = mergedVcfFile

        commandArguments +=
            repeat(" -vcf ", vcfFileList) +
            required(" -R ", referenceFile) +
            required(" -filterVariants ", "true")
    }

    class CNVFinalQCAnnotations(vcfFile: File, outFile: File)
        extends CNVDefaultQCAnnotations(vcfFile, outFile) {
        commandArguments +=
            required(" -A ", "CopyNumberClass")
    }

    class CreatePipelineResults(dependsOn: File) extends InProcessFunction
        with LogDirectory with UnifiedInputOutput {

        this.outputFile = new File(new File(outputDirectory, "results"), "gs_cnv.genotypes.vcf.gz")
        this.dependsOnFile :+= dependsOn

        def run = {
            val resultsOutputDir = new File(outputDirectory, "results")
            createDirectory(resultsOutputDir)
            copyFile(new File(runDirectory, "gs_cnv.genotypes.vcf.gz"), new File(resultsOutputDir, "gs_cnv.genotypes.vcf.gz"))
            copyFile(new File(runDirectory, "gs_cnv.genotypes.vcf.gz.tbi"), new File(resultsOutputDir, "gs_cnv.genotypes.vcf.gz.tbi"))
            copyDirectory(new File(runDirectory, "eval"), new File(resultsOutputDir, "eval"))

            val midStageSourceDir = new File(runDirectory.getParent(), "cnv_stage5")
            val midStageOutputDir = new File(outputDirectory, "midstage")
            createDirectory(midStageOutputDir)
            copyDirectory(new File(midStageSourceDir, "eval"), midStageOutputDir)

            // Temporary:  We create a paritition map to the stage10 output directory.
            // This is a temporary measure to allow plotting of genotyped sites (via PlotGenotypingResults).
            // We should adapt the plotting code to work off of the output VCF directly.
            val stage10RunDir = new File(runDirectory.getParent(), "cnv_stage10")
            val partitionMap = new File(resultsOutputDir, "partition.genotypes.map.dat")
            createPartitionMapFile(stage10RunDir, partitionMap)
        }
    }

    def copyFile(sourceFile: File, destFile: File) {
        val commandRunner = new CommandRunner()
        commandRunner.setMergeOutput(true)
        val status = commandRunner.runCommand("cp %s %s".format(sourceFile.getPath(), destFile.getPath()))
        if (status != 0) {
            throw new RuntimeException("Copy failed: " + commandRunner.getStandardOutputString())
        }
    }

    def copyDirectory(sourceFile: File, destFile: File) {
        val commandRunner = new CommandRunner()
        commandRunner.setMergeOutput(true)
        val status = commandRunner.runCommand("cp -rT %s %s".format(sourceFile.getPath(), destFile.getPath()))
        if (status != 0) {
            throw new RuntimeException("Copy failed: " + commandRunner.getStandardOutputString())
        }
    }

    def createPartitionMapFile(stage10RunDir: File, outFile: File) {
        val partitionMapWriter = new ErrorCheckingPrintWriter(outFile)
        partitionMapWriter.println("CNP\tPARTITION\tCHR\tSTART\tEND")
        for (seqDirName <- stage10RunDir.list()) {
            val seqDir = new File(stage10RunDir, seqDirName)
            val inputMapFile = new File(seqDir, "partition.genotypes.map.dat")
            val partition = seqDir.getAbsolutePath() + "/" + seqDirName + ".adjacent_merged"
            System.out.println("#DBG: inFile " + inputMapFile + " partition " + partition)
            for (line <- Source.fromFile(inputMapFile).getLines()) {
                val fields = line.split("\t")
                if (fields(0) != "CNP") {
                    fields(1) = partition
                    partitionMapWriter.println(fields.mkString("\t"))
                }
            }
        }
        partitionMapWriter.close
    }
}
