
import org.broadinstitute.sv.qscript.SVQScript

import org.broadinstitute.sv.queue.ComputeVCFPartitions

class SVGenotyper extends SVQScript {

    @Input(fullName="vcfFile", shortName="vcf", doc="The vcf file of events to genotype.")
    var vcfFile: File = null

    @Output(shortName="O", doc="The output vcf file.")
    var outputFile: File = null

    @Argument(shortName="partition", required=false, doc="Specific partitions to rerun")
    var partitionList: List[String] = null

    @Argument(fullName="genotypeFilterFile", shortName="genotypeFilterFile", doc="File(s) conatining per-sample regions to be excluded from genotyping", required=false)
    var genotypeFilterFileList: List[File] = Nil

    @Argument(shortName="skipAnnotator", required=false, doc="Genotyping annotator(s) that should not be run as part of the default annotators set")
    var skipAnnotatorsList: List[String] = Nil

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        val suffix = if (outputFile.endsWith(".gz") ) ".gz" else ""
        val outputFileBaseName = outputFile.getName.stripSuffix(".gz").stripSuffix(".vcf")
        val runFilePrefix = outputFileBaseName.stripSuffix(".genotypes")
        val unfilteredOutFile = new File(runDirectory, outputFileBaseName + ".unfiltered.vcf" + suffix)
        val annotatedOutFile = new File(runDirectory, outputFileBaseName + ".annotated.vcf" + suffix)
        val partitions = computeVCFPartitions(vcfFile)
        if (partitions.isEmpty) {
            addCommand(new SVGenotyper(vcfFile, unfilteredOutFile, runFilePrefix, genotypeFilterFileList))
        } else {
            var gtPartFiles: List[File] = Nil
            for ((partitionName, partitionArg) <- partitions) {
                if (partitionList == null || partitionList.contains(partitionName)) {
                    val pGenotyper = new SVParallelGenotyper(vcfFile, partitionName, partitionArg, genotypeFilterFileList)
                    pGenotyper.jobArrayName = "Genotypers"
                    gtPartFiles :+= addCommand(pGenotyper)
                }
            }
            addCommand(new MergeGenotyperOutput(vcfFile, unfilteredOutFile, gtPartFiles, runFilePrefix))
        }
        // Added CopyNumberClassAnnotator to GenotyperDefaultFilterAnnotations as to
        addCommand(new GenotyperDefaultFilterAnnotations(unfilteredOutFile, annotatedOutFile, Nil, skipAnnotatorsList))
        addCommand(new SVGenotyperDefaultFilter(annotatedOutFile, outputFile))
        val qcAnnotationsReport = addCommand(new GenotyperDefaultQCAnnotations(outputFile, null))
        addCommand(new TouchSentinelFile(new File(runDirectory, "SVGenotyper2.sent"), qcAnnotationsReport))
        //addCommand(new CopyNumberClassAnnotator(outputFile))
    }
}
