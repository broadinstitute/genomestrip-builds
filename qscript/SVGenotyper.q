
import org.broadinstitute.sv.qscript.SVQScript
import org.broadinstitute.sv.queue.ComputeVCFPartitions
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature

//@QscriptDocStart
@DocumentedGATKFeature(groupName = "Queue Scripts")
class SVGenotyper extends SVQScript {

    @Input(fullName="vcfFile", shortName="vcf", doc="The vcf file of events to genotype.")
    var vcfFile: File = null

    @Output(shortName="O", doc="The output vcf file.")
    var outputFile: File = null

    @Argument(shortName="partition", required=false, doc="Specific partitions to rerun")
    var partitionList: List[String] = null

    @Argument(shortName="genotypeFilterFile", fullName="genotypeFilterFile", doc="File of filter regions", required=false)
    var genotypeFilterFileList: List[File] = Nil


//@QscriptDocEnd

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        var gtOutFile: File = null
        val partitions = computeVCFPartitions(vcfFile)
        if (partitions.isEmpty) {
            gtOutFile = addCommand(new SVGenotyper(vcfFile, outputFile, null, genotypeFilterFileList))
        } else {
            var gtPartFiles: List[File] = Nil
            for ((partitionName, partitionArg) <- partitions) {
                if (partitionList == null || partitionList.contains(partitionName)) {
                    val pGenotyper = new SVParallelGenotyper(vcfFile, partitionName, partitionArg, genotypeFilterFileList)
                    pGenotyper.jobArrayName = "Genotypers"
                    gtPartFiles :+= addCommand(pGenotyper)
                }
            }
            gtOutFile = addCommand(new MergeGenotyperOutput(vcfFile, outputFile, gtPartFiles))
        }
    }
}
