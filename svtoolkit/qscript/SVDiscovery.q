
import org.broadinstitute.sv.qscript.SVQScript
import org.broadinstitute.sv.queue.ComputeDiscoveryPartitions
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature

//@QscriptDocStart
@DocumentedGATKFeature(groupName = "Queue Scripts")
class SVDiscovery extends SVQScript {

    @Output(fullName="outputFile", shortName="O", required=true, doc="The output vcf file")
    var outputFile: File = null

    @Argument(fullName="partition", shortName="partition", required=false, doc="Specific partitions to rerun")
    var partitionList: List[String] = null

    @Argument(shortName="storeReadPairFile", required=false, doc="Indicates whether to save reads into BAM file (defaults to true)")
    var storeReadPairFileArg: String = null
    var storeReadPairFile: Boolean = true

//@QscriptDocEnd

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        if (storeReadPairFileArg != null) {
            storeReadPairFile = parseBoolean(storeReadPairFileArg)
        }

        val suffix = if (outputFile.endsWith(".gz") ) ".gz" else ""
        var runFilePrefix = outputFile.getName.stripSuffix(".gz").stripSuffix(".vcf").stripSuffix(".discovery")
        var unfilteredOutFile = new File(runDirectory, runFilePrefix + ".unfiltered.vcf" + suffix)
        val partitions = computeDiscoveryPartitions()
        if (partitions.size < 2) {
            addCommand(new SVDiscovery(unfilteredOutFile, runFilePrefix, storeReadPairFile))
        } else {
            var discoPartFiles: List[File] = Nil
            for ((partitionName, partitionArgs) <- partitions) {
                if (partitionList == null || partitionList.contains(partitionName)) {
                    val pDiscovery = new SVParallelDiscovery(partitionName, partitionArgs, storeReadPairFile)
                    pDiscovery.jobArrayName = "Discovery"
                    discoPartFiles :+= addCommand(pDiscovery)
                }
            }
            addCommand(new MergeDiscoveryOutput(unfilteredOutFile, discoPartFiles))
            if (storeReadPairFile) {
                val mergedBamFile = new File(runDirectory, outputFile.getName.stripSuffix(".gz").stripSuffix(".vcf").stripSuffix(".sites") + ".pairs.bam")
                addCommand(new MergePartitionPairBamFiles(mergedBamFile, discoPartFiles))
                addCommand(new CreateBamIndex(mergedBamFile))
                addCommand(new DeletePartitionPairBamFiles(mergedBamFile))
            }
        }
        addCommand(new SVDiscoveryDefaultFilter(unfilteredOutFile, outputFile))
    }

    class MergePartitionPairBamFiles(mergedBamFile: File, discoPartFiles: List[File]) extends SimpleCommand {
        this.dependsOnFile ++= discoPartFiles
        this.outputFile = mergedBamFile

        val builder = new StringBuilder
        builder.append("samtools merge -c -p %s".format(mergedBamFile))
        discoPartFiles.foreach(f =>
            builder.append(" %s".format(new File(runDirectory, f.getName.stripSuffix(".gz").stripSuffix(".vcf") + ".pairs.bam")))
        )
        commandArguments = builder.toString
    }

    class DeletePartitionPairBamFiles(mergedBamFile: File) extends SimpleCommand {
        this.dependsOnFile +:= mergedBamFile
        this.outputFile = new File(runDirectory, "DeletePartitionPairBamFiles.sent")
        commandArguments =
            "find " + runDirectory +
            " -maxdepth 1 -type f -name 'P[0-9]*.discovery.pairs.bam' " +
            " -exec rm {} \\;"
    }
}
