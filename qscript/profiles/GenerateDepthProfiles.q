/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2015 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */

import scala.collection.JavaConverters._
import org.broadinstitute.sv.qscript.SVQScript
import org.broadinstitute.sv.util.GenomeInterval
import org.broadinstitute.sv.util.bed.BedFileLine
import org.broadinstitute.sv.util.bed.BedFileReader

class GenerateDepthProfiles extends SVQScript {

    @Argument(shortName="profileBinSize", required=false, doc="Size of profile bins")
    var profileBinSize: java.lang.Integer = 100;

    @Argument(shortName="maximumReferenceGapLength", required=false, doc="Max reference gap length")
    var maximumReferenceGapLength: java.lang.Integer = null


    class ComputeDepthProfile(profileDir: File, sequenceName: String, intervalList: List[GenomeInterval]) extends JavaCommand with BAMInputOutput {
        this.javaMainClass = "org.broadinstitute.sv.apps.ComputeDepthProfiles"
        this.outputFile = new File(profileDir, "profile_seq_%s_%d.dat.gz".format(sequenceName, profileBinSize))
        this.jobArrayName = "Profiles"

        commandArguments +=
            repeat(" -I ", bamLocations) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            required(" -R ", referenceFile) +
            repeat(" -L ", intervalList) +
            repeat(" -md ", metaDataLocationList) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            required(" -profileBinSize ", profileBinSize) +
            optional(" -maximumReferenceGapLength ", maximumReferenceGapLength)
    }

    // Generate a tabix index for the profile file
    class IndexDepthProfile(profileFile: File)  extends SimpleCommand {
        val indexFile = new File(profileFile.getPath() + ".tbi")
        this.inputFile :+= profileFile
        this.outputFile = indexFile

        this.jobArrayName = "Indices"

        commandArguments = "tabix -f -S 1 -s 2 -b 3 -e 4 %s".format(profileFile)
    }

    def profileIntervalMap = {
        // Currently, we use one of three methods, in priority order.
        // If -L was supplied, generate one profile per chromosome (for every one mentioned in -L).
        // If we have reference metadata interval list (as a bed file), use the name field to build the interval map.
        // Otherwise, fall back to the historical method of one profile per reference sequence.
        if (genomeIntervalList != null) {
            expandListFiles(genomeIntervalList).map { interval => GenomeInterval.parse(interval) }
                .map { genomeInterval => genomeInterval.getSequenceName() -> genomeInterval }
                .groupBy(_._1).map { case (k,v) => (k, v.map(_._2)) }
        } else if (getReferenceMetadata.profileIntervalBed != null) {
            parseBedFileIntervalMap(getReferenceMetadata.profileIntervalBed)
        } else {
            computeReferenceSequenceNames.map { x => x -> List(GenomeInterval.parse(x)) }.toMap
        }
    }

    def parseBedFileIntervalMap(bedFile: File) = {
        (new BedFileReader(bedFile) : java.lang.Iterable[BedFileLine]).asScala.toList
            .map { line => line.getField(3) -> line.getInterval }
            .groupBy(_._1).map { case (k,v) => (k, v.map(_._2)) }
    }

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        val profileDir = runDirectory
        var profileIndexFiles: List[File] = Nil
        for ((sequenceName, intervalList) <- profileIntervalMap) {
            var profileFile = addCommand(new ComputeDepthProfile(profileDir, sequenceName, intervalList))
            profileIndexFiles :+= addCommand(new IndexDepthProfile(profileFile))
        }
    }
}
