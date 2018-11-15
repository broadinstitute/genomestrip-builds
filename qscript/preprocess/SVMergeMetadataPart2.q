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
package org.broadinstitute.sv.qscript.preprocess

import org.broadinstitute.sv.qscript.SVQScript
import org.broadinstitute.sv.util.GenomeInterval
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature

// This script runs the last part of metadata merging that generates 100kb profiles and does gender calling.
// This is a start towards a complete Queue solution for metadata merging.
// Unfortunately, doing this in Queue requires either extensive refactoring or lots of duplicated code.
// Currently, unless we fix Queue, you cannot subclass a Queue class that defines "script".
// So, all of the code in SVPreprocess either has to be copied or refactored and moved to some other class.

//@QscriptDocStart
//@DocumentedGATKFeature(groupName = "Queue Scripts")
class SVMergeMetadataPart2 extends SVQScript {

    @Argument(shortName="profileBinSize", required=false, doc="Size of profile bins to use for depth profiles")
    var profileBinSize: java.lang.Integer = 100000;

    @Argument(shortName="maximumReferenceGapLength", required=false, doc="Max reference gap in depth profiles")
    var maximumReferenceGapLength: java.lang.Integer = 10000

//@QscriptDocEnd

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        createProfilesCallGenders
    }

    def createProfilesCallGenders = {
        var profilesDirName: String = "profiles_"
        if (profileBinSize < 1000) {
            profilesDirName += profileBinSize
        } else {
            profilesDirName += (profileBinSize / 1000) + "Kb"
        }
        val profilesDir = new File(metaDataLocation, profilesDirName)
        createDirectory(profilesDir)
        var profileIndexFiles: List[File] = Nil
        for ((sequenceName, intervalList) <- intervalMap) {
            var profileFile = addCommand(new ComputeDepthProfile(profilesDir, sequenceName, intervalList))
            profileIndexFiles :+= addCommand(new IndexDepthProfile(profileFile))
        }
        addCommand(new ComputeSampleReadDepth(profilesDir, profileIndexFiles))

        val sampleGenderReport = addCommand(new CallSampleGender)
        val plotPdf = new File(metaDataLocation, "chrY_vs_chrX.pdf")
        addCommand(new PlotChrYvsChrX(sampleGenderReport, plotPdf))
    }

    def intervalMap = {
        if (genomeIntervalList == null) {
            computeReferenceSequenceNames.map { x => x -> List(GenomeInterval.parse(x)) }.toMap
        }  else {
            expandListFiles(genomeIntervalList).map { interval => GenomeInterval.parse(interval) }
                .map { genomeInterval => genomeInterval.getSequenceName() -> genomeInterval }
                .groupBy(_._1).map {case (k,v) => (k, v.map(_._2))}
        }
    }

    class CallSampleGender() extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.CallSampleGender"
        // this.dependsOnFile = mergeOutputFiles
        this.outputFile = new File(metaDataLocation, "sample_gender.report.txt")

        commandArguments +=
            repeat(" -I ", bamLocations) +
            required(" -R ", referenceFile) +
            required(" -md ", metaDataLocation) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -genomeInterval ", genomeIntervalList) +
            optional(" -ploidyMapFile ", getReferenceMetadata.ploidyMap) +
            required(" -genderBedFile ", getReferenceMetadata.genderMaskBed) +
            optional(" -minMapQ ", depthMinimumMappingQuality)
    }

    class ComputeDepthProfile(profilesDir: File, sequenceName: String, intervalList: List[GenomeInterval]) extends JavaCommand with BAMInputOutput {
        this.javaMainClass = "org.broadinstitute.sv.apps.ComputeDepthProfiles"
        // this.dependsOnFile = mergeOutputFiles
        this.outputFile = new File(profilesDir, "profile_seq_%s_%d.dat.gz".format(sequenceName, profileBinSize))

        commandArguments +=
            repeat(" -I ", bamLocations) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            required(" -R ", referenceFile) +
            repeat(" -L ", intervalList) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -md ", metaDataLocationList) +
            repeat(" -genderMapFile ", getGenderMapFileList) +
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
            repeat(" -sequence ", intervalMap.keySet)
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
}
