/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */

import java.io.File

import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.sv.qscript.SVQScript

class UnarchiveData extends SVQScript {

    @Argument(shortName="archive", required=true, doc="File to be unarchived")
    var archivesList: List[File] = Nil

    @Argument(shortName="deleteArchives", required=false, doc="Delete the archive(s)")
    var deleteArchivesArg: String = null
    var deleteArchives: Boolean = false
	    
    @Output(shortName="O", fullName="outputDir", required=true, doc="The output directory.")
    var outputDir: File = null
	
    def script = {
      if (deleteArchivesArg != null) {
        deleteArchives = parseBoolean(deleteArchivesArg)
      }
      createDirectory(outputDir)
      var idx: Int = 0
      for (archive <- expandFileListFiles(archivesList)) {
        addCommand(new UnarchiveCommand(idx, new File(archive), outputDir))
        idx += 1
      }
    }

    class UnarchiveCommand(idx: Int, archive: File, outputDir: File) extends CommandLineFunction
        with FileInputOutput {

        this.jobOutputFile = new File(jobLogDirectory, "Unarchive_" + idx + ".out")
        this.outputFile = jobOutputFile

        this.jobArrayName = "Unarchive"

        var archiveOutputDir = new File(outputDir.getPath, "archive_" + idx)
        createDirectory(archiveOutputDir)

        val unarchiveCmd = archive.getName() match {
            case name if name.endsWith(".zip") =>
                "unzip -d " + archiveOutputDir + " " + archive
            case name if name.endsWith(".tar.gz") => 
                " tar -xvzf " + archive +
                " -C " + archiveOutputDir.getPath +
                " --strip-components 1"
            case _ => throw new RuntimeException("Unsupported archive type: " + archive)
        }

        val deleteCmd = if (deleteArchives) " &&  rm " + archive else ""
        override def commandLine =
            unarchiveCmd +
            deleteCmd
    }
}
