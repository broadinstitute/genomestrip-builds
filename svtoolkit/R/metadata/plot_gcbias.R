library(reshape2)
library(ggplot2)
library(directlabels)

printUsage <- function(args) {
    cat("Usage: Rscript plotGCBias.R outfile_prefix referenceProfile -project project_name1 project_file1 ... [-project project_namen project_filen]\n" )
    stop(paste("invalid arguments: ", args))
}

isZipFile <- function(fileName) {
    substring(fileName, nchar(fileName)-3) %in% c(".zip", ".ZIP")
}

readGroup <- function (groupName, fileName) {
    GCBiasGroup <- read.table(fileName, header=T, sep="\t")
    GCBiasGroup$group <- groupName
    return(GCBiasGroup)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    printUsage(args)
}

outFile <- args[1]
referenceFile <- args[2]
restOfArgs <- args[3:length(args)]

if (isZipFile(referenceFile)) {
    referenceProfile <- read.table(unz(referenceFile, "gcprofile.200.txt"), skip=1)
} else {
    referenceProfile <- read.table(referenceFile, skip=1)
}
referenceProfile$bin = 1:length(referenceProfile)

annotateOutliers <- FALSE

allGCBias <- data.frame()
while (length(restOfArgs) > 0) {
  if (restOfArgs[1] == "-project") {
    groupName <- restOfArgs[2]
    fileName <- restOfArgs[3]
    cat(sprintf("Loading project %s ...\n", groupName))
    allGCBias <- rbind(allGCBias, readGroup(groupName, fileName))
    if (length(restOfArgs) > 3) {
      restOfArgs <- restOfArgs[4:length(restOfArgs)]
    } else {
      restOfArgs <- c()
    }
  } else if (restOfArgs[1] == "-annotateOutliers") {
    annotateOutliers <- TRUE
    if (length(restOfArgs) > 1) {
      restOfArgs <- restOfArgs[2:length(restOfArgs)]
    } else {
      restOfArgs <- c()
    }
  } else {
    printUsage(args)
  }
}

cat(sprintf("Done loading data.\n"))

allData <- melt(allGCBias, id=c("SAMPLE", "LIBRARY", "READGROUP", "NBINS", "group"))
allData$SAMPLE_LIBRARY <- paste(allData$SAMPLE, allData$LIBRARY, sep="/")

plotDetail <- function(outFile, allData) {
  pdf(outFile, height=10, width=12)
  p <- ggplot(allData, aes(x=as.numeric(variable)/max(as.numeric(variable)), y=value, group=SAMPLE, color=group)) + 
    xlab("GC Fraction") + geom_line() + 
    geom_dl(aes(label=SAMPLE), method=dl.combine("first.qp", "last.qp"))
  print(p)
  dev.off()
}

allData$alphaVal <- .99

plotLarge <- function(outFile, allData, referenceProfile, limitY=TRUE) {
  allData[allData$SAMPLE_LIBRARY != allData[1,"SAMPLE_LIBRARY"],]$windowPortions <- NA
  allData$group2 <- factor(allData$group, as.character(unique(allData$group)))
  maxY <- ifelse(limitY, 2.5, max(2.5, max(allData$value, na.rm = TRUE) + .5))
  p <- ggplot(allData, aes(x=gcValues, y=value, group=SAMPLE_LIBRARY, colour=group2)) + 
    theme_grey(base_size = 18) +
    geom_line(aes(alpha=alphaVal)) + 
    geom_text(aes(label=ANNOTATION), colour="black", hjust=0, vjust=0, size=4) +
    scale_y_continuous(limits=c(-0.51,maxY), breaks=seq(from=0, to=maxY, by=.5)) + xlim(.2,.8) +
    geom_linerange(aes(x=gcValues, ymin=-.51, ymax=(-.51 + windowPortions)), size=.5, colour="black") + 
    geom_hline(aes(x=0)) +  
    guides(size=FALSE, alpha=FALSE) + labs(colour="") + 
    ylab("Normalized Coverage (fragments)") + 
    xlab("GC% of 400bp windows") + # scale_colour_manual(values=scale_colour_hue()$palette(3)[2:3]) +
    annotate("text", x=.7, y=-.25, label="Portion of Reference Genome", colour="black") 
  
  pdf(paste0(outFile, ".pdf"), width=12, height=8)
  print(p)
  dev.off()
  
  #png(paste0(outFile, ".png"), width=1200, height=800)
  #print(p)
  #dev.off()
}

allData$windowPortions <- 0.5 * referenceProfile[as.numeric(allData$variable), "V3"] / max(referenceProfile$V3)  

allData$gcValues <- as.numeric(allData$variable)/max(as.numeric(allData$variable))

allData$correctedValue <- ifelse(is.na(allData$value), 0, allData$value)

allData$ANNOTATION <- ""
if (annotateOutliers) {
  for (sample in unique(allData$SAMPLE_LIBRARY)) {
    x <- allData[allData$SAMPLE_LIBRARY == sample & allData$gcValues > .2 & allData$gcValues < .8,]
    allData[allData$SAMPLE_LIBRARY == sample, "ANNOTATION"] <- NA
  
    maxRow <- x[x$correctedValue == max(x$correctedValue),]
    binVals <- allData[allData$variable == maxRow$variable & allData$group == maxRow$group,"correctedValue"]

    if (maxRow$correctedValue > median(binVals) + 5 * mad(binVals)) {
      allData[allData$SAMPLE_LIBRARY == maxRow$SAMPLE_LIBRARY & allData$variable == maxRow$variable, "ANNOTATION"] <- maxRow$SAMPLE_LIBRARY
    }
  
    minRow <- x[x$correctedValue == min(x$correctedValue),]
    binVals <- allData[allData$variable == minRow$variable & allData$group == minRow$group,"correctedValue"]
    if (minRow$correctedValue < median(binVals) - 5 * mad(binVals)) {
      allData[allData$SAMPLE_LIBRARY == minRow$SAMPLE_LIBRARY & allData$variable == minRow$variable, "ANNOTATION"] <- minRow$SAMPLE_LIBRARY
    }
  
  }
}

plotLarge(outFile, allData, referenceProfile)
plotLarge(paste0(outFile,"_noMaxY"), allData, referenceProfile, limitY = FALSE)

