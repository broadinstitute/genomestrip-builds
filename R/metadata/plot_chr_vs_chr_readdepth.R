args <- commandArgs(TRUE)

if (length(args) < 5) {
    cat("Usage: Rscript plot_chr_vs_chr_readdepth.R <dataFile> <pdfFile> <pdfTitle> <seq1> <seq2>")
    q(save="no", status=1)
}

dataFile <- args[1]
pdfFile <- args[2]
pdfTitle <- args[3]
seq1 <- args[4]
seq2 <- args[5]

data <- read.table(dataFile, row.names=1, header=TRUE, sep = "\t")

if (!(seq1 %in% colnames(data))) {
    cat(paste("Sequence", seq1, "is not in the read depth data file\n"))
    q(save="no", status=1)
}
if (!(seq2 %in% colnames(data))) {
    cat(paste("Sequence", seq2, "is not in the read depth data file\n"))
    q(save="no", status=1)
}
if (nrow(data) == 0)  {
    cat("The data file is empty")
    q(save="no", status=0)
}

pdf(file=pdfFile)

limits = range(c(0, 2, data[,seq1], data[,seq2]), na.rm=TRUE)
plot(data[, seq1], data[, seq2], xlim=limits, ylim=limits,
     pch=1, col="red", xlab=paste(seq1, " normalized read depth"), ylab=paste(seq2, " normalized read depth"), main=pdfTitle)

abline(0, 1, untf=FALSE)

dev.off()
