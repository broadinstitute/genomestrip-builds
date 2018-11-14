args = commandArgs(TRUE)
vpsReport = args[1]
madRange = as.numeric(args[2])
samplesFile = args[3]

data = read.table(vpsReport, header=T, stringsAsFactors=F, sep="\t")

discoverySamples = data[data$VARIANTS <= median(data$VARIANTS) + madRange * mad(data$VARIANTS), "SAMPLE"]

write.table(discoverySamples, file=samplesFile, quote=F, col.names=F, row.names=F, sep="\t")
