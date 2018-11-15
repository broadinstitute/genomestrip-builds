computeLogPValue <- function(row) {
    mean = row["MEAN"]
    variance = row["VARIANCE"]
    x = row["COUNT"]

    lowerTail = (x <= mean)
    return(log(2) + pnorm(x, mean=mean, sd=sqrt(variance), lower.tail = lowerTail, log.p = TRUE))
}

args = commandArgs(TRUE)
dataFile = args[1]
if (is.na(dataFile)) {
    cat("Usage: compute_cn_qual.R <dataFile>")
    q(save="no", status=1)
}

data <- read.table(dataFile, sep="\t", header=TRUE, colClasses=c("numeric", "numeric", "numeric", "logical"))

data$LOG_P_VALUE <- apply(data, 1, computeLogPValue)
varData = data[data$IS_VARIANT, ]

LOG10_E = log10(exp(1))
minScore = -LOG10_E * min(data$LOG_P_VALUE)

logSumAll = -2 * sum(data$LOG_P_VALUE)
fisherAllScore = pchisq(logSumAll, df = 2 * nrow(data), lower.tail = FALSE, log.p = TRUE)
fisherAllScore = -LOG10_E * fisherAllScore

fisherVarScore = NA
if (nrow(varData) > 0) {
    logSumVar = -2 * sum(varData$LOG_P_VALUE)
    fisherVarScore = pchisq(logSumVar, df = 2 * nrow(varData), lower.tail = FALSE, log.p = TRUE)
    fisherVarScore = -LOG10_E * fisherVarScore
}

cat(sprintf("%f\t%f\t%f\t%f\t%f\n", minScore, fisherAllScore, fisherVarScore, nrow(data), logSumAll))
