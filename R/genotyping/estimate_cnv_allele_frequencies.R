
# Globals

cmdArguments = NULL
cmdErrorCount = 0

cmdDebug = FALSE
cmdVerbose = FALSE

AF_EPSILON = 0.00001
MAX_ITERATIONS = 100
CNL_ZERO = -10000
COPY_NUMBER_CNQ_THRESHOLD = 13

# Inputs:
# Data File
#  tab-separated file with header line and three columns:
#    SAMPLE: sample ID for this sample
#    PLOIDY: 0,1,2 depending on the ploidy of this sample at this site
#    CNL: A comma-separated list of diploid copy number log10 likelihoods (starting from zero).  Shorter lists are padded (with small likelihoods) to the length of the longest one.
# Population File (optional)
#  tab-separated file with header line and two columns:
#    SAMPLE: sample ID for this sample
#    POP: population identifier
#
# Allele frequencies are estimated separately for each population.
# Populations can be overlapping (or subset/superset).
# If no population file is supplied, all samples are assumed to be in a single population.

main <- function() {
    cmdArguments <<- parseProgramArguments()
    if (is.null(cmdArguments$dataFile)) {
        cat("Usage: estimate_cnv_allele_frequencies [options]\n")
        cat("Options:\n")
        cat(" --dataFile file           File containing sample ploidy and input copy-number likelihoods [required].\n")
        cat(" --populationFile file     File containing sample population map (not necessarily unique) [optional].\n")
        cat(" --minimumCopyNumber N     Minimum haploid copy number to allow [optional, default is to calculate based on CNLs].\n")
        cat(" --maximumCopyNumber N     Maximum haploid copy number to allow [optional, default is to calculate based on CNLs].\n")
        cat(" --debug true/false        Enable debug output [optional, default false].\n")
        cat(" --verbose true/false      Enable verbose output [optional, default false].\n")
        q(save="no", status=1)
    }
    cmdDebug <<- asBoolean(cmdArguments$debug)
    cmdVerbose <<- asBoolean(cmdArguments$verbose)
    if (cmdDebug) {
      cmdVerbose <<- TRUE
    }
    dataFile = cmdArguments$dataFile
    popFile = cmdArguments$populationFile
    result = processDataFile(dataFile, popFile)
    if (cmdErrorCount == 0) {
        printResults(result)
    }
    if (cmdErrorCount > 0) {
        q(save="no", status=1)
    }
}

printResults <- function(result) {
    cat(sprintf("SUCCESS\n"))
    write.table(result, quote=FALSE, sep="\t", col.names=FALSE)
}

processDataFile <- function(dataFile, popFile) {
    data = parseDataFile(dataFile)

    ploidyData = data
    if (!all(ploidyData$PLOIDY %in% 0:2)) {
        reportFatalError("Invalid PLOIDY values found")
    }

    if (!is.null(popFile)) {
        popTable = parsePopulationFile(popFile)
    } else {
        popTable = data.frame(SAMPLE=rownames(ploidyData), POP=rep("NA", nrow(ploidyData)), stringsAsFactors=FALSE)
    }

    haploidSamples = rownames(ploidyData)[ploidyData$PLOIDY == 1]
    diploidSamples = rownames(ploidyData)[ploidyData$PLOIDY == 2]

    haploidCNLs = parseCNLs(data, haploidSamples)
    diploidCNLs = parseCNLs(data, diploidSamples)
    maxhcn = 0
    maxdcn = 0
    if (!is.null(cmdArguments$maximumCopyNumber)) {
        maxhcn = as.integer(cmdArguments$maximumCopyNumber)
    } else if (!is.null(haploidCNLs)) {
        maxhcn = ncol(haploidCNLs)-1
    }
    if (!is.null(diploidCNLs)) {
        maxdcn = ncol(diploidCNLs)-1
    }
    maxdcn = 2 * trunc(maxdcn/2 + 0.5)
    maxdcn = pmax(maxdcn, 2*maxhcn)
    maxhcn = maxdcn/2
    haploidCNLs = padCNLs(haploidCNLs, maxhcn+1)
    diploidCNLs = padCNLs(diploidCNLs, maxdcn+1)

    minhcn = cmdArguments$minimumCopyNumber
    if (is.null(minhcn)) {
        minhcn = computeMinimumHaploidCopyNumber(haploidCNLs, diploidCNLs)
        if (cmdVerbose) {
            cat(sprintf("#DBG: computed minhcn = %d\n", minhcn))
        }
    } else {
        minhcn = as.integer(minhcn)
        if (cmdVerbose) {
            cat(sprintf("#DBG: using minhcn = %d\n", minhcn))
        }
    }

    afsData = c()
    pops = sort(unique(popTable$POP))
    for (pop in pops) {
        afs = runEMForPopulation(haploidCNLs, diploidCNLs, minhcn, ploidyData, popTable, pop)
        afsData = c(afsData, afs)
    }
    afsMatrix = matrix(afsData, nrow=length(pops), byrow=TRUE)
    rownames(afsMatrix) = pops
    return(afsMatrix)
}

computeMinimumHaploidCopyNumber <- function(haploidCNLs, diploidCNLs) {
    if (cmdVerbose) {
        cat(sprintf("#DBG: Computing mininum hcn with %d+%d samples...\n", nrow(haploidCNLs), nrow(diploidCNLs)))
    }
    minhcn = NA
    cns = apply(haploidCNLs, 1, callCopyNumber, COPY_NUMBER_CNQ_THRESHOLD)
    if (!all(is.na(cns))) {
        minFromHaploid = min(cns, na.rm=T)
        minhcn = minFromHaploid
    }
    cns = apply(diploidCNLs, 1, callCopyNumber, COPY_NUMBER_CNQ_THRESHOLD)
    if (!all(is.na(cns))) {
        minFromDiploid = floor(min(cns, na.rm=T)/2)
        if (is.na(minhcn) || minFromDiploid < minhcn) {
            minhcn = minFromDiploid
        }
    }
    if (is.na(minhcn)) {
        minhcn = 0
    }
    return(minhcn)
}

callCopyNumber <- function(cnls, qualThreshold) {
    maxIndex = which.max(cnls)
    if (length(maxIndex) == 0) {
        return(NA)
    }
    secondIndex = which.max(cnls[-maxIndex])
    if (length(secondIndex) == 0) {
        return(maxIndex - 1)
    }
    secondIndex = ifelse(secondIndex < maxIndex, secondIndex, secondIndex+1)
    qual = 10*(cnls[maxIndex]-cnls[secondIndex])
    if (qual >= qualThreshold) {
        return(maxIndex - 1)
    }
    return(NA)
}

runEMForPopulation <- function(haploidCNLs, diploidCNLs, minhcn, ploidyData, popTable, pop) {
    popSamples = popTable$SAMPLE[popTable$POP == pop]
    samples = intersect(rownames(ploidyData)[ploidyData$PLOIDY > 0], popSamples)
    if (length(samples) == 0) {
        nAlleles = ncol(haploidCNLs)
        return(c(NSAMPLES=0,NCHROMOSOMES=0,rep(NA, nAlleles)))
    }

    haploidSamples = intersect(rownames(haploidCNLs), popSamples)
    diploidSamples = intersect(rownames(diploidCNLs), popSamples)
    nHaploidSamples = length(haploidSamples)
    nDiploidSamples = length(diploidSamples)
    haploidCNLs = subsetCNLs(haploidCNLs, haploidSamples)
    diploidCNLs = subsetCNLs(diploidCNLs, diploidSamples)

    maxhcn = ncol(haploidCNLs)-1
    maxdcn = ncol(diploidCNLs)-1
    nAlleles = maxhcn+1
    nDCN = maxdcn+1

    if (cmdVerbose) {
        cat(sprintf("#DBG: Running EM for population %s with %d+%d samples...\n", pop, nHaploidSamples, nDiploidSamples))
    }
    if (cmdDebug) {
        cat(sprintf("#DBG: maxhcn = %d\n", maxhcn))
        cat(sprintf("#DBG: maxdcn = %d\n", maxdcn))
        cat(sprintf("#DBG: nSamples = %d+%d\n", nHaploidSamples, nDiploidSamples))
        cat(sprintf("#DBG: nAlleles = %d\n", nAlleles))
        cat(sprintf("#DBG: nDCN = %d\n", nDCN))
    }

    if (nHaploidSamples + nDiploidSamples == 0) {
        afs = rep(NA, nAlleles)
        return(afs)
    }

    if (nHaploidSamples == 0) {
        obsHaploid = rep(0, nAlleles)
    } else {
        probMatrix = t(apply(haploidCNLs, 1, function(cnl) { llToProbs(cnl) }))
        obsHaploid = apply(probMatrix, 2, sum)
        if (cmdDebug) {
            cat(sprintf("#DBG: obsHaploid:\n"))
            print(obsHaploid)
            print(as.integer(obsHaploid))
        }
    }

    if (nDiploidSamples == 0) {
        afs = obsHaploid/sum(obsHaploid)
        names(afs) = paste("CN", (1:length(afs))-1, sep="")
        result = c(NSAMPLES=length(samples),
                   NCHROMOSOMES=nHaploidSamples + 2*nDiploidSamples,
                   afs)
        return(result)
    }

    probMatrix = t(apply(diploidCNLs, 1, function(cnl) { llToProbs(cnl) }))
    obsDiploid = apply(probMatrix, 2, sum)
    if (cmdDebug) {
        cat(sprintf("#DBG: obsDiploid:\n"))
        print(obsDiploid)
        print(as.integer(obsDiploid))
    }

    dcnMap = outer(0:maxhcn, 0:maxhcn, FUN="+")
    obsDiploidMatrix = apply(dcnMap, 1:2, function(dcn) { obsDiploid[dcn+1] })
    if (cmdDebug) {
        cat(sprintf("#DBG: dcnMap:\n"))
        print(dcnMap)
        cat(sprintf("#DBG: obsDiploidMatrix:\n"))
        print(obsDiploidMatrix)
    }

    # initial conditions
    # using a uniform distribution on the initial conditions leads to instability and bogus estimates in some cases
    # instead we use an initial afs estimate based on the diagonal of the obsDiplolidMatrix, essentially the ratio of likely homozygotes.
    # to this we add a smoothing factor of 1/nAlleles to avoid initial estimates of zero
    # flatafs = rep(1/nAlleles, nAlleles)
    # flatafs = c(rep(0,minhcn), rep(1/(nAlleles-minhcn), nAlleles-minhcn))
    obsEstimate = sapply(1:nAlleles, function(hcn) { obsDiploidMatrix[hcn,hcn] })
    if (minhcn > 0) {
        obsEstimate[1:minhcn] = 0
    }
    smoothingFactor = c(rep(0,minhcn), rep(1/(nAlleles-minhcn), nAlleles-minhcn))
    obsEstimate = obsEstimate + smoothingFactor
    afs = obsEstimate/sum(obsEstimate)
    if (nHaploidSamples > 0) {
        obsTotal = obsHaploid + 2 * afs * nDiploidSamples
        if (minhcn > 0) {
            obsTotal[1:minhcn] = 0
        }
        afs = obsTotal/sum(obsTotal)
    }
    if (cmdDebug) {
        cat(sprintf("#DBG: iteration 0, delta NA, afs = %s\n", formatAFS(afs)))
    }

    for (iteration in 1:MAX_ITERATIONS) {
        # estimation step
        afMatrix = outer(afs, afs, FUN="*")
        diagonalSums = sapply(0:maxdcn, function(dcn) { sum(afMatrix[dcnMap == dcn]) })
        if (cmdDebug) {
            cat(sprintf("#DBG: diagonalSums:\n"))
            print(diagonalSums)
        }

        dsumMatrix = apply(dcnMap, 1:2, function(dcn) { diagonalSums[dcn+1] })
        if (cmdDebug) {
            cat(sprintf("#DBG: dsumMatrix (sum = %f):\n", sum(dsumMatrix)))
            print(dsumMatrix)
        }

        proportions = afMatrix/dsumMatrix
        proportions[is.nan(proportions)] = 0
        if (cmdDebug) {
            cat(sprintf("#DBG: proportions:\n"))
            print(proportions)
        }

        estimates = proportions * obsDiploidMatrix
        if (cmdDebug) {
            cat(sprintf("#DBG: estimates:\n"))
            print(estimates)
        }

        # maximization step
        newafs = sapply(1:nAlleles, function(n) { sum(estimates[n,]) + sum(estimates[,n]) })
        newafs = newafs + obsHaploid
        if (cmdDebug) {
            cat(sprintf("#DBG: sum(newafs) = %f\n", sum(newafs)))
        }
        newafs = newafs/sum(newafs)
        delta = max(abs(newafs - afs))
        afs = newafs
        if (cmdDebug) {
            cat(sprintf("#DBG: iteration %d, delta %f, afs = %s\n", iteration, delta, formatAFS(afs)))
        }
        if (delta < AF_EPSILON) {
            break
        }
    }
    if (cmdVerbose) {
        cat(sprintf("#DBG: EM for population %s: afs = %s\n", pop, formatAFS(afs)))
    }

    names(afs) = paste("CN", (1:length(afs))-1, sep="")
    result = c(NSAMPLES=length(samples),
               NCHROMOSOMES=nHaploidSamples + 2*nDiploidSamples,
               afs)
    return(result)
}

llToProbs <- function(lls) {
    temp = exp(lls-max(lls))
    temp[is.na(lls)] = 0
    return(temp/sum(temp))
}

formatAFS <- function(afs) {
    return(paste(sprintf("%1.5f", afs), collapse=" "))
}

parseCNLs <- function(cnlData, samples) {
    if (length(samples) == 0) {
        return(matrix(nrow=0, ncol=0))
    }
    cnls = sapply(strsplit(cnlData[samples,]$CNL, ",", fixed=TRUE), as.numeric, simplify=FALSE)
    maxlength = max(sapply(cnls, length))
    cnls = t(sapply(cnls, function(row) { c(row, rep(CNL_ZERO, maxlength - length(row))) }))
    cnls[is.na(cnls)] = CNL_ZERO  # handle missing input data
    rownames(cnls) = samples
    return(cnls)
}

padCNLs <- function(cnls, ncols) {
    while (ncol(cnls) < ncols) {
        cnls = cbind(cnls, rep(CNL_ZERO, nrow(cnls)))
    }
    return(cnls)
}

subsetCNLs <- function(cnls, samples) {
    result = cnls[samples,]
    if (length(samples) == 1) {
        result = matrix(result, nrow=1)
        rownames(result) = samples
    }
    return(result)
}

parseDataFile <- function(dataFile) {
    return(read.table(dataFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1))
}

parsePopulationFile <- function(popFile) {
    return(read.table(popFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, na.strings=c()))
}

reportFatalError <- function(message) {
    cat(sprintf("ERROR\t%s\n", message, file=stderr()))
    q(save="no", status=1)
}

parseProgramArguments <- function() {
    result = list()
    positional = list()
    result[[""]] = positional
    args = commandArgs()
    if (length(args) == 0) {
        return(result)
    }
    for (i in 1:length(args)) {
        if (i == length(args)) {
            return(result)
        } else if (args[i] == "--args") {
            argpos = i+1
            break
        }
    }
    while (argpos <= length(args)) {
        arg = args[argpos]
        argpos = argpos + 1
        keyword = NULL
        if (nchar(arg) > 2 && substr(arg,1,2) == "--") {
            keyword = substr(arg,3,nchar(arg))
        } else if (nchar(arg) > 1 && substr(arg,1,1) == "-") {
            keyword = substr(arg,2,nchar(arg))
        } else {
            positional = c(positional, arg)
        }
        #cat(sprintf("pos %d kw %s arg %s\n", argpos, keyword, args[argpos]))
        if (!is.null(keyword) && argpos <= length(args)) {
            keyword = as.character(keyword)
            arg = as.character(args[[argpos]])
            argpos = argpos + 1
            result[[keyword]] = c(result[[keyword]], arg)
        }
    }
    result[[1]] = positional
    return(result)
}

asBoolean <- function(arg) {
    if (is.null(arg)) {
        return(FALSE)
    }
    if (is.na(arg)) {
        return(FALSE)
    }
    if (is.logical(arg)) {
        return(arg)
    }
    return(arg == "true")
}

main()
