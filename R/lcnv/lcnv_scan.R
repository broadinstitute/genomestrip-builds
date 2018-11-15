library(reshape2)

# Experimental
USE_LOCAL_VARIANCE <- FALSE

MAX_BIN_SCORE = 1000
MAX_MEDIAN_DEV = 0.5

SCORE_THRESHOLD = 1000

cmdChrom = NULL
cmdMaxDepth = NULL
cmdTargetSamples = NULL
cmdBgSamples = NULL
cmdStartBin = NULL
cmdNumBins = NULL

cndProfileObserved = NULL
cmdProfileExpected = NULL
cmdProfileBinStarts = NULL
cmdProfileBinEnds = NULL

cmdPloidy = NULL
cmdGenderMap = NULL
cmdMinDistanceRatio = 0.5

cmdFinestRunData = NULL

parseSampleListFile <- function(sampleListFile) {
    sampleListData = read.table(sampleListFile, header=F, sep="\t")
    sampleList = as.character(sampleListData$V1)
    return(sampleList)
}

loadGenderMap <- function(dataFile) {
    if (is.null(dataFile)) {
        return(NULL)
    }
    cmdGenderMap <<- read.table(dataFile, header=TRUE, colClasses=c("character", "numeric", "numeric", "character", "character"), sep="\t", stringsAsFactors=FALSE)
}

loadPloidyMap <- function(ploidyMapFile) {
    if (is.null(ploidyMapFile)) {
        return(NULL)
    }
    ploidyMap = read.table(ploidyMapFile, sep=" ", stringsAsFactors=FALSE)
    colnames(ploidyMap) = c("CHROM", "START", "END", "GENDER", "PLOIDY")
    ploidyMap = ploidyMap[ploidyMap$CHROM != "*", ]
    ploidyMap$GENDER = ifelse(ploidyMap$GENDER == "M", "Male", "Female")
    ploidyMap$START = as.numeric(ploidyMap$START)
    ploidyMap$END = as.numeric(ploidyMap$END)
    ploidyMap
}

getPloidyIntervals <- function(ploidyMap, chrom, startPos, endPos) {

    intervals = do.call(rbind, lapply(c("Male", "Female"), function(gender) {
        genderPloidy = ploidyMap[ploidyMap$CHROM == chrom & ploidyMap$GENDER == gender, ]
        genderPloidy = genderPloidy[order(genderPloidy$START), ]
        if (nrow(genderPloidy) == 0 | nrow(genderPloidy) == 1 & genderPloidy$PLOIDY == 2) {
            genderPloidy = c(chrom, startPos, endPos, gender, 2)
        } else {
            minStart = min(genderPloidy$START)
            maxEnd = max(genderPloidy$END)
            if (minStart > startPos) {
                genderPloidy = rbind(genderPloidy, c(chrom, startPos, minStart - 1, gender, 2))
            }
            if (maxEnd < endPos) {
                genderPloidy = rbind(genderPloidy, c(chrom, maxEnd + 1, endPos, gender, 2))
            }
        }
        genderPloidy
    }))
    intervals[, c("START", "END", "PLOIDY")] = apply(intervals[, c("START", "END", "PLOIDY")], 2, function(x) {as.numeric(x)})
    intervals = intervals[with(intervals, order(GENDER, START)), ]
    intervals
}

getSamplesByGender <- function(gender) {
    # For now, we return samples with GENDER=NA together with Female samples
    if (gender == "Male") {
        selected = !is.na(cmdGenderMap$GENDER) & cmdGenderMap$GENDER == "Male"
    } else {
        selected = is.na(cmdGenderMap$GENDER) | cmdGenderMap$GENDER == "Female"
    }
    cmdGenderMap[selected, ]$SAMPLE
}

computeBinScores <- function(binRange, sample, v, runData) {
    if (is.na(v)) {
        return(rep(NA, length(binRange)))
    }
    if (v == 0) {
        #return(rep(MAX_BIN_SCORE, length(binRange)))
        return(rep(NA, length(binRange)))
    }
    obsCount = runData$observed[binRange, sample]
    expCount = runData$expected[binRange, sample]
    #mb = runData$binMedians[binRange]/2
    mb = runData$binMedians[binRange] / round(runData$binMedians[binRange])
    mb2 = mb
    if (USE_LOCAL_VARIANCE) {
        mb2 = mb * runData$binVarCoverageScale[binRange]
    }

    muAlt = expCount * mb * v
    sdAlt = sqrt(expCount * mb2 * v)
    muNull = expCount * mb * cmdPloidy
    sdNull = sqrt(expCount * mb2 * cmdPloidy)

    #cat(sprintf("muNull = %s\n", paste(muNull, collapse=" ")))
    #cat(sprintf("obs = %s\n", paste(obsCount, collapse=" ")))
    #cat(sprintf("o1/muNull = %s\n", paste(obsCount/muNull, collapse=" ")))
    testCount = ifelse(rep(v < cmdPloidy, length(binRange)), muAlt + abs(muAlt - obsCount), muAlt - abs(muAlt - obsCount))

    #cat(sprintf("muAlt = %s\n", paste(muAlt, collapse=" ")))
    #cat(sprintf("sdAlt = %s\n", paste(sdAlt, collapse=" ")))
    #cat(sprintf("muNull = %s\n", paste(muNull, collapse=" ")))
    #cat(sprintf("sdNull = %s\n", paste(sdNull, collapse=" ")))

    logDensityAlt = dnorm(testCount, mean=muAlt, sd=sdAlt, log=T)
    logDensityNull = dnorm(testCount, mean=muNull, sd=sdNull, log=T)
    #cat(sprintf("logDensityAlt = %s\n", paste(logDensityAlt, collapse=" ")))
    #cat(sprintf("logDensityNull = %s\n", paste(logDensityNull, collapse=" ")))
    return(logDensityAlt - logDensityNull)
}

computeIntervalScore <- function(startBin, endBin, sample, intervalCoverage, intervalMedianCoverage, runData) {
    #cat(sprintf("computeIntervalScore %d %d %s\n", startBin, endBin, sample))
    intervalSampleCoverage = intervalCoverage[endBin - startBin + 1, sample]
    #cat(sprintf("computeIntervalScore %d %d %.2f\n", startBin, endBin, intervalSampleCoverage))
    adjFactor = round(intervalMedianCoverage) / intervalMedianCoverage
    v = adjFactor * intervalSampleCoverage
    #v = 2 * intervalSampleCoverage / intervalMedianCoverage

    binScores = computeBinScores(startBin:endBin, sample, v, runData)
    #cat(sprintf("scores = %s\n", paste(binScores, collapse=" ")))
    return(sum(binScores, na.rm=TRUE))
}

computeIntervalCoverage <- function(startBin, endBin, sample) {
    endBin = ifelse(endBin <= cmdNumBins, endBin, cmdNumBins)

    enum = sum(cmdFinestRunData$observed[startBin:endBin, sample, drop=F])
    denom = sum(cmdFinestRunData$expected[startBin:endBin, sample, drop=F])
    
    return(enum/denom)
}

vComputeIntervalCoverage <- Vectorize(computeIntervalCoverage, c("startBin", "endBin"))

computeScoresForRange <- function(startBin, endBin, sample, runData) {
    #cat(sprintf("#DBG computeScoresForRange %d %d\n", startBin, endBin))
    cumCovNums = apply(runData$observed[startBin:endBin, sample, drop=F], 2, cumsum)
    cumCovDenoms = apply(runData$expected[startBin:endBin, sample, drop=F], 2, cumsum)
    intervalCoverage = cumCovNums / cumCovDenoms
    if (is.null(dim(intervalCoverage))) {
        dim(intervalCoverage) = c(1,1)
        colnames(intervalCoverage) = sample
    }

    # interval medians are indexed by start, length
    #intervalMedians = runData$intervalMedians[startBin, 1:nBins]

    #print(intervalCoverage)
    #print(medianCoverage)

    nBins = endBin - startBin + 1
    scores = sapply(1:nBins, function(idx) {
        eb = startBin + idx - 1
        intervalMedianCoverage = runData$intervalMedians[startBin, idx]
        computeIntervalScore(startBin, eb, sample, intervalCoverage, intervalMedianCoverage, runData)
    })
    return(scores)
}

computeScoresFor <- function(startBin, sample, runData, cpData) {
    nEndBins = length(cpData$endBinRange)
    endBin = cpData$endBinRange[nEndBins]
    #cat(sprintf("#DBG computeScoresFor %d %d\n", startBin, endBin))
    cumCovNums = apply(runData$observed[startBin:endBin, sample, drop=F], 2, cumsum)
    cumCovDenoms = apply(runData$expected[startBin:endBin, sample, drop=F], 2, cumsum)
    intervalCoverage = cumCovNums / cumCovDenoms
    if (is.null(dim(intervalCoverage))) {
        dim(intervalCoverage) = c(1,1)
        colnames(intervalCoverage) = sample
    }

    scores = sapply(cpData$endBinRange, function(eb) {
        intervalMedianCoverage = cpData$intervalMedians[startBin - cpData$startBinRange[1] + 1, eb - cpData$endBinRange[1] + 1]
        computeIntervalScore(startBin, eb, sample, intervalCoverage, intervalMedianCoverage, runData)
    })
    return(scores)
}

# Computes scores for the cross product of startBinRange x endBinRange
computeCrossProductScores <- function(startBinRange, endBinRange, targetSample, runData) {
    cpMedians = computeCrossProductMedians(startBinRange, endBinRange, cmdFinestRunData)
    cpData = list(intervalMedians=cpMedians,
                  startBinRange=startBinRange,
                  endBinRange=endBinRange)

    scores = t(sapply(startBinRange,
        function(sb) {
            scores = computeScoresFor(sb, targetSample, runData, cpData)
            return(scores)
    }))

    meltData = melt(scores, varnames=c("startBin", "endBin"), value.name="score", na.rm=T)
    meltData$startBin = meltData$startBin + startBinRange[1] - 1
    meltData$endBin = meltData$endBin + endBinRange[1] - 1

    return(meltData)
}

runIteration <- function(runData, targetSample) {
    maxDepth = runData$maxDepth
    #observed = runData$observed
    #expected = runData$expected
    #print(c(nrow(observed), ncol(observed)))
    #print(c(nrow(expected), ncol(expected)))

    nBins = length(runData$binMedians)
    scores = t(sapply(1:nBins,
        function(sb) {
            eb = pmin(sb+maxDepth-1, nBins)
            scores = computeScoresForRange(sb, eb, targetSample, runData)
            padding = sb+maxDepth-nBins-1
            if (padding > 0) {
                scores = c(scores, rep(NA, padding))
            }
            return(scores)
    }))

    meltData = melt(scores, varnames=c("start", "length"), na.rm=T)
    nMerged = runData$mergedBins
    startBins = (meltData$start-1)*nMerged + 1
    endBins = startBins + (meltData$length*nMerged) -1
    scoreData = data.frame(startBin=startBins, endBin=endBins, score=meltData$value)
    return(scoreData)
}

mergeBins <- function(valueArray, merge, calls) {
    nRows = nrow(valueArray)
    resultRows = ceiling(nRows/merge)
    result = matrix(data=0, nrow=resultRows, ncol=ncol(valueArray))
    colnames(result) = colnames(valueArray)
    for (r in 1:nRows) {
        result[ceiling(r/merge),] = result[ceiling(r/merge),] + valueArray[r,]
    }
    return(result)
}

makeCalls <- function(runData, targetSample, merge, isFullScan, scoreData, calls) {

    if (nrow(scoreData) == 0) {
        return(calls)
    }
    scoreThreshold = SCORE_THRESHOLD #* merge
    callMargin = ifelse(isFullScan, 0, 10)
    maxCallBins = cmdMaxDepth - callMargin

    df = scoreData
    df$nBins <- (df$endBin - df$startBin + 1) / merge
    df = df[df$nBins <= maxCallBins, ]

    while (nrow(df) > 0) {
        maxScore = max(df$score, na.rm=TRUE)
        if (maxScore < scoreThreshold) {
            break
        }
        call = df[df$score == maxScore, ][1, ]
        maxStartBin = call$startBin
        maxEndBin = call$endBin

        # Remove the scores for all the intervals that overlap the maxScore interval
        overlappingIntervals = !(df$endBin < maxStartBin | df$startBin > maxEndBin)
        df = df[!overlappingIntervals, ]

        numBins = (maxEndBin - maxStartBin + 1) / merge
        if (numBins >= 5) {
            overlappingIntervals = !(scoreData$endBin < maxStartBin | scoreData$startBin > maxEndBin)
            largerCalls = scoreData[overlappingIntervals & scoreData$score > maxScore, ]

            stripCallCoverage = computeIntervalCoverage(maxStartBin, maxEndBin, targetSample)
            if (nrow(largerCalls) > 1) {
                largerCalls$coverage <- vComputeIntervalCoverage(largerCalls$startBin, largerCalls$endBin, targetSample)
            }

            # Make this call if one of these is true
            # - The strip call score is the max among all the ovelapping intervals
            # - The strip call coverage is sufficiently different than the coverage for each of the overlapping intervals with the scores larger than its score
            if (nrow(largerCalls) == 0 | all(abs(largerCalls$coverage - stripCallCoverage) >= 0.2)) {
                #scoreData = scoreData[!overlappingIntervals, ]
                if (merge > 1) {
                    call = refineCall(targetSample, call, merge)
                    if (nrow(call) == 0) {
                        next
                    }
                }
                call = refineCallEdges(call, targetSample)
                isNewCall = !any(calls$startBin == call$startBin & calls$endBin == call$endBin & calls$sample == targetSample)
                if (isNewCall & call$score >= scoreThreshold) {

                    print("Adding call")
                    print(call)
                    calls = rbind(calls, data.frame(
                        startBin=call$startBin,
                        endBin=call$endBin,
                        score=call$score,
                        sample=targetSample,
                        merge=merge,
                        mSD=sd(runData$binMedians[call$startBin:call$endBin], na.rm=T),
                        stringsAsFactors=F))
                }
            }
        }
    }

    return(calls)
}

createRunData <- function(startBin, endBin, merge, calls) {
    runData = list(observed=mergeBins(cmdProfileObserved[startBin:endBin, , drop=F], merge, calls),
                   expected=mergeBins(cmdProfileExpected[startBin:endBin, , drop=F], merge, calls),
                   mergedBins=merge,
                   maxDepth=cmdMaxDepth)

    binCoverage = runData$observed[ , cmdBgSamples, drop=F] / runData$expected[ , cmdBgSamples, drop=F]
    runData$binMedians = apply(binCoverage, 1, median)
    runData$binMedians = ifelse(runData$binMedians - cmdPloidy <= MAX_MEDIAN_DEV, runData$binMedians, NA)
    runData$observed[is.na(runData$binMedians), ] = 0
    runData$expected[is.na(runData$binMedians), ] = 0
    runData$intervalMedians = computeStripIntervalMedians(runData)

    if (nrow(calls) > 0) {
        for (idx in 1:nrow(calls)) {
            call = calls[idx, ]
            callStartBin = ceiling(call$startBin / merge)
            callEndBin = ceiling(call$endBin / merge)
            runData$observed[callStartBin:callEndBin, call$sample] = 0
            runData$expected[callStartBin:callEndBin, call$sample] = 0
        }
    }

    return(runData)
}

computeIntervalMedianCoverage <- function(startBin, endBin) {
    #cat(sprintf("computeIntervalMedianCoverage: %d %d\n", startBin, endBin))
    if (endBin == startBin) {
        obs = cmdFinestRunData$observed[startBin, cmdBgSamples, drop=F]
        exp = cmdFinestRunData$expected[startBin, cmdBgSamples, drop=F]
    } else {
        obs = apply(cmdFinestRunData$observed[startBin:endBin, cmdBgSamples, drop=F], 2, sum, na.rm=T)
        exp = apply(cmdFinestRunData$expected[startBin:endBin, cmdBgSamples, drop=F], 2, sum, na.rm=T)
    }
    cov = obs / exp
    return(median(cov, na.rm=T))
}

computeStripIntervalMedians <- function(runData) {
    nBins = length(runData$binMedians)
    intervalMedians = matrix(nrow=nBins, ncol=runData$maxDepth)
    for (idx in 1:nBins) {
        obs = rep(0, length(cmdBgSamples))
        exp = rep(0, length(cmdBgSamples))
        for (jdx in 1:runData$maxDepth) {
            lastBin = idx + jdx - 1
            if (lastBin <= nBins) {
                obs = obs + runData$observed[lastBin, cmdBgSamples]
                exp = exp + runData$expected[lastBin, cmdBgSamples]
                intervalMedians[idx, jdx] = median(obs/exp)
            }
        }
    }

    return(intervalMedians)
}

computeCrossProductMedians <- function(startBinRange, endBinRange, runData) {
    nStartBins = length(startBinRange)
    nEndBins = length(endBinRange)
    intervalMedians = matrix(nrow= nStartBins, ncol=nEndBins)
    for (startBin in startBinRange) {
        obs = rep(0, length(cmdBgSamples))
        exp = rep(0, length(cmdBgSamples))
        for (endBin in startBin:endBinRange[nEndBins]) {
            obs = obs + runData$observed[endBin, cmdBgSamples]
            exp = exp + runData$expected[endBin, cmdBgSamples]
            if (endBin >= endBinRange[1]) {
                intervalMedians[startBin - startBinRange[1] + 1, endBin - endBinRange[1] + 1] = median(obs/exp)
            }
        }
    }

    return(intervalMedians)
}

# For a given call in the matrix and a margin, finds the max calls in the region constructed around this call
findRegionMaxCalls <- function(targetSample, call, margin, maxCalls = NULL) {
    cat(sprintf("%s: Finding maxCalls for %d-%d, score=%.2f...\n", Sys.time(), call$startBin, call$endBin, call$score))

    if (is.null(maxCalls)) {
        maxCalls = data.frame()
    }

    halfDistance = floor((call$endBin - call$startBin) / 2)
    innerMargin = min(halfDistance, margin)
    startBinRange = max(1, call$startBin - margin):(call$startBin + innerMargin)
    endBinRange = (call$endBin - innerMargin):min(call$endBin + margin, cmdNumBins)

    regionCalls = computeCrossProductScores(startBinRange, endBinRange, targetSample, cmdFinestRunData)

    maxScore = max(regionCalls$score)
    # Keep only the calls that do not start or end in a bin with the sample median coverage of NA
    regionCalls = regionCalls[regionCalls$score == maxScore & !is.na(cmdFinestRunData$binMedians[regionCalls$startBin]) & !is.na(cmdFinestRunData$binMedians[regionCalls$endBin]), ]

    cat(sprintf("%s: There are %d maxCalls\n", Sys.time(), nrow(regionCalls)))
    print(regionCalls)

    if (maxScore < SCORE_THRESHOLD || nrow(regionCalls) == 0) {
        return(maxCalls)
    }

    for (idx in 1:nrow(regionCalls)) {
        regionCall = regionCalls[idx, ]
        if (nrow(maxCalls[maxCalls$startBin == regionCall$startBin & maxCalls$endBin == regionCall$endBin, ]) > 0) {
            next
        }
        maxCalls = rbind(maxCalls, regionCall)
        maxCalls = findRegionMaxCalls(targetSample, regionCall, margin, maxCalls)
    }
    return(maxCalls)
}

includeBin <- function(call, binIdx, targetSample) {
    intervalMedianCoverage = computeIntervalMedianCoverage(call$startBin, call$endBin)
    sampleBinCoverage = computeIntervalCoverage(binIdx, binIdx, targetSample)
    medianBinCoverage = cmdFinestRunData$binMedians[binIdx]

    distratio = (sampleBinCoverage - medianBinCoverage)  / (call$coverage - intervalMedianCoverage)
    includeBin =
        (call$coverage > intervalMedianCoverage & sampleBinCoverage > medianBinCoverage | call$coverage < intervalMedianCoverage & sampleBinCoverage < medianBinCoverage) & distratio >= cmdMinDistanceRatio 

    includeBin = ifelse(is.na(includeBin), FALSE, includeBin)

    return(includeBin)
}

refineCallEdges <- function(call, targetSample) {
    callCoverage = computeIntervalCoverage(call$startBin, call$endBin, targetSample)
    call$coverage <- callCoverage
    startBin = call$startBin
    while (startBin > 1) {
        if (!includeBin(call, startBin - 1, targetSample)) {
            break
        }
        startBin = startBin - 1
    }

    endBin = call$endBin
    while (endBin < cmdNumBins) {
        if (!includeBin(call, endBin + 1, targetSample)) {
            break
        }
        endBin = endBin + 1
    }

    if (startBin != call$startBin | call$endBin != endBin) {
        call = computeCrossProductScores(startBin:startBin, endBin:endBin, targetSample, cmdFinestRunData)
    }
    
    return(call)
}

createCallData <- function(call) {

    startPos = cmdProfileBinStarts[call$startBin + cmdStartBin - 1]
    endPos = cmdProfileBinEnds[call$endBin + cmdStartBin - 1]
    data.frame(
        CHROM=cmdChrom,
        START=startPos,
        END=endPos,
        SAMPLE=call$sample,
        SCORE=round(call$score, 2),
        LENGTH=(endPos - startPos + 1),
        COVERAGE=round(computeIntervalCoverage(call$startBin, call$endBin, call$sample), 2),
        NBINS=(call$endBin - call$startBin + 1),
        STARTBIN=call$startBin + cmdStartBin - 1,
        ENDBIN=call$endBin + cmdStartBin - 1,
        MERGE=call$merge,
        MEDIAN_SD=call$mSD
    )
}

refineCall <- function(targetSample, call, merge) {
    cat(sprintf("%s: Refining %d-%d, score=%.2f, merge=%d ...\n", Sys.time(), call$startBin, call$endBin, call$score, merge))
    call$endBin = ifelse(call$endBin <= cmdNumBins, call$endBin, cmdNumBins)
    margin = 10
    maxCalls = findRegionMaxCalls(targetSample, call, margin)
    if (nrow(maxCalls) == 0) {
        return(maxCalls)
    }

    maxScore = max(maxCalls$score)
    maxCalls = maxCalls[maxCalls$score  == maxScore, ]

    if (nrow(maxCalls) > 1) {
        cat(sprintf("Warning: there are %d maxCalls found\n", nrow(maxCalls))) 
        maxCalls = maxCalls[1, ]
    }

    return(maxCalls)
}

main <- function() {
    cmdArguments <<- parseProgramArguments()
    if (is.null(cmdArguments$profileFile) ||
        is.null(cmdArguments$targetSample) ||
        is.null(cmdArguments$outputFile)) {
        cat("Usage: lcnv_scan.R [options]\n")
        cat("Options:\n")
        cat(" --profileFile file              Profile file(s).\n")
        cat(" --targetSample                List of samples to be scanned for LCNV events\n")
        cat(" --backgroundSample file        File containing a list of samples to be used as the background.\n")
        q(save="no", status=1)
    }

    if (length(cmdArguments$targetSample) == 1 & grepl("\\.list$", cmdArguments$targetSample[1])) {
        targetSamples <<- parseSampleListFile(cmdArguments$targetSample)
    } else {
        targetSamples <<- cmdArguments$targetSample
    }

    profileData = NULL
    for (profileFile in cmdArguments$profileFile) {
        pData = readProfileData(profileFile)
        if (is.null(profileData)) {
            profileData = pData
        } else {
            profileData = mergeProfileData(profileData, pData)
        }
    }

    # Make sure the loaded profiles contain data for all the targetSamples
    missingDataSamples = setdiff(targetSamples, colnames(profileData))
    if (length(missingDataSamples) > 0) {
        cat("Error: profile data is missing for the following samples: ", paste(missingDataSamples, collapse=","), "\n")
        q(save="no", status=1)
    }

    if (is.null(cmdArguments$backgroundSample)) {
        bgSamples = colnames(profileData)[-c(1:6)]
    } else {
        bgSamples = parseSampleListFile(cmdArguments$backgroundSample)

        if (length(setdiff(bgSamples, colnames(profileData))) > 0) {
            cat("Error: profile data is missing for some of the samples specified in --backgroundSample\n")
            q(save="no", status=1)
        }
    }

    allSamples = unique(union(targetSamples, bgSamples))

    cmdChrom <<- profileData[1, "CHR"]
    profileBinCount = nrow(profileData)

    startBin = 1
    endBin = profileBinCount
    cmdMaxDepth <<- profileBinCount
    if (!is.null(cmdArguments$startBin)) {
        startBin = as.numeric(cmdArguments$startBin)
    }
    if (!is.null(cmdArguments$endBin)) {
        endBin = as.numeric(cmdArguments$endBin)
    }
    if (endBin <= 0 | endBin > profileBinCount) {
        endBin = profileBinCount
    }
    if (!is.null(cmdArguments$maxDepth)) {
        cmdMaxDepth <<- as.numeric(cmdArguments$maxDepth)
    }
    if (cmdMaxDepth <= 0 | cmdMaxDepth > profileBinCount) {
        cmdMaxDepth <<- profileBinCount
    }

    cmdProfileBinStarts <<- profileData$START
    cmdProfileBinEnds <<- profileData$END

    profileData = profileData[ , allSamples, drop=F]
    cmdProfileObserved <<- apply(profileData, c(1,2), function(v) {
        as.numeric(strsplit(v,",")[[1]][1])
    })
    cmdProfileExpected <<- apply(profileData, c(1,2), function(v) {
        as.numeric(strsplit(v,",")[[1]][2])
    })

    ploidyMap = loadPloidyMap(cmdArguments$ploidyMapFile)
    if (is.null(ploidyMap) | is.na(match(cmdChrom, ploidyMap$CHROM))) {
        cmdPloidy <<- 2
        cmdTargetSamples <<- targetSamples
        cmdBgSamples <<- bgSamples
        calls = findLargeEvents(startBin, endBin)
    } else {
        loadGenderMap(cmdArguments$genderMapFile)
        if (is.null(cmdGenderMap)) {
            cat("The option genderMapFile must be specified to make calls on a sex chromosome\n")
            q(save="no", status=1)
        }

        ploidyIntervals = getPloidyIntervals(ploidyMap, cmdChrom, 1, cmdProfileBinEnds[profileBinCount])
        calls = do.call(rbind, lapply(1:nrow(ploidyIntervals), function(idx) {
            ploidyInterval = ploidyIntervals[idx, ]
            cmdPloidy <<- ploidyInterval$PLOIDY
            samplesByGender = getSamplesByGender(ploidyInterval$GENDER)
            cmdTargetSamples <<- intersect(targetSamples, samplesByGender)
            if (length(cmdTargetSamples) == 0) {
                return(NULL)
            }

            cmdBgSamples <<- intersect(bgSamples, samplesByGender)
            startBin = which.max(cmdProfileBinStarts >= ploidyInterval$START)
            endBin = which.max(cmdProfileBinEnds >= ploidyInterval$END | c(rep(FALSE, profileBinCount - 1), TRUE))
            cat("targetSamples: ", paste(cmdTargetSamples, collapse=" "), "\n")
            cat(paste("startBin", startBin, "\n"))
            cat(paste("endBin", endBin, "\n"))
            findLargeEvents(startBin, endBin)
        }))
    }

    if (nrow(calls) > 0) {
        write.table(calls, file = cmdArguments$outputFile, quote = FALSE, row.names = FALSE, sep="\t")
    }

    cat(sprintf("%s Done successfully.\n", Sys.time()))
    q(save="no")
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
            key = as.character(keyword)
            pos = ifelse(is.null(result[[key]]), 1, length(result[[key]]) + 1)
            result[[key]][pos] = args[[argpos]]
            argpos = argpos + 1
        }
    }
    result[[1]] = positional
    result$debug = asBoolean(result$debug)
    result$verbose = asBoolean(result$verbose)
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

readProfileData <- function(profileFile) {
    read.table(profileFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=F)
}

# Method that allows to read the data more efficiently for a subset of samples
readProfileDataSubset <- function(profileFile, samples) {
    header = unlist(read.table(profileFile, sep="\t", nrows=1, header=FALSE, stringsAsFactors=FALSE, check.names=F))

    # If this profile contains no data for the target/background samples, skip it
    if (length(intersect(header, samples)) == 0) {
        return(NULL)
    }

    colClasses = rep("NULL", length(header))
    sampleIndices = match(samples, header)
    sampleIndices = sampleIndices[!is.na(sampleIndices)]
    colClasses[1:6] <- c(rep("character", 2), rep("numeric", 3), "character")
    colClasses[sampleIndices] <- rep("character", length(sampleIndices))

    #cat(sprintf("%s Reading profile %s for %d samples...\n", Sys.time(), profileFile, length(sampleIndices)))
    profileData = read.table(profileFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=F, colClasses=colClasses)
}

mergeProfileData <- function(pData1, pData2) {
    chr1 = unique(pData1$CHR)
    chr2 = unique(pData2$CHR)
    if (chr2 != chr1) {
        cat(sprintf("Error: profiles contain different chromosomes %s and %s\n", chr1, chr2))
        q(save="no", status=1)
    }
    if (sum(pData2$START != pData1$START) > 0 | sum(pData2$END != pData1$END) > 0) {
        cat(sprintf("Error: profiles' bin start or end positions do not match\n"))
        q(save="no", status=1)
    }
    pData = merge(pData1, pData2[, -c(2:6)], by="BIN", sort=FALSE)
    pData[order(pData$START), ]
}

findLargeEvents <- function(startBin, endBin) {
    cmdStartBin <<- startBin
    cmdNumBins <<- (endBin - startBin + 1)
    endMerge = ceiling(cmdNumBins / cmdMaxDepth)
    iterations = 1 + ceiling(log2(endMerge))

    calls = data.frame()
    for (iter in 1:iterations) {
        cat(sprintf("%s: Running iteration %d of %d ...\n", Sys.time(), iter, iterations))
        merge = pmin(2**(iter-1), endMerge)
        runData = createRunData(startBin, endBin, merge, calls)
        if (iter == 1) {
            cmdFinestRunData <<- runData
        }

        isFullScan = (iter == iterations)
        for (targetSample in cmdTargetSamples) {
            cat(sprintf("%s: Computing scores for %s...\n", Sys.time(), targetSample))
            scoreData = runIteration(runData, targetSample)

            cat(sprintf("%s: Making calls for %s...\n", Sys.time(), targetSample))
            calls = makeCalls(runData, targetSample, merge, isFullScan, scoreData, calls)
            cat(sprintf("%s: Total # of calls is %d...\n", Sys.time(), nrow(calls)))
        }
    }

    callData = data.frame()
    if (nrow(calls) > 0) {
        callData = do.call(rbind, lapply(1:nrow(calls), function(idx) {
            createCallData(calls[idx, ])
        }))
    }

    return(callData)
}

main()

warnings()


