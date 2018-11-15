
# Globals

cmdArguments = NULL
cmdErrorCount = 0

site.data <- NULL
coverage.data <- NULL
genotype.data <- NULL
genotype.qual.data <- NULL
pair.data <- NULL
truth.data <- NULL
sample.group.data <- NULL

main <- function() {
    options(warn = 2)
    cmdArguments <<- parseProgramArguments()
    siteList = if (is.null(cmdArguments$site)) NULL else unlist(strsplit(cmdArguments$site, ","))
    outputFile = if (is.null(cmdArguments$outputFile)) cmdArguments$O else cmdArguments$outputFile
    if (is.null(cmdArguments$dataFile) || is.null(cmdArguments$site) || is.null(outputFile)) {
        cat("Usage: plot_cnps_from_vcf [options]\n")
        cat("Options:\n")
        cat(" --dataFile file    File containing the genotyping data [required].\n")
        cat(" --site siteList            List of site IDs to plot [required].\n")
        cat(" --outputFile file          Output file name [required].\n")
        cat(" --outputFileFormat format  Output file type (PDF or PNG) [optional, default PDF].\n")
        cat(" --genderMapFile file       Map file providing the gender of each sample [optional].\n")
        cat(" --truthDataFile file       File containing truth data genotypes [optional].\n")
        cat(" --sampleFile file     File containing samples to be plotted [optional].\n")
        cat(" --sampleGroupFile file     File containing samples/groups to highlight [optional].\n")
        cat(" --excludedSiteSampleFile file     File containing site/samples mapping specifying what sample(s) are not to be plotted for the given site [optional].\n")
        cat(" --debug true/false         Enable debug output [optional, default false].\n")
        cat(" --verbose true/false       Enable verbose output [optional, default false].\n")
        cat(" --pretty true/false        Produce prettier plots but with less detail [optional, default false].\n")
        q(save="no", status=1)
    }
    plotCnps(siteList, outputFile)
    if (cmdErrorCount > 0) {
        q(save="no", status=1)
    }
}

plotCnps <- function(siteList, outputFile) {
    outputFormat = cmdArguments$outputFileFormat
    if (is.null(outputFormat) || outputFormat == "PDF") {
        if (length(siteList) == 1) {
            pdf(outputFile, width=8, height=5)
            on.exit(dev.off(),add=T)
        } else {
            pdf(outputFile, width=8, height=10)
            on.exit(dev.off(),add=T)
            layout(matrix(1:2,2,1))
            on.exit(layout(1),add=T)
        }
    } else if (outputFormat == "PNG") {
        if (length(siteList) != 1) {
            reportFatalError("Output type PNG can only be used with a single site")
        }
        png(outputFile, type="cairo", width=8, height=5, units="in", res=300)
        on.exit(dev.off(),add=T)
    } else {
        reportFatalError(paste("Unrecognized outputFileFormat:", outputFormat, collapse=" "))
    }
    plotSites(siteList)
}

plotSites <- function(siteList, threshold=13, ymax=NULL, impMethod=NULL) {
    loadDataFile()
    gender.data <- loadGenderMap()
    sample.data <- loadSamples()
    excluded.site.samples <- loadExcludedSiteSamples()
    truth.data <<- loadTruthData()
    sample.group.data <<- loadSampleGroups()

    for (cnp in siteList) {
        if (!(cnp %in% rownames(genotype.data))) {
            #cmdErrorCount <<- cmdErrorCount + 1
            cat(sprintf("Warning: No data found for site %s\n", cnp))
            next
        }
        if (asBoolean(cmdArguments$verbose)) {
            cat(sprintf("Plotting site %s ...\n", cnp))
        }
        plotCnpInternal(cnp,
                        site.data[cnp,],
                        coverage.data[cnp,],
                        genotype.data[cnp,],
                        genotype.qual.data[cnp,],
                        pair.data[cnp,],
                        truth.data[cnp,],
                        sample.group.data,
                        sample.data,
                        excluded.site.samples[excluded.site.samples$ID == cnp, ]$SAMPLE,
                        gender.data,
                        threshold=threshold,
                        ymax=ymax)
    }
}

plotCnpInternal <- function(cnp, site.data, coverage.data, genotype.data, genotype.qual.data, pair.data, truth.data, sample.group.data=NULL, sample.data=NULL, excluded.site.samples=NULL, gender.data=NULL, threshold=NULL, ymax=NULL) {
    # genotype.data is *our* genotype calls
    # optionally, we can take in gold standard reference calls and display those as well
    # make sure we plot cna genotypes (e.g. light blue)

    model = parseModelData(site.data)

    pretty <- asBoolean(cmdArguments$pretty)
    plot.legend <- TRUE
    if (pretty) {
        plot.legend <- FALSE
    }
    break.width <- 0.05
    cn.colors <- c("green", "orange", "blue", "red", "pink", "yellow", "purple")
    if (pretty) {
        cn.colors <- c("green", "orange", "blue", "red", "pink", "yellow", "purple", "darkgreen", "lightblue", "yellow3", "red3", "steelblue", "green3", "pink3", "brown", "blue3")
    }
    cna.color <- "lightgray"

    samples <- getSamples(genotype.data, sample.data, excluded.site.samples)

    covdata <- asNumericInf(coverage.data[samples])
    names(covdata) <- samples
    if (!pretty && !is.null(model)) {
        m1 = model$MEANS[2]
        covdata = m1 * covdata
    }

    gtdata <- as.numeric(genotype.data[samples])
    names(gtdata) <- samples

    if (!is.null(threshold) && !is.null(genotype.qual.data)) {
        confdata <- as.numeric(t(genotype.qual.data[samples]))
        gtdata[confdata < threshold] <- NA
    }

    cnmax <- max(c(2,gtdata), na.rm=T)
    dmax <- ceiling(max(c(0,covdata[is.finite(covdata)]))+0.5)
    breaks <- seq(0,dmax,break.width)
    nbreaks <- length(breaks)

    hdata <- c()
    for (cn in 0:cnmax) {
        cnx <- covdata[!is.na(gtdata) & gtdata == cn]
        hn <- makeHist(seq(0,nbreaks,1), 0, cnx, round(cnx/break.width), length)
        hdata <- c(hdata, hn, 0)
    }
    cna <- covdata[is.na(gtdata)]
    names(cna) <- samples[is.na(gtdata)]
    hn <- makeHist(seq(0,nbreaks,1), 0, cna, round(cna/break.width), length)
    hdata <- c(hdata, hn, 0)
    mat <- matrix(hdata, nrow=cnmax+2, byrow=TRUE)
    mat.height <- max(apply(mat,2,sum))

    cnp.length <- max(0, site.data$END - site.data$START + 1)

    ncalls <- sum(!is.na(gtdata))
    call.rate <- ncalls / length(gtdata)
    ncalls <- sum(mat[1:3,])
    nnonref <- 2*sum(mat[1,]) + sum(mat[2,])
    if (cnmax >= 3) {
        ncalls <- ncalls + sum(mat[4,])
        nnonref <- nnonref + sum(mat[4,])
        if (cnmax >= 4) {
            ncalls <- ncalls + sum(mat[5,])
            nnonref <- nnonref + 2*sum(mat[5,])
        }
    }
    maf <- ifelse(ncalls == 0, NA, nnonref/(2*ncalls))
    eff.length <- site.data$ELENGTH
    cn.qual = as.numeric(site.data$CNQUAL)

    title1 <- sprintf("%s", cnp)
    if (!is.null(truth.data) && ("ALT" %in% names(truth.data))) {
        title1 <- sprintf("%s / %s", title1, truth.data$ALT)
    }
    if (substr(site.data$CHR,1,3) == "chr") {
        siteInterval <- sprintf("%s %d-%d", site.data$CHR, site.data$START, site.data$END)
    } else {
        siteInterval <- sprintf("chr%s %d-%d", site.data$CHR, site.data$START, site.data$END)
    }
    title2 <- siteInterval
    if (cnp.length >= 1000) {
        title2 <- sprintf("%s %1.1fKb", title2, cnp.length/1000)
    } else {
        title2 <- sprintf("%s %db", title2, cnp.length)
    }
    #bkpt.id <- get_breakpoint_id(cnp)
    #if (!is.null(bkpt.id)) {
    #    title2 <- sprintf("%s %s", title2, bkpt.id)
    #}
    title3 <- sprintf("CR: %1.1f%%", call.rate * 100)
    if (!is.na(cn.qual)) {
        title3 <- sprintf("%s CNQUAL: %1.2f", title3, cn.qual)
    }
    if (!is.null(truth.data)) {
        truth.samples <- intersect(samples, names(truth.data))
        truth <- truth.data[truth.samples]
        gtcalls <- gtdata[truth.samples]
        miscalls <- sum(!is.na(gtcalls) & !is.na(truth) & gtcalls != truth)
        denom <- sum(!is.na(gtcalls) & !is.na(truth))
        if (denom > 0) {
            accuracy <- 1 - (miscalls / denom)
            title3 <- sprintf("%s ACC: %1.1f%%", title3, accuracy * 100)
        } else {
            title3 <- sprintf("%s ACC: NA", title3)
        }
        if (sum(!is.na(gtcalls) & !is.na(truth)) > 0) {
            denom2 <- sum(!is.na(gtcalls) & !is.na(truth) & (gtcalls != 2 | truth != 2))
            discord.rate <- ifelse(denom2 == 0, 0, miscalls/denom2)
            title3 <- sprintf("%s DR: %1.1f%%", title3, discord.rate * 100)
        } else {
            title3 <- sprintf("%s DR: NA", title3)
        }
    }
    title4 <- sprintf("MAF: %1.2f", maf)
    if (!is.null(eff.length)) {
        cnp.conf.length <- max(0, site.data$RIGHTSTART - site.data$LEFTEND - 1)
        eff.fraction <- eff.length / cnp.conf.length
        if (eff.length > 1000) {
            title4 <- sprintf("%s EL: %1.1fKb", title4, eff.length/1000)
        } else {
            title4 <- sprintf("%s EL: %db", title4, eff.length)
         }
         title4 <- sprintf("%s %1.1f%%", title4, eff.fraction*100)
    }

    if (is.null(ymax)) {
        ymax <- mat.height
    }
    bar.colors <- c(cn.colors[1:pmin(cnmax+1,length(cn.colors))])
    if (length(bar.colors) < cnmax + 1) {
        bar.colors <- c(bar.colors, rep(bar.colors[length(bar.colors)], cnmax + 1 - length(bar.colors)))
    }
    bar.colors <- c(bar.colors, cna.color)
    ### hack to show only truth data
    ### bar.colors <- rep(cna.color, length(bar.colors))
    titles <- c(title1, title2, title3, title4)
    if (pretty) {
        titles <- c(title1, title2)
    }
    barplot(mat, space=0, col=bar.colors,
            ylim=c(0,pmin(pmax(1,mat.height),ymax)),
            xlab="normalized read depth", ylab="samples",
            main=titles)

    plot.x.scale = 1/break.width

    dseq = seq(0, dmax, 1)
    atPos = 0.5 + plot.x.scale * dseq
    axis(1, at=atPos, labels=dseq)

    if (plot.legend) {
        legend("topright",
               fill=c(cn.colors,cna.color),
               legend=c(sprintf("CN%d", 0:(length(cn.colors)-2)),
                        sprintf("CN%d+", length(cn.colors)-1),
                        "NC"),
               inset=0.02, cex=0.8)
    }

    if (mat.height == 0) {
        return(NULL)
    }

    # build position map for samples in bar chart
    # order all samples by depth and genotype call and determine correct bin and offset
    # print(mat)
    colsums <- apply(mat,2,sum,na.rm=T)
    cumsums <- cumsum(colsums)
    b <- c(mapply(function(n) { rep(n, colsums[n]) }, 1:length(colsums)), recursive=T)
    s <- names(covdata)[order(covdata,na.last=T)]
    s <- setdiff(s,names(covdata)[is.na(covdata)])
    g <- gtdata[s]
    s <- s[order(b,g,na.last=T)]
    g <- gtdata[s]
    d <- covdata[s]
    offsets <- 1:length(b) - c(0,head(cumsums,-1))[b]
    bin.data <- data.frame(d=d,g=g,b=b,off=offsets)
    rownames(bin.data) <- s

    # print(d)
    # print(b)
    # print(l)
    # print(cumsums)
    # print(offsets)
    # print(length(d))
    # print(length(l))
    # print(length(b))
    # print(length(offsets))

    if (!is.null(truth.data)) {
        truth.colors <- sapply(rownames(bin.data),
                               function(sample) { v <- truth.data[cnp,sample];
                                                  if (is.null(v) || is.na(v)) {
                                                      return(NA)
                                                      #return(cna.color)
                                                  } else if (v > length(cn.colors)) {
                                                      return(cn.colors[length(cn.colors)])
                                                  } else {
                                                      return(cn.colors[v+1])
                                                  }})
        xvals <- bin.data$b[!is.null(truth.colors)]-0.5
        yvals <- bin.data$off[!is.null(truth.colors)]-0.3
        col <- truth.colors[!is.null(truth.colors)]
        points(xvals, yvals, pch=21, col=col, bg=col)
    }
    if (!is.null(sample.group.data)) {
        sample.colors <- sample.group.data[rownames(bin.data),]$COLOR
        xvals <- bin.data$b[!is.null(sample.colors)]-0.5
        yvals <- bin.data$off[!is.null(sample.colors)]-0.3
        col <- sample.colors[!is.null(sample.colors)]
        points(xvals, yvals, pch=21, col=col, bg=col)
    }

    if (site.data$CHR %in% c("X", "Y", "chrX", "chrY") && !is.null(gender.data)) {
        sample.genders <- gender.data[rownames(bin.data)]
        gender.colors <- c("lightblue3", "darkred")[sample.genders]
        xvals <- bin.data$b-0.5
        yvals <- bin.data$off-0.6
        points(xvals, yvals, pch=22, col=gender.colors, bg=gender.colors, cex=0.8)
    }

    if (!is.null(pair.data)) {
        # for all > 0, plot the pair count at the correct depth bin and offset
        pair.counts <- sapply(rownames(bin.data),
                              function(sample) { v <- pair.data[cnp,sample]; ifelse(is.null(v),NA,v) })
        xvals <- bin.data$b[pair.counts > 0]-0.5
        yvals <- bin.data$off[pair.counts > 0]-1
        labels <- pair.counts[pair.counts > 0]
        if (length(labels) > 0) {
            text(xvals,yvals,pos=3,labels=labels,font=2,cex=0.8)
        }
    }

    if (!is.null(model) && !pretty) {
        # plot model curves and dashed vertical lines at the means
        variance.factor = site.data$EXPMEAN

        # print(c("weights",model$WEIGHTS))
        # print(c("means",model$MEANS))
        # print(c("variances",model$VARIANCES))
        # print(c("std",sqrt(model$VARIANCES)))
        # print(c("vfactor", variance.factor))

        mapply(function (m,s,n) {
                   yscale <- n/plot.x.scale
                   range <- seq(m-4*s, m+4*s, length.out=100)
                   lines(0.5+plot.x.scale*range,yscale*dnorm(range,m,s),col="black", lwd=2.5)
                   plotVertLine(0.5+plot.x.scale*m,col="black",lty="dashed")
               },
               model$MEANS,
               sqrt(model$VARIANCES/variance.factor),
               model$WEIGHTS*length(samples))
    }
}

makeHist <- function(domain, value, data, idx, func) {
    result <- rep(value, length(domain))
    names(result) <- domain
    tapply.result <- tapply(data, idx, func)
    result[names(tapply.result)] <- tapply.result
    # remove NaNs from the input data, also infinite values that can happen in pathological cases (zero expectation)
    result <- result[!(names(result) %in% c("NaN","Inf"))]
    return(result)
}

asNumericInf <- function(v) {
    v[v == "Infinity"] <- Inf
    v[v == "-Infinity"] <- -Inf
    return(as.numeric(v))
}

plotHorizLine <- function(coord, color="black", lty="solid") {
    width <- 2
    if (names(dev.cur()) == "pdf") {
        width <- 1
    }
    lines(par("usr")[1:2],rep(coord,2),col=color,lwd=width,lty=lty)
}

plotVertLine <- function(coord, color="black", lty="solid") {
    width <- 2
    if (names(dev.cur()) == "pdf") {
        width <- 1
    }
    lines(rep(coord,2),par("usr")[3:4],col=color,lwd=width,lty=lty)
}

parseModelData <- function(site.data) {
    if (is.na(site.data$M1) | is.na(site.data$MWEIGHTS)) {
        return(NULL)
    }

    weights = as.numeric(strsplit(site.data$MWEIGHTS,",")[[1]])

    nClusters = length(weights)
    means = as.numeric(site.data$M1) * seq(0, nClusters-1)
    variances = as.numeric(strsplit(site.data$M2,",")[[1]])
    if (is.na(variances) || length(variances) != 2) {
        cat("The column M2 must contain 2 variance values\n")
	return(NULL)
    }
    if (nClusters > 2) {
        variances = c(variances, variances[2] * seq(2, nClusters-1))
    }

    model = list(NCLUSTERS=nClusters,
                 WEIGHTS=weights,
                 MEANS=means,
                 VARIANCES=variances)
    return(model)
}

loadDataFile <- function() {

    data = read.table(cmdArguments$dataFile, sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)

    if (nrow(data) == 0) {
        return()
    }

    numInfoFields = 13
    site.data <<- data[, 1:numInfoFields]
    sites = rownames(data)
    samples = colnames(data)[(numInfoFields+1):ncol(data)]

    # Parse out samples' coverage:genotype:genotype_qual:pair_count data
    data = data[, samples]
    mtx = apply(data, c(1,2), function(s) {
        v = strsplit(s, ":")[[1]]
	v = sub("NA", NA, v)
	return(as.numeric((v)))
    })
    coverage.data <<- createDataFrame(mtx[1, , ], sites, samples)
    genotype.data <<- createDataFrame(mtx[2, , ], sites, samples)
    genotype.qual.data <<- createDataFrame(mtx[3, , ], sites, samples)
    pair.data <<- createDataFrame(mtx[4, , ], sites, samples)
}

createDataFrame <- function(mtx, sites, samples) {
    if (length(sites) == 1) {
        mtx = t(mtx)
    }
    df = as.data.frame(mtx)
    rownames(df) = sites
    colnames(df) = samples
    return(df)
}

loadTruthData <- function() {
    if (is.null(cmdArguments$truthDataFile)) {
        return(NULL)
    }
    file.data = read.table(cmdArguments$truthDataFile, sep="\t", header=T, row.names=1, stringsAsFactors=F)
    return(file.data)
}

loadSamples <- function() {
    if (is.null(cmdArguments$sampleFile)) {
        return(NULL)
    }
    file.data = read.table(cmdArguments$sampleFile, sep="\t", stringsAsFactors=F)
    return(file.data$V1)
}

loadSampleGroups <- function() {
    if (is.null(cmdArguments$sampleGroupFile)) {
        return(NULL)
    }
    file.data = read.table(cmdArguments$sampleGroupFile, sep="\t", header=T, stringsAsFactors=F)
    rownames(file.data) = file.data$SAMPLE
    return(file.data)
}

loadExcludedSiteSamples <- function() {
    if (is.null(cmdArguments$excludedSiteSampleFile)) {
        return(NULL)
    }
    file.data = read.table(cmdArguments$excludedSiteSampleFile, sep="\t", header=T, stringsAsFactors=F)
    return(file.data)
}

loadGenderMap <- function() {
    if (is.null(cmdArguments$genderMapFile)) {
        return(NULL)
    }
    file.data = read.table(cmdArguments$genderMapFile, sep="\t", stringsAsFactors=F)
    if ("SAMPLE" %in% names(file.data) && "GENDER" %in% names(file.data)) {
        gender.data = file.data$GENDER
        names(gender.data) = file.data$SAMPLE
    } else {
        gender.data = file.data[[2]]
        names(gender.data) = file.data[[1]]
        gender.data = setdiff(gender.data, "GENDER")
    }
    gender.data[toupper(substr(gender.data,1,1)) %in% c("M","Male","1")] = 1
    gender.data[toupper(substr(gender.data,1,1)) %in% c("F","Female","2")] = 2
    gender.map = as.numeric(gender.data)
    names(gender.map) = names(gender.data)
    return(gender.map)
}

getSamples <- function(genotype.data, sample.data, excluded.site.samples) {
    samples = names(genotype.data)
    if (!is.null(sample.data)) {
        included.samples = make.names(sample.data)
        samples = intersect(samples, included.samples)
    }
    if (!is.null(excluded.site.samples)) {
        excluded.site.samples = make.names(excluded.site.samples)
        samples = setdiff(samples, excluded.site.samples)
    }
    return(samples)
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
            result[[as.character(keyword)]] = args[[argpos]]
            argpos = argpos + 1
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
