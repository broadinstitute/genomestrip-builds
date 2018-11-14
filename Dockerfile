# Dockerfile for Genome STRiP
FROM ubuntu:16.04
MAINTAINER Seva Kashin
MAINTAINER Bob Handsaker

RUN apt-get -qq update && apt-get install -qqy \
    build-essential \
    lbzip2 \
    openjdk-8-jdk \
    r-base-core \
    unzip \
    vim-common \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# google cloud sdk
RUN apt-get -qq update && apt-get install -qqy curl
RUN curl -sSL https://sdk.cloud.google.com | bash

ENV PATH=${PATH}:/root/google-cloud-sdk/bin

# R
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('argparser')" && \
    Rscript -e "install.packages('reshape2', dependencies=TRUE)" && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("qvalue");'

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar -xf samtools-1.3.1.tar.bz2 && rm samtools-1.3.1.tar.bz2 && cd samtools-1.3.1 && make && make install && make clean

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar -xf htslib-1.3.2.tar.bz2 && rm htslib-1.3.2.tar.bz2 && cd htslib-1.3.2 && make && make install && make clean

# Genome STRiP
#ADD svtoolkit.tar.gz /opt
COPY svtoolkit /opt/svtoolkit
ENV SV_DIR=/opt/svtoolkit
ENV SV_CLASSPATH="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar"
ENV SV_CLASSPATH="${SV_CLASSPATH}:${SV_DIR}/lib/depend/aws-java-sdk-1.3.23.jar"
ENV SV_CLASSPATH="${SV_CLASSPATH}:${SV_DIR}/lib/depend/shaded-google-nio-all.jar"
ENV LD_LIBRARY_PATH="${SV_DIR}/bwa:${LD_LIBRARY_PATH}"
ENV GATK_JAR="${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar"
ENV PATH="${PATH}:${SV_DIR}/scripts/firecloud/util"

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
