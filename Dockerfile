# base image:
FROM ubuntu:22.04 as learn_doker
ARG DEBIAN_FRONTEND=noninteractive

# install packages:
RUN apt-get update && apt-get install -y \
    gcc-9 \
    g++-9 \
    make \
    git \
    wget \
    zlib1g-dev \
    build-essential autoconf automake libtool \
    libncurses5-dev libbz2-dev liblzma-dev \
    openjdk-11-jdk \
    unzip \
    curl \
    ca-certificates	

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 100

# install dependencies:

# install python 3.8:
WORKDIR /dependencies
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3.8 python3-pip python3.8-distutils && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 100 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.8 100

# install fastp (https://github.com/OpenGene/fastp?tab=readme-ov-file#or-download-the-latest-prebuilt-binary-for-linux-users):
WORKDIR /dependencies
RUN wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp && \
    mv fastp /usr/local/bin/fastp

# install fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/):
WORKDIR /dependencies
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod a+x ./FastQC/fastqc && \
    mv FastQC /usr/local/bin/FastQC && \
    ln -s /usr/local/bin/FastQC/fastqc /usr/local/bin/fastqc

# install multiqc (https://multiqc.info/docs/getting_started/installation/#pip--pypi):
WORKDIR /dependencies
RUN pip install multiqc

## intall bwa (https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2):
WORKDIR /dependencies
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    mv bwa /usr/local/bin/bwa

## install samtools (https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2):
WORKDIR /dependencies
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xvf samtools-1.9.tar.bz2 && \
    rm samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

## install picard (https://github.com/broadinstitute/picard/releases/download/2.25.5/picard.jar):
WORKDIR /dependencies
RUN wget https://github.com/broadinstitute/picard/releases/download/2.25.5/picard.jar && \
    mv picard.jar /usr/local/bin/picard.jar && \
    echo "java -jar /usr/local/bin/picard.jar" > /usr/local/bin/picard && \
    chmod a+x /usr/local/bin/picard

# install gatk 3.7 (https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2):
WORKDIR /dependencies
RUN wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2 && \
    tar -xvf GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2 && \
    rm GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2 && \
    mv GenomeAnalysisTK.jar /usr/local/bin/GenomeAnalysisTK.jar && \
    echo "java -jar /usr/local/bin/GenomeAnalysisTK.jar" > /usr/local/bin/gatk3.7 && \
    chmod a+x /usr/local/bin/gatk3.7

# install qualimap (https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip):
WORKDIR /dependencies
RUN wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip && \
    unzip qualimap_v2.3.zip && \
    mv qualimap_v2.3 /usr/local/bin/qualimap_v2.3 && \
    ln -s /usr/local/bin/qualimap_v2.3/qualimap /usr/local/bin/qualimap

# install gatk 4.5.0.0 (https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip):
WORKDIR /dependencies
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip && \
    unzip gatk-4.5.0.0.zip && \
    mv gatk-4.5.0.0 /usr/local/bin/gatk-4.5.0.0 && \
    ln -s /usr/local/bin/gatk-4.5.0.0/gatk /usr/local/bin/gatk

# install bcftools (https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2)
WORKDIR /dependencies
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xvf bcftools-1.19.tar.bz2 && \
    rm bcftools-1.19.tar.bz2 && \
    cd bcftools-1.19 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

# install bedtools (https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz)
WORKDIR /dependencies
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -xvf bedtools-2.31.1.tar.gz && \
    rm bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin/

# install stacks (https://catchenlab.life.illinois.edu/stacks/source/stacks-2.6b.tar.gz)
WORKDIR /dependencies
RUN wget https://catchenlab.life.illinois.edu/stacks/source/stacks-2.66.tar.gz && \
    tar -xvf stacks-2.66.tar.gz && \
    rm stacks-2.66.tar.gz && \
    cd stacks-2.66 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

# install plink2 (https://www.cog-genomics.org/plink/2.0/)
WORKDIR /dependencies
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20240318.zip && \
    unzip plink2_linux_x86_64_20240318.zip && \
    mv plink2 /usr/local/bin/plink2

# install plink1.9 (https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip)
WORKDIR /dependencies
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && \
    unzip plink_linux_x86_64_20231211.zip && \
    mv plink /usr/local/bin/plink

# install R:
WORKDIR /dependencies
RUN apt-get update && \
    apt-get install -y r-base r-base-core r-recommended r-base-dev && \
    apt-get clean

# install R packages:
RUN apt install -y r-cran-adegenet
RUN apt install -y r-cran-factominer
# RUN ln -s /usr/lib/aarch64-linux-gnu/libgfortran.so.5.0.0 /usr/lib/libgfortran.so
# # RUN Rscript -e 'install.packages("pegas", repos="https://cran.rstudio.com")'
# RUN apt install -y r-cran-devtools
# RUN Rscript -e 'devtools::install_github("pievos101/PopGenome")'
# 
# # RUN apt-get install libudunits2-dev
# RUN apt-get install libgdal-dev

# https://www-users.york.ac.uk/~dj757/popgenomics/workshop5.html
# install vcftools (https://github.com/vcftools/vcftools/tarball/master):
WORKDIR /dependencies
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xvf vcftools-0.1.16.tar.gz && \
    rm vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

# https://speciationgenomics.github.io/ADMIXTURE/
# https://www.popgen.dk/software/index.php/EvalAdmix
# install admixture (https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz)
WORKDIR /dependencies
RUN wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz && \
    tar -xvf admixture_linux-1.3.0.tar.gz && \
    rm admixture_linux-1.3.0.tar.gz && \
    mv dist/admixture_linux-1.3.0/admixture /usr/local/bin/admixture

# install evaladmix (https://github.com/GenisGE/evalAdmix.git)
WORKDIR /dependencies
RUN git clone https://github.com/GenisGE/evalAdmix.git && \
    cd evalAdmix && \
    make && \
    mv evalAdmix /usr/local/bin/evalAdmix

# install gone Linux (https://github.com/esrud/GONE/tree/master/Linux)
WORKDIR /dependencies
RUN git clone https://github.com/esrud/GONE.git && \
    cp -r GONE/Linux /root/GONE-Linux

# install tabix (https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2)
WORKDIR /dependencies
RUN wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 && \
    tar -xvf htslib-1.19.1.tar.bz2 && \
    rm htslib-1.19.1.tar.bz2 && \
    cd htslib-1.19.1 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

# RUN apt-get update
# RUN apt-get install gcc-11.3 g++-11.3
# RUN Rscript -e 'install.packages("pegas", repos="https://cran.rstudio.com")'
# RUN wget https://cran.r-project.org/src/contrib/Archive/pegas/pegas_1.2.tar.gz && \
#     R CMD INSTALL pegas_1.2.tar.gz && \
#     rm pegas_1.3.tar.gz

# # install igv (https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_2.16.0.zip)
# RUN apt-get update && apt-get install -y \
#     xauth \ 
#     libxrender1 \
#     libxtst6 \
#     libxi6
# 
# WORKDIR /dependencies
# RUN wget https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_2.16.0.zip && \
#     unzip IGV_2.16.0.zip && \
#     mv IGV_2.16.0 /usr/local/bin/IGV_2.16.0 && \
#     ln -s /usr/local/bin/IGV_2.16.0/igv.sh /usr/local/bin/igv

# # clone the course repository:
# WORKDIR /dependencies
# RUN git clone https://github.com/Enricobazzi/EBD-PopGenIntro.git && \
#     mv EBD-PopGenIntro /root/EBD-PopGenIntro

WORKDIR /root
