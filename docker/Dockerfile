# Base image
FROM ubuntu:latest

# Install necessary libraries and tools
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    git \
    mercurial \
    subversion

# Download and install Miniconda
RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda clean -afy

# Add Miniconda to the PATH environment variable
ENV PATH /opt/conda/bin:$PATH

# Set up channels for package installation
RUN conda config --set ssl_verify no && \
    conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

# Install required packages
RUN conda install --yes --insecure bowtie2 samtools bedtools

# Install additional tools
RUN conda install --yes --insecure -c conda-forge mamba

# Set up the directory structure and copy necessary files
RUN mkdir /workflow
COPY Snakefile /workflow/Snakefile

# Create a directory in the container to store data
RUN mkdir /workflow/data

# Copy the necessary input data into the container
COPY data/GRCh38.fa /workflow/data/GRCh38.fa
COPY data/library.fa /workflow/data/library.fa
COPY data/gencode.v44.annotation.gtf /workflow/data/gencode.v44.annotation.gtf

# Set the working directory
WORKDIR /workflow
