# sgRNALib
sgRNA library characterization small project

---

## Contents

> Directory details

- [Snakemake pipeline for Tasks 1-3](#Snake)
- [Create Docker image for the Nextflow pipeline](#Docker)
- [R Jupyter Notebook to perform Task 4](#RJupyter)

---

## Snakemake pipeline for Tasks 1-3

- To map the sgRNA sequences provided in the multifasta file library.fa (stored in the [data](data/) folder), the human reference GRCh38 was accessed from [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
- In the context of mapping sgRNA sequences to the human genome, the necessity of using a splice-aware aligner depends on the specific requirements of the analysis. sgRNA sequences are typically short sequences designed to target specific genomic loci for CRISPR-Cas9 genome editing. These sequences are not inherently related to splicing events, which involve the removal of introns from pre-mRNA during the process of gene expression. As in the scope of this short we will not be working on assessing potential off-target effects or evaluating the specificity of sgRNA sequences (not enough information provided), it is not necessarily beneficial to use a splice-aware aligner. We will use **Bowtie2** which is very efficient on short DNA sequences, including sgRNAs [Haeussler  *et al*. *Genome Biol* 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2)
- GENCODE v44 was used for comprehensive gene annotation, the gtf file can be found [here](https://www.gencodegenes.org/human/)
- To wrap all the steps required from data processing, a simple snakemake pipeline was set up to:
	1. run the mapping (**Task 1**) with Bowtie2 after creating the reference index,
		- *technical info*: mapping summary
		```shell
  		77736 reads; of these:
		77736 (100.00%) were unpaired; of these:
 		499 (0.64%) aligned 0 times
 		53202 (68.44%) aligned exactly 1 time
 		24035 (30.92%) aligned >1 times
  		99.36% overall alignment rate
  		```
	3. get chromosome name, start-end positions and strand for each mapped sgRNA sequence (**Task 2**),
	4. assign gene names based on mapping results and compare the obtained annotation with gene names provided within description lines of fasta file (**Task 3**)
	5. create a final output of unique, mapped gene names (official symbols) to be used in **Task 4**
- The **Snakefile** is stored in the [docker](docker/) folder
- Main steps wrapped in the snakemake pipeline consist of:
```shell
#!/bin/bash

## Set the paths and filenames
reference_genome="data/GRCh38/GRCh38.fa"
input_sequences="data/library.fa"
output_sam="sgRNAlib.sam"
output_with_gene_names="sequence_mapping_info_with_gene_names.txt"
gene_annotation_file="data/gencode.v44.annotation.gtf"
unique_genes_output="unique_gene_names.txt"

## Perform the alignment with Bowtie2
bowtie2-build "$reference_genome" GRCh38_index
bowtie2 -x GRCh38_index -f "$input_sequences" -S "$output_sam"

## Process the SAM file to get a sorted bam file 
# Convert SAM to BAM
samtools view -bS "$output_sam" > output.bam
# Sort the BAM file
samtools sort output.bam -o sorted_output.bam
# Index the sorted BAM file
samtools index sorted_output.bam

## Extract gene names from reference GTF (needed for gene annotation by mapping)
gawk -F'\t' 'BEGIN {OFS="\t"} $3 == "gene" { if (match($9, /gene_name "([^"]+)"/, m)) { print $1, $4, $5, m[1] } }' "$gene_annotation_file" > genes_with_names.bed

## Use BEDTools to assign gene names to the mapped sequences
bedtools bamtobed -i sorted_output.bam > temp.bed

## Generate a first output with Coordinates and gene names from both the library.fa and gene mapping from the GTF
bedtools intersect -a temp.bed -b genes_with_names.bed -wa -wb | gawk -F '\\t' 'BEGIN {OFS="\t"; print "Chromosome", "Start_Position", "End_Position", "Gene_Name", "Mapping_Info", "Gene_Info"} {split($4, a, "|"); split(a[3], b, "_"); print $1, $2, $3, b[1], $5, $10}' > "$output_with_gene_names"
rm genes_with_names.bed temp.bed

## Add a column with TRUE/FALSE when the genes from the 2 different sources match=TRUE (or not=FALSE)
awk -F'\t' 'BEGIN {OFS="\t"} {print $0, ($4 == $6) ? "TRUE" : "FALSE"}' "$output_with_gene_names" > temp_with_true_false.txt
mv temp_with_true_false.txt "$output_with_gene_names"

## Create the desired final list with unique gene names provided within description lines of library.fa that
## were also found by name assignation after mapping
awk -F'\t' '$7 == "TRUE" {print $4}' ${output_with_gene_names} | sort | uniq > ${unique_genes_output}

```

**Important note: please keep in mind that all the code was created, run and tested on a MacOS Ventura (13.6.1) device**

## Create Docker image for the pipeline

- If running on a MacOS, first make sure to install [Docker Desktop](https://docs.docker.com/desktop/install/mac-install/)
- Create a Dockerfile (stored here in the [docker](docker/) folder) with the following content (where data is your local input data directory - please modify accordingly):
```shell
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

```
- Build the Docker Image
```shell
docker build -t sgrnalib_snakemake_image .

```
- Run the Docker Container
```shell
docker run --rm -it -v /path/to/local/data/directory:/workflow/data sgrnalib_snakemake_image

```
- Once inside the Docker container, execute the Snakemake pipeline
```shell
snakemake --cores all
```

**Important note: please keep in mind that all the code was created, run and tested on a MacOS Ventura (13.6.1) device**

## R Jupyter Notebook to perform Task 4

- It would have been possible to include another process within the Nextflow pipeline to run an R script allowing to retrieve a matrix of gene expression for the gene IDs from Task 3 for samples "TCGA-A7-A13D-01A-13R-A12P-07" and "TCGA-E9-A1RH-11A-34R-A169-07" of the TCGA-BRCA cohort (Task 4). However, to showcase another tool for collaborative work, an R Jupyter Notebook was created to show the steps required to accomplish the task.
- The notebook as well as the outputs can be found in the folder [Task4_Jupyter](Task4_Jupyter/). Please see the **TCGA-BRCA_sgRNA_gene_expression_matrix.html** output file for an interactive html table with the final expression matrix (columns can be sorted, specific genes can be looked for in the search bar) - useful result format to share with non-computational colleagues.
- Expression data accessed from the GDC are stored in the [data](data/) folder (along with the manifest file).

