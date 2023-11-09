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
	
**Snakefile** (stored in the [docker](docker/) folder)
```shell
# Define the paths and filenames
reference_genome = "data/GRCh38.fa"
input_sequences = "data/library.fa"
output_sam = "sgRNAlib.sam"
output_with_gene_names = "sequence_mapping_info_with_gene_names.txt"
gene_annotation_file = "data/gencode.v44.annotation.gtf"
unique_genes_output = "unique_gene_names.txt"

# Rule for performing the alignment
rule bowtie2_align:
    input:
        reference_genome=reference_genome,
        input_sequences=input_sequences
    output:
        output_sam=output_sam
    shell:
        """
        bowtie2-build "{input.reference_genome}" GRCh38_index
        bowtie2 -x GRCh38_index -f "{input.input_sequences}" -S "{output.output_sam}"
        """


```

**Important note: please keep in mind that all the code was created, run and tested on a MacOS Ventura (13.6.1) device**

## Create Docker image for the pipeline

- If running on a MacOS, first make sure to install [Docker Desktop](https://docs.docker.com/desktop/install/mac-install/)
- Create a Dockerfile (stored here in the [docker](docker/) folder) with the following content (where data is your local input data directory - please modify accordingly):
```shell
FROM continuumio/miniconda3

RUN conda install -c bioconda bowtie2 samtools bedtools

RUN conda install -c conda-forge mamba

RUN mamba install -c bioconda snakemake

RUN mkdir /workflow
COPY Snakefile /workflow/Snakefile

# Create a directory in the container to store data
RUN mkdir /workflow/data

# Copy the necessary input data into the container
COPY data/GRCh38.fa /workflow/data/GRCh38.fa
COPY data/library.fa /workflow/data/library.fa
COPY data/gencode.v44.annotation.gtf /workflow/data/gencode.v44.annotation.gtf

WORKDIR /workflow

```
- Build the Docker Image
```shell
docker build -t sgRNAlib_snakemake_image .

```
- Run the Docker Container
```shell
docker run --rm -it -v /path/to/local/data/directory:/workflow/data sgRNAlib_snakemake_image

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

