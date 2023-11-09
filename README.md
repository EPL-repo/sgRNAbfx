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
	1. run the mapping (Task 1) with Bowtie2 after creating reference index,
	2. get chromosome name, start-end positions and strand for each mapped sgRNA sequence (Task 2),
	3. assign gene names based on mapping results and compare the obtained annotation with gene names provided within description lines of fasta file (Task 3)
	4. create a final output of unique, mapped gene names (official symbols) to be used in Task 4
	
**Snakefile**
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

## Create Docker image for the pipeline

- First step is to create a Dockerfile with the following content:


**Important note: please keep in mind that all the code was created, run and tested on a MacOS Ventura (13.6.1) device**

## R Jupyter Notebook to perform Task 4

- It would have been possible to include another process within the Nextflow pipeline to run an R script allowing to retrieve a matrix of gene expression for the gene IDs from Task 3 for samples "TCGA-A7-A13D-01A-13R-A12P-07" and "TCGA-E9-A1RH-11A-34R-A169-07" of the TCGA-BRCA cohort (Task 4). However, to showcase another tool for collaborative work, an R Jupyter Notebook was created to show the steps required to accomplish the task.
- The notebook as well as the outputs can be found in the folder [Task4_Jupyter](Task4_Jupyter/).
- Expression data accessed from the GDC are stored in the [data](data/) folder (along with the manifest file).

