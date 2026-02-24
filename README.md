# Comp_Bio_project

This project provides a Snakemake pipeline for analyzing the transcriptome of Human Cytomegalovirus (HCMV) across two conditions (2 days post-infection and 6 days post-infection).

# Dependencies
To run this pipeline, you must have the following tools installed: Python 3.x, Snakemake, Biopython, Kallisto, Sleuth (R package), Bowtie2, SPAdes, NCBI Datasets CLI, BLAST+.

# Installation
If using Conda, you can install most dependencies with:

```
conda install -c bioconda snakemake kallisto bowtie2 spades blast-plus biopython
```

To install R, use the following link: [R website](https://www.r-project.org/)
# How to run the code
1. Clone this repository:
```
git clone https://github.com/sofiaorellanar-boop/Comp_Bio_project.git
cd Comp_Bio_project
```
2. Place your raw fastq files in a directory named data/ (test data set is already included).
3. Ensure your metadata.txt file is present in the main directory. It should look something like this with tab delimited spaces:

sample	condition   path
SRR5660030	2dpi    results/SRR5660030
SRR5660033	6dpi    results/SRR5660033
SRR5660044	2dpi    results/SRR5660044
SRR5660045	6dpi    results/SRR5660045

4. Run the pipeline using Snakemake:
```
snakemake --cores 4
```