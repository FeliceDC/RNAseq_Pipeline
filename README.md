<h1 align="center">
   FeliceDC/RNAseq_Pipeline
</h1>

<p align="center">
  <img src="https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg" alt="Nextflow">
  <img src="https://img.shields.io/badge/status-active-success.svg" alt="Status">
  <img src="https://img.shields.io/badge/Bioinformatics-RNA--seq-blue" alt="Bioinformatics">
</p>

## Introduction
**FeliceDC/RNAseq_Pipeline** is a bioinformatics analysis pipeline used for RNA sequencing data. Developed in Nextflow, it automates the entire workflow from raw FASTQ reads to Differential Expression analysis, ensuring reproducibility and scalability.

The pipeline is built using Docker containers, meaning you don't need to install any bioinformatics tools manually.

This pipeline can automatically handles datasets with very small genomes (e.g., viral genomes or targeted panels) by bypassing standard dispersion curves in DESeq2, preventing crashes common in other public workflows.

## Pipeline Summary
1. Raw read QC (`FastQC`)
2. Adapter and quality trimming (`Trim Galore!`)
3. Read alignment and indexing (`STAR`)
4. Gene-level quantification (`featureCounts`)
5. Pipeline QC report (`MultiQC`)
6. Differential Expression Analysis & Visualization (`DESeq2`)


## Usage

To run the pipeline on your own samples, you need to provide:
1. Your raw fastq.gz files
2. A reference genome
3. An annotation file
4. A design matrix (named "samplesheet"). The samplesheet must be a comma-separated values file (.csv). The first column (called "sample") must match the FASTQ file names (excluding the _1.fastq.gz suffix).The second column should contain the variable you want to use for differential analysis. Optionally, you can perform differential analysis using two variables if needed.
Example:

**samplesheet.csv**
```bash
sample,condition,age,library_selection
SRR8518319,normal_adiacent,52,cDNA
SRR8518327,tumor,37,cDNA
SRR8518335,normal_adiacent,62,cDNA
SRR8518360,tumor,54,cDNA
```

Once the samplesheet has been created, make sure you have ready the samplesheet, FASTQ files, genome, and GTF files paths.

Now you should be ready to run the pipeline.
>[!NOTE]
>An example running code is
>```bash
>nextflow run Filic03/RNAseq_Pipeline --input_reads "/apps/Felice/GSE/prova_per_nextflow/*_{1,2}.fastq.gz" --fasta "/apps/Felice/GSE/nuovo_genoma/GRCh38.primary_assembly.genome.fa" --gtf "/apps/Felice/GSE/nuovo_genoma/gencode.v49.primary_assembly.annotation.gtf" --design "condition" --samplesheet "./samplesheet.csv" --max_cpus 20
>```
>
>If you want, you can run Deseq2 with two variables. Then you have to write --design "variable1 + variable2"
>```bash
>nextflow run Filic03/RNAseq_Pipeline --input_reads "/apps/Felice/GSE/prova_per_nextflow/*_{1,2}.fastq.gz" --fasta "/apps/Felice/GSE/nuovo_genoma/GRCh38.primary_assembly.genome.fa" --gtf "/apps/Felice/GSE/nuovo_genoma/gencode.v49.primary_assembly.annotation.gtf" --design "condition + age" --samplesheet "./samplesheet.csv" --max_cpus 12
>```

>[!WARNING]
>Running the pipeline on full human datasets requires significant computational resources. It is highly recommended to check your machine and specify an appropriate --max_cpus limit.

## Output Structure
By default, the pipeline creates a results/ directory containing the following sub-directories:

- fastqc/ and multiqc/: Interactive HTML quality reports.

- star/: Sorted .bam files ready for IGV visualization.

- featurecounts/: Raw count matrices.

- deseq2/: CSV tables with statistically significant Differentially Expressed Genes (DEGs) and related plots (MA plot, PCA, Volcano plot, Heatmap ecc.).



