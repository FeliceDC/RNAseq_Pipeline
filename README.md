<h1 align="center">
   FeliceDC/RNAseq_Pipeline
</h1>

<p align="center">
  <img src="https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg" alt="Nextflow">
  <img src="https://img.shields.io/badge/status-active-success.svg" alt="Status">
  <img src="https://img.shields.io/badge/Bioinformatics-RNA--seq-blue" alt="Bioinformatics">
</p>

## Introduction
**FeliceDC/RNAseq_Pipeline** is a comprehensive, modular bioinformatics analysis pipeline used for RNA sequencing data. Developed in Nextflow (DSL2), it automates the entire workflow from raw FASTQ reads to advanced downstream analysis (Differential Expression, Splicing, Fusions, and Deconvolution), ensuring reproducibility and scalability.

The pipeline is built using Docker containers, meaning you don't need to install any bioinformatics tools manually.

## Pipeline Summary
1. Raw read QC (`FastQC`)
2. Adapter and quality trimming (`Trim Galore!`)
3. Read alignment and indexing (`STAR`)
4. Gene-level quantification (`featureCounts`)
5. Pipeline QC report (`MultiQC`)
6. Differential Expression Analysis (`DESeq2`), followed by pathway enrichment (`EnrichR`)
7. Tumor Deconvolution: Immune and stromal cell infiltration estimation ((`ImmuCellAI`) and (`ImSig`)).
8. Alternative Splicing: Classical statistical splicing analysis (`rMATS`) compared with Bayesian Deep Learning predictions (`DARTS`), complete with automated Volcano, Bar, and Sashimi plots.
9. Gene Fusions: Structural variant detection (`Arriba`)


## Usage

To run the pipeline on your own samples, you need to provide:
1. Your raw fastq.gz files
2. A reference genome
3. An annotation file
4. A design matrix (named "samplesheet").

The samplesheet must be a comma-separated values file (.csv). The first column (called "sample") must match the FASTQ file names (excluding the _1.fastq.gz suffix).The second column should contain the variable you want to use for differential analysis. Optionally, you can perform differential analysis using two variables if needed.
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
>nextflow run Filic03/RNAseq_Pipeline --input_reads "/Your/Files/Path/*fastq.gz" --fasta "/Your/Genome/Path/GRCh38.primary_assembly.genome.fa" --gtf "/Your/Path/gencode.v49.primary_assembly.annotation.gtf" --design "condition" --samplesheet "/Your/File/Path/samplesheet.csv"
>```
>
>If you want, you can run Deseq2 with two variables. Then you have to write --design "variable1 + variable2"
>```bash
>nextflow run Filic03/RNAseq_Pipeline --input_reads "/Your/Files/Path/*fastq.gz" --fasta "/Your/Genome/Path/GRCh38.primary_assembly.genome.fa" --gtf "/Your/Path/gencode.v49.primary_assembly.annotation.gtf" --design "condition + age" --samplesheet "/Your/File/Path/samplesheet.csv"
>```

 | Parametro | Descrizione |
| :--- | :--- |
| `--skip_fusions` | Salta l'analisi delle fusioni geniche (Arriba) |
| `--skip_deconvolution` | Salta la stima dell'infiltrazione cellulare |
| `--skip_differential` | Salta l'analisi differenziale (DESeq2/EnrichR) |

</select>

>[!WARNING]
>Running the pipeline on full human datasets requires significant computational resources. It is highly recommended to check your machine and specify an appropriate --max_cpus limit.

## Output Structure
By default, the pipeline creates a results/ directory containing the following sub-directories:

- fastqc/ and multiqc/: Interactive HTML quality reports.

- star/: Sorted .bam files ready for IGV visualization.

- featurecounts/: Raw count matrices.

- deseq2/: CSV tables with statistically significant Differentially Expressed Genes (DEGs) and related plots (MA plot, PCA, Volcano plot, Heatmap ecc.).



