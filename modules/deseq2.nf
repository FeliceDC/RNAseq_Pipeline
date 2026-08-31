process DESEQ2 {
    tag "Differential Analysis"
    label 'process_high'
   
    container 'quay.io/biocontainers/bioconductor-deseq2:1.50.2--r45ha27e39d_0'

    input:
    path counts     
    path samplesheet  

    output:
output:
    path "filtered_results_*.txt", emit: filtered_results
    path "deseq2_results_*.txt", emit: raw_results
    path "complete_table_*.txt", emit: complete_tables
    path "deseq2_plots.pdf", emit: results_pdf
    path "deseq2_volcano_*_mqc.png", emit: multiqc_png
    path "deseq2_pca_*_mqc.png", emit: pca_png

    script:
    """
    Rscript ${projectDir}/bin/run_deseq2.R $counts $samplesheet "${params.design}" ${params.deseq2_pvalue} ${params.deseq2_logfc}
    """
}
