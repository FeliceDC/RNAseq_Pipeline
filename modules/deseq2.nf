process DESEQ2 {
    tag "Differential Analysis"
    label 'process_high'
   
    container 'quay.io/biocontainers/bioconductor-deseq2:1.50.2--r45ha27e39d_0'

    input:
    path counts     
    path samplesheet  

    output:
    path "*.{csv,txt}", emit: results_tables
    path "*.pdf", emit: results_pdf

    script:
    """
    Rscript ${projectDir}/bin/run_deseq2.R $counts $samplesheet "${params.design}"
    """
}
