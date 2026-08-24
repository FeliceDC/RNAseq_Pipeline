process DESEQ2 {
    tag "Differential Analysis"
    label 'process_high'
   
    container 'quay.io/biocontainers/mulled-v2-8849acf39a43ddc833ad691d6ff400e193fa167b:a92e4a61c5638d3d93e624ddec11559648908476-0'

    input:
    path counts     
    path samplesheet  

    output:
    path "*.{csv,txt}", emit: results_tables
    path "*.pdf", emit: results_pdf
    path "*_mqc.png", emit: multiqc_png, optional: true

    script:
    """
    Rscript ${projectDir}/bin/run_deseq2.R $counts $samplesheet "${params.design}" ${params.deseq2_pvalue} ${params.deseq2_logfc}
    """
}
