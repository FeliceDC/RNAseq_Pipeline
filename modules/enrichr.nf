process ENRICHR {
    tag "Pathway Analysis"
    label 'process_low'
    
    container 'rocker/geospatial:4.3.1'
    input:
    path deseq2_results

    output:
    path "*.{csv,pdf}", emit: enrichr_results

    script:
    """
    Rscript ${projectDir}/bin/run_enrichr.R \\
        --input filtered_results.txt \\
        --databases "${params.enrichr_database}" \\
        --outdir .
    """
}
