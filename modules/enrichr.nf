process ENRICHR {
    tag "Pathway Analysis su ${deseq2_results.baseName}"
    label 'process_low'
    
    container 'rocker/geospatial:4.3.1'
    
    input:
    path deseq2_results

    output:
    path "*.{csv,pdf}", emit: enrichr_results, optional: true
    path "*_mqc.png", emit: multiqc_png, optional: true

    script:
    """
    Rscript ${projectDir}/bin/run_enrichr.R \\
        --input ${deseq2_results} \\
        --databases "${params.enrichr_database}" \\
        --outdir .
    """
}
