process IMSIG {
    tag "ImSig Deconvolution"
    label 'process_medium'

    container 'rocker/geospatial:4.3.1'

    input:
    path counts_file

    output:
    path "ImSig_results.csv", emit: results
    path "ImSig_plot.pdf"   , emit: plot
    path "*_mqc.png", emit: multiqc_png, optional: true

    script:
    """
    Rscript ${projectDir}/bin/run_imsig.R ${counts_file} || 
    touch ImSig_results.csv 
    touch ImSig_plot.pdf
    """
}
