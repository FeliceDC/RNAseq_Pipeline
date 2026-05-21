process IMSIG {
    tag "ImSig Deconvolution"
    label 'process_medium'

    container 'rocker/geospatial:4.3.1'

    input:
    path counts_file

    output:
    path "ImSig_results.csv", emit: results
    path "ImSig_plot.pdf"   , emit: plot

    script:
    """
    Rscript ${projectDir}/bin/run_imsig.R ${counts_file}
    """
}
