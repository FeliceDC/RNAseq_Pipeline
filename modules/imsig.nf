process IMSIG {
    tag "ImSig Deconvolution"
    label 'process_medium'

    container 'rocker/r-ver:4.3.1'

    input:
    path counts_file

    output:
    path "ImSig_results.csv", emit: results

    script:
    """
    Rscript ${projectDir}/bin/run_imsig.R ${counts_file}
    """
}
