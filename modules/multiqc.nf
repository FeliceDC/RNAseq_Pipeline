process MULTIQC {
    label 'process_high'

    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'

    input:
    path multiqc_files 

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data

    script:
    """
    multiqc .
    """
}
