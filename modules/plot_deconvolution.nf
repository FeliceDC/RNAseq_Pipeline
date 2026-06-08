process PLOT_DECONVOLUTION {
    tag "Plot TME"
    label 'process_low'

    container 'python:3.10-slim'

    input:
    path immucellai_results

    output:
    path "*.pdf", emit: plots
    path "*_mqc.png", emit: multiqc_png, optional: true

    script:
    """
   
    pip install --no-cache-dir pandas openpyxl matplotlib seaborn scipy

   
    python ${projectDir}/bin/plot_deconvolution.py ${immucellai_results}
    """
}
