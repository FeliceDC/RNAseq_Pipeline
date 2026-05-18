process IMMUCELLAI {
    tag "TME Deconvolution v2"
    label 'process_high'

    container 'python:3.10-slim'

    input:
    path featurecounts_output

    output:
    path "tpm_matrix.txt"           , emit: tpm_matrix
    path "ImmuCellAI2_results.xlsx"  , emit: fractions

    script:
    """
  pip install --no-cache-dir --default-timeout=1000 pandas numba scipy tqdm joblib scikit-learn dask distributed immucellai2

python ${projectDir}/bin/run_immucellai2.py ${featurecounts_output} ${task.cpus}
    """
}
