process RMATS {
    tag "Alternative Splicing"
    label 'process_high'

    container 'quay.io/biocontainers/rmats:4.1.2--py39h6d91be2_2'

    input:
    path bams_cond1
    path bams_cond2
    path gtf

    output:
    path "rmats_out/*", emit: splicing_results

    script:
    """
    echo \${bams_cond1.join(',')} > b1.txt
    echo \${bams_cond2.join(',')} > b2.txt

    rmats.py \\
        --b1 b1.txt \\
        --b2 b2.txt \\
        --gtf \${gtf} \\
        -t paired \\
        --readLength 100 \\
        --nthread \${task.cpus} \\
        --od rmats_out \\
        --tmp rmats_tmp
    """
}
