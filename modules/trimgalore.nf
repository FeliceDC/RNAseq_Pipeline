process TRIMGALORE {
    tag "Trimming on $sample_id"
    label 'process_high'

    container 'quay.io/biocontainers/trim-galore:0.6.11--hdfd78af_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    
    tuple val(sample_id), path("*{trimmed,val}*.fq.gz"), emit: reads
    path "*_trimming_report.txt"                        , emit: log

    script:
    
    def single_end_flag = params.single_end ? "" : "--paired"
    
    """
    trim_galore $single_end_flag --cores ${task.cpus} $reads
    """
}
