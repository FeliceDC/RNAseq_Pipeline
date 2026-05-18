process FEATURECOUNTS {
    tag "All Samples"
    label 'process_high'
    container 'quay.io/biocontainers/subread:2.0.6--he4a0461_0'

    input:
    path gtf
    path bams

    output:
    path "*.txt"        , emit: counts
    path "*.txt.summary", emit: summary

    script:
    
    def paired = params.single_end ? "" : "-p"
    """
    featureCounts $paired -a $gtf -o counts.txt -T ${task.cpus} -g ${params.g} -s ${params.strandedness} $bams
    """
}
