process ARRIBA {
    tag "Arriba on ${sample_id}"
    label 'process_high'


    container 'uhrigs/arriba:2.4.0'

    input:

    tuple val(sample_id), path(bam)
    path fasta
    path gtf

    output:

    path "*_fusions.tsv", emit: fusions
    path "*_fusions.discarded.tsv", emit: discarded

    script:
    """
    /arriba_v2.4.0/arriba \\
        -x ${bam} \\
        -a ${fasta} \\
        -g ${gtf} \\
        -o ${sample_id}_fusions.tsv \\
        -O ${sample_id}_fusions.discarded.tsv
    """
}
