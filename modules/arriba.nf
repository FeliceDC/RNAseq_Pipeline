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
    path "*.pdf", emit: plots
    script:
    """
    /arriba_v2.4.0/arriba \\
        -x ${bam} \\
        -a ${fasta} \\
        -g ${gtf} \\
        -b /arriba_v2.4.0/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz \\
        -k /arriba_v2.4.0/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz \\
        -p /arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 \\
        -o ${sample_id}_fusions.tsv \\
        -O ${sample_id}_fusions.discarded.tsv > ${sample_id}.arriba.log 2>&1



    /arriba_v2.4.0/draw_fusions.R \\
        --fusions=${sample_id}_fusions.tsv \\
        --alignments=${bam} \\
        --annotation=${gtf} \\
        --cytobands=/arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv \\
        --proteinDomains=/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 \\
        --output=${sample_id}_fusions.pdf || touch ${sample_id}_fusions.pdf
    """

}
