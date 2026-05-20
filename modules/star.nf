process STAR_INDEX {
    tag "$fasta"
    label 'process_high'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    """
    mkdir star_index
    STAR --runMode genomeGenerate \\
         --genomeDir star_index \\
         --genomeFastaFiles $fasta \\
         --sjdbGTFfile $gtf \\
         --runThreadN ${task.cpus} \\
         --chimSegmentMin 10 \\
         --chimOutType WithinBAM SoftClip
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'

    input:
    path index
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.bam"), emit: bam
    path "*.Log.final.out"             , emit: log

    script:
    """
    STAR --genomeDir $index \\
         --readFilesIn $reads \\
         --readFilesCommand zcat \\
         --outFileNamePrefix ${sample_id}. \\
         --outSAMtype BAM SortedByCoordinate \\
         --runThreadN ${task.cpus} \\
         --chimSegmentMin 10 \\
         --chimOutType WithinBAM SoftClip \\
         --chimJunctionOverhangMin 10 \\
         --chimScoreMin 1 \\
         --chimScoreDropMax 30 \\
         --chimScoreJunctionNonGTAG 0 \\
         --chimScoreSeparation 1 \\
         --alignSJstitchMismatchNmax 5 -1 5 5 \\
         --chimSegmentReadGapMax 3 \\
         --outSAMattributes NH HI AS NM MD SA
    """
}
