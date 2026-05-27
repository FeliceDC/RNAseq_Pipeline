nextflow.enable.dsl=2

if (!params.design) {
    exit 1, """
    ❌ WAIT! The analysis cannot start: '--design' missing!
You must specify the column(s) of the samplesheet to use for differential analysis.
    
    Example:
    nextflow run main.nf --design "condition" -resume
    nextflow run main.nf --design "treatment + age" -profile test
    =======================================================
    """
}

if (!params.input_reads || !params.fasta || !params.gtf) {
    exit 1, """
    ❌ WAIT!: Missing input parameters!
    
    You must specify the file paths
    
    Required parameters:
      --input_reads   : Path to FASTQ files (es. "data/*_{1,2}.fastq.gz")
      --fasta         : Reference genome in FASTA format (es. "genome.fa")
      --gtf           : Annotation file in GTF format (es. "annotation.gtf")
      
    Usage example:
    nextflow run main.nf \\
      --input_reads "my_data/*_{1,2}.fastq.gz" \\
      --fasta "ref/genome.fa" \\
      --gtf "ref/annotation.gtf" \\
      --design "condition" \\
      --samplesheet "my_samplesheet.csv"
    =======================================================
    """
}

workflow.onComplete {
    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    println msg

    if (workflow.success) {
        println "✅ Analysis completed successfully! The results are in: ${params.outdir}"
    } else {
        println "❌ Oops... the pipeline stopped due to an error."
    }
}

include {RNA_SEQ_ANALYSIS} from './workflows/rnaseq_pipeline'

workflow {
RNA_SEQ_ANALYSIS()

publish:
        fastqc_results        = RNA_SEQ_ANALYSIS.out.fastqc_results
        trimgalore_results    = RNA_SEQ_ANALYSIS.out.trimgalore_results
        star_index_results    = RNA_SEQ_ANALYSIS.out.star_index_results
        star_align_results    = RNA_SEQ_ANALYSIS.out.star_align_results
        featurecounts_results = RNA_SEQ_ANALYSIS.out.featurecounts_results
        multiqc_results       = RNA_SEQ_ANALYSIS.out.multiqc_results
        deseq2_results        = RNA_SEQ_ANALYSIS.out.deseq2_results
        enrichr_results       = RNA_SEQ_ANALYSIS.out.enrichr_results
        immucellai_results    = RNA_SEQ_ANALYSIS.out.immucellai_results
        deconvolution_plots   = RNA_SEQ_ANALYSIS.out.deconvolution_plots
        imsig_results       = RNA_SEQ_ANALYSIS.out.imsig_results
        imsig_plot          = RNA_SEQ_ANALYSIS.out.imsig_plot
        arriba_fusions      = RNA_SEQ_ANALYSIS.out.arriba_fusions
        arriba_discarded    = RNA_SEQ_ANALYSIS.out.arriba_discarded
        arriba_plots        = RNA_SEQ_ANALYSIS.out.arriba_plots
        rmats_results = RNA_SEQ_ANALYSIS.out.rmats_results
        rmats_plots = RNA_SEQ_ANALYSIS.out.rmats_plots

}

output {
    fastqc_results        { path "fastqc"; mode 'copy' }
    trimgalore_results    { path "trimgalore"; mode 'copy' }
    star_index_results    { path "star/index"; mode 'copy' }
    star_align_results    { path "star/alignment"; mode 'copy' }
    featurecounts_results { path "featureCounts"; mode 'copy' }
    multiqc_results       { path "multiqc"; mode 'copy' }
    deseq2_results        { path "deseq2"; mode 'copy' }
    enrichr_results       { path "enrichr"; mode 'copy' }
    immucellai_results    { path "deconvolution/immucellai_results/immucellai"; mode 'copy' }
    deconvolution_plots   { path "deconvolution/immucellai_results/plots"; mode 'copy' }
    imsig_results       { path "deconvolution/imsig_results/imsig"; mode 'copy' }
    imsig_plot          { path "deconvolution/imsig_results/plot"; mode 'copy'}
    arriba_fusions      { path "fusions/arriba"; mode 'copy' }
    arriba_discarded    { path "fusions/arriba/discarded"; mode 'copy' }
    arriba_plots        { path "fusions/arriba/plots"; mode 'copy' }
    rmats_results { path "splicing/rmats"; mode 'copy' }
    rmats_plots { path "splicing/plots"; mode 'copy' }
}
