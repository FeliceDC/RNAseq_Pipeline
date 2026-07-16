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
        rmats_results         = RNA_SEQ_ANALYSIS.out.rmats_results
        rmats_plots           = RNA_SEQ_ANALYSIS.out.rmats_plots
        rmats_sashimi         = RNA_SEQ_ANALYSIS.out.rmats_sashimi
        majiq_results         = RNA_SEQ_ANALYSIS.out.majiq_results
        leafcutter_results    = RNA_SEQ_ANALYSIS.out.leafcutter_results
        leafcutter_plots      = RNA_SEQ_ANALYSIS.out.leafcutter_plots
        leafcutter_annotated  = RNA_SEQ_ANALYSIS.out.leafcutter_annotated
}

output {
    
    fastqc_results        { path "${params.outdir}/FastQC"; mode 'copy' }
    trimgalore_results    { path "${params.outdir}/TrimGalore"; mode 'copy' }
    star_index_results    { path "${params.outdir}/STAR/Index"; mode 'copy' }
    star_align_results    { path "${params.outdir}/STAR/Alignment"; mode 'copy' }
    featurecounts_results { path "${params.outdir}/featureCounts"; mode 'copy' }
    multiqc_results       { path "${params.outdir}/MultiQC"; mode 'copy' }
    deseq2_results        { path "${params.outdir}/DESeq2"; mode 'copy' }
    enrichr_results       { path "${params.outdir}/Enrichr"; mode 'copy' }
    immucellai_results    { path "${params.outdir}/Deconvolution/ImmucellAI_results/ImmucellAI"; mode 'copy' }
    deconvolution_plots   { path "${params.outdir}/Deconvolution/ImmucellAI_results/Plots"; mode 'copy' }
    imsig_results         { path "${params.outdir}/Deconvolution/ImSig_results/ImSig"; mode 'copy' }
    imsig_plot            { path "${params.outdir}/Deconvolution/ImSig_results/Plot"; mode 'copy'}
    arriba_fusions        { path "${params.outdir}/Fusions/Arriba"; mode 'copy' }
    arriba_discarded      { path "${params.outdir}/Fusions/Arriba/Discarded"; mode 'copy' }
    arriba_plots          { path "${params.outdir}/Fusions/Arriba/Plots"; mode 'copy' }
    rmats_results         { path "${params.outdir}/Splicing/rMATS-turbo"; mode 'copy' }
    rmats_plots           { path "${params.outdir}/Splicing/rMATS-turbo/Plots"; mode 'copy' }
    rmats_sashimi         { path "${params.outdir}/Splicing/rMATS-turbo/Plots/Sashimi"; mode 'copy' }
    majiq_results         { path "${params.outdir}/Splicing/MAJIQ"; mode 'copy' }
    leafcutter_results    { path "${params.outdir}/Splicing/LeafCutter"; mode 'copy' }
    leafcutter_plots      { path "${params.outdir}/Splicing/LeafCutter/Plots"; mode 'copy' }
    leafcutter_annotated  { path "${params.outdir}/Splicing/LeafCutter"; mode 'copy' }

}
