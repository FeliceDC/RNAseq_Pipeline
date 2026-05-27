include { FASTQC } from '../modules/fastqc'
include {TRIMGALORE} from '../modules/trimgalore'
include { STAR_INDEX; STAR_ALIGN } from '../modules/star'
include { FEATURECOUNTS } from '../modules/subread'
include { MULTIQC } from '../modules/multiqc'
include { DESEQ2 } from '../modules/deseq2'
include { ENRICHR } from '../modules/enrichr'
include { IMMUCELLAI } from '../modules/immucellai'
include { PLOT_DECONVOLUTION } from '../modules/plot_deconvolution'
include { IMSIG } from '../modules/imsig'
include { ARRIBA } from '../modules/arriba'
include { RMATS } from '../modules/rmats'
include { RMATS_PLOT } from '../modules/rmats_plot'

workflow RNA_SEQ_ANALYSIS {
log.info "RNA-seq analysis started..."

if (params.single_end) {
ch_reads = Channel.fromPath(params.input_reads, checkIfExists: true)
                  .map {file -> tuple(file.simpleName, [file]) }
} else {
ch_reads = Channel.fromFilePairs(params.input_reads, checkIfExists: true)
}

FASTQC(ch_reads)
TRIMGALORE(ch_reads)
ch_fasta = file(params.fasta)
ch_gtf   = file(params.gtf)
STAR_INDEX(ch_fasta, ch_gtf)
STAR_ALIGN(STAR_INDEX.out.index, TRIMGALORE.out.reads)
ch_bams_raccolti = STAR_ALIGN.out.bam.map { it[1] }.collect()
FEATURECOUNTS(ch_gtf, ch_bams_raccolti)
ARRIBA(STAR_ALIGN.out.bam, ch_fasta, ch_gtf)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC.out.zip,
        TRIMGALORE.out.log,
        STAR_ALIGN.out.log,
        FEATURECOUNTS.out.summary
    )
MULTIQC( ch_multiqc_files.collect() )
DESEQ2(FEATURECOUNTS.out.counts, file(params.samplesheet))
ENRICHR(DESEQ2.out.results_tables)
IMMUCELLAI(FEATURECOUNTS.out.counts)
PLOT_DECONVOLUTION( IMMUCELLAI.out.fractions )
IMSIG(FEATURECOUNTS.out.counts)
RMATS(ch_bams_raccolti, file(params.samplesheet), ch_gtf)
RMATS_PLOT(RMATS.out.splicing_results)

emit:
        fastqc_results        = FASTQC.out.html.mix(FASTQC.out.zip)
        trimgalore_results    = TRIMGALORE.out.reads.mix(TRIMGALORE.out.log)
        star_index_results    = STAR_INDEX.out.index
        star_align_results    = STAR_ALIGN.out.bam.mix(STAR_ALIGN.out.log)
        featurecounts_results = FEATURECOUNTS.out.counts.mix(FEATURECOUNTS.out.summary)
        multiqc_results       = MULTIQC.out.report.mix(MULTIQC.out.data)
        deseq2_results        = DESEQ2.out.results_tables.mix(DESEQ2.out.results_pdf)
        enrichr_results       = ENRICHR.out.enrichr_results
        immucellai_results    = IMMUCELLAI.out.tpm_matrix.mix(IMMUCELLAI.out.fractions)
        deconvolution_plots   = PLOT_DECONVOLUTION.out.plots
        imsig_results       = IMSIG.out.results
        imsig_plot          = IMSIG.out.plot
        arriba_fusions      = ARRIBA.out.fusions
        arriba_discarded    = ARRIBA.out.discarded
        arriba_plots        = ARRIBA.out.plots
        rmats_results = RMATS.out.splicing_results
        rmats_plots = RMATS_PLOT.out.plots
}
