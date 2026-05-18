include { FASTQC } from '../modules/fastqc'
include {TRIMGALORE} from '../modules/trimgalore'
include { STAR_INDEX; STAR_ALIGN } from '../modules/star'
include { FEATURECOUNTS } from '../modules/subread'
include { MULTIQC } from '../modules/multiqc'
include { DESEQ2 } from '../modules/deseq2'
include { ENRICHR } from '../modules/enrichr'
include { IMMUCELLAI } from '../modules/immucellai'
include { PLOT_DECONVOLUTION } from '../modules/plot_deconvolution'

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
}

publish:



output {

fastqc_results {
        path FASTQC.out.html
        path FASTQC.out.zip
        mode 'copy'
    }

trimgalore_results {
        path TRIMGALORE.out.reads
        path TRIMGALORE.out.log
        mode 'copy'
    }

star_index_results {
        path STAR_INDEX.out.index
        mode 'copy' 
    }

star_align_results {
        path STAR_ALIGN.out.bam
        path STAR_ALIGN.out.log
        mode 'copy'
    }

featurecounts_results {
        path FEATURECOUNTS.out.counts
        path FEATURECOUNTS.out.summary
        mode 'copy'
    }

multiqc_results {
        path MULTIQC.out.report
        path MULTIQC.out.data
        mode 'copy'
    }

deseq2_results {
        path DESEQ2.out.results_tables
        path DESEQ2.out.results_pdf
        mode 'copy'
    }

enrichr_results {
        path ENRICHR.out.enrichr_results
        mode 'copy'
    }


immucellai_results {
        path IMMUCELLAI.out.tpm_matrix
        path IMMUCELLAI.out.fractions
        mode 'copy'
    }

deconvolution_plots {
        path PLOT_DECONVOLUTION.out.plots
        mode 'copy'
    }

}
