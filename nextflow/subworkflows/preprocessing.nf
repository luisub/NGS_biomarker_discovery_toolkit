/*
 * Preprocessing Subworkflow
 * Quality control and read trimming
 */

include { FASTQC } from '../modules/local/fastqc'
include { FASTP  } from '../modules/local/fastp'

workflow PREPROCESSING {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    skip_fastqc     // boolean
    skip_trimming   // boolean

    main:
    ch_fastqc_reports = Channel.empty()
    ch_fastp_reports  = Channel.empty()
    ch_trimmed_reads  = reads

    // Run FastQC on raw reads
    if (!skip_fastqc) {
        FASTQC(reads)
        ch_fastqc_reports = FASTQC.out.reports
    }

    // Trim reads with fastp
    if (!skip_trimming) {
        FASTP(reads)
        ch_trimmed_reads  = FASTP.out.reads
        ch_fastp_reports  = FASTP.out.json
    }

    emit:
    reads   = ch_trimmed_reads      // channel: [ val(meta), [ reads ] ]
    fastqc  = ch_fastqc_reports     // channel: [ reports ]
    fastp   = ch_fastp_reports      // channel: [ json ]
}
