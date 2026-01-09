/*
 * Variant Calling Subworkflow
 * LoFreq variant calling and SnpEff annotation
 */

include { LOFREQ_CALL      } from '../modules/local/lofreq'
include { SNPEFF_ANNOTATE  } from '../modules/local/snpeff'

workflow VARIANT_CALLING {
    take:
    bam_bai         // channel: [ val(meta), bam, bai ]
    fasta           // channel: [ fasta ]
    skip_annotation // boolean

    main:
    // Call variants with LoFreq
    LOFREQ_CALL(bam_bai, fasta)

    // Annotate variants with SnpEff (optional)
    if (!skip_annotation) {
        SNPEFF_ANNOTATE(LOFREQ_CALL.out.vcf.join(LOFREQ_CALL.out.tbi))
        ch_vcf = SNPEFF_ANNOTATE.out.vcf
        ch_stats = SNPEFF_ANNOTATE.out.stats
    } else {
        ch_vcf = LOFREQ_CALL.out.vcf
        ch_stats = Channel.empty()
    }

    emit:
    vcf     = ch_vcf        // channel: [ val(meta), vcf.gz ]
    stats   = ch_stats      // channel: [ stats ]
}
