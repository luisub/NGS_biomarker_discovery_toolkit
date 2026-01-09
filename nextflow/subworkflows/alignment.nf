/*
 * Alignment Subworkflow
 * BWA-MEM alignment and BAM processing
 */

include { BWA_MEM          } from '../modules/local/bwa'
include { SAMTOOLS_SORT    } from '../modules/local/samtools'
include { SAMTOOLS_INDEX   } from '../modules/local/samtools'
include { SAMTOOLS_MARKDUP } from '../modules/local/samtools'
include { SAMTOOLS_FLAGSTAT} from '../modules/local/samtools'

workflow ALIGNMENT {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    bwa_index       // channel: [ index_files ]
    fasta           // channel: [ fasta ]
    skip_markdup    // boolean

    main:
    // Align reads with BWA-MEM
    BWA_MEM(reads, bwa_index, fasta)

    // Sort BAM
    SAMTOOLS_SORT(BWA_MEM.out.bam)

    // Mark duplicates (optional)
    if (!skip_markdup) {
        SAMTOOLS_MARKDUP(SAMTOOLS_SORT.out.bam)
        ch_bam = SAMTOOLS_MARKDUP.out.bam
    } else {
        ch_bam = SAMTOOLS_SORT.out.bam
    }

    // Index BAM
    SAMTOOLS_INDEX(ch_bam)

    // Combine BAM with index
    ch_bam_indexed = ch_bam.join(SAMTOOLS_INDEX.out.bai)

    // Generate alignment statistics
    SAMTOOLS_FLAGSTAT(ch_bam_indexed)

    emit:
    bam       = ch_bam              // channel: [ val(meta), bam ]
    bam_bai   = ch_bam_indexed      // channel: [ val(meta), bam, bai ]
    stats     = SAMTOOLS_FLAGSTAT.out.stats  // channel: [ val(meta), stats ]
}
