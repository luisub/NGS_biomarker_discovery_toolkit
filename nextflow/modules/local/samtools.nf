/*
 * Samtools - Utilities for SAM/BAM file manipulation
 */

process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/aligned", mode: params.publish_dir_mode,
        pattern: "*.sorted.bam"
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    
    script:
    def prefix = "${meta.id}"
    def memory_per_thread = (task.memory.toMega() / task.cpus).intValue()
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -m ${memory_per_thread}M \\
        -o ${prefix}.sorted.bam \\
        ${bam}
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.sorted.bam
    """
}

process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'
    
    publishDir "${params.outdir}/aligned", mode: params.publish_dir_mode,
        pattern: "*.bai"
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.bai"), emit: bai
    
    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
    
    stub:
    """
    touch ${bam}.bai
    """
}

process SAMTOOLS_MARKDUP {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/aligned", mode: params.publish_dir_mode,
        pattern: "*.{bam,txt}"
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.markdup.bam")    , emit: bam
    tuple val(meta), path("*.markdup_stats.txt"), emit: stats
    
    script:
    def prefix = "${meta.id}"
    """
    # First, fix mate information and add ms/MC tags
    samtools fixmate -@ ${task.cpus} -m ${bam} fixmate.bam
    
    # Sort by coordinate (required for markdup)
    samtools sort -@ ${task.cpus} -o sorted_fixmate.bam fixmate.bam
    
    # Mark duplicates
    samtools markdup \\
        -@ ${task.cpus} \\
        -f ${prefix}.markdup_stats.txt \\
        -r \\
        sorted_fixmate.bam \\
        ${prefix}.markdup.bam
    
    # Cleanup intermediate files
    rm -f fixmate.bam sorted_fixmate.bam
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.markdup.bam
    touch ${prefix}.markdup_stats.txt
    """
}

process SAMTOOLS_FLAGSTAT {
    tag "$meta.id"
    label 'process_single'
    
    publishDir "${params.outdir}/stats", mode: params.publish_dir_mode
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*.flagstat.txt"), emit: stats
    
    script:
    def prefix = "${meta.id}"
    """
    samtools flagstat \\
        -@ ${task.cpus} \\
        ${bam} \\
        > ${prefix}.flagstat.txt
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.flagstat.txt
    """
}

process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'
    
    publishDir "${params.outdir}/stats", mode: params.publish_dir_mode
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*.stats.txt"), emit: stats
    
    script:
    def prefix = "${meta.id}"
    """
    samtools stats \\
        -@ ${task.cpus} \\
        ${bam} \\
        > ${prefix}.stats.txt
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.stats.txt
    """
}

process SAMTOOLS_FAIDX {
    tag "faidx"
    label 'process_single'
    
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode,
        pattern: "*.fai"
    
    input:
    path(fasta)
    
    output:
    path("*.fai"), emit: fai
    
    script:
    """
    samtools faidx ${fasta}
    """
    
    stub:
    """
    touch ${fasta}.fai
    """
}
