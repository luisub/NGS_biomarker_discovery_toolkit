/*
 * mosdepth - Fast BAM/CRAM depth calculation
 * Generates coverage statistics for quality assessment
 */

process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/mosdepth:0.3.5--hd299d5a_0'
    
    publishDir "${params.outdir}/coverage", mode: params.publish_dir_mode
    
    input:
    tuple val(meta), path(bam), path(bai)
    path(target_bed)
    
    output:
    tuple val(meta), path("*.mosdepth.global.dist.txt")  , emit: global_dist
    tuple val(meta), path("*.mosdepth.summary.txt")      , emit: summary
    tuple val(meta), path("*.per-base.bed.gz")           , emit: per_base, optional: true
    tuple val(meta), path("*.regions.bed.gz")            , emit: regions, optional: true
    path("*.mosdepth.summary.txt")                       , emit: summary_only
    
    script:
    def prefix = "${meta.id}"
    def bed_arg = target_bed ? "--by ${target_bed}" : ""
    def thresholds = "1,10,20,50,100,200,500,1000"
    
    """
    mosdepth \\
        --threads ${task.cpus} \\
        --fast-mode \\
        --thresholds ${thresholds} \\
        ${bed_arg} \\
        ${prefix} \\
        ${bam}
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.mosdepth.global.dist.txt
    touch ${prefix}.mosdepth.summary.txt
    touch ${prefix}.per-base.bed.gz
    """
}
