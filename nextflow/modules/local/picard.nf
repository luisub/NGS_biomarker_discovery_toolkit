/*
 * Picard - BAM/SAM manipulation and metrics collection
 */

process PICARD_COLLECTINSERTSIZEMETRICS {
    tag "$meta.id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'
    
    publishDir "${params.outdir}/stats", mode: params.publish_dir_mode
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*.insert_size_metrics.txt"), emit: metrics
    tuple val(meta), path("*.insert_size_histogram.pdf"), emit: histogram
    path("*.insert_size_metrics.txt")                  , emit: metrics_only
    
    script:
    def prefix = "${meta.id}"
    def mem_gb = task.memory ? task.memory.toGiga() : 8
    
    """
    picard CollectInsertSizeMetrics \\
        -Xmx${mem_gb}g \\
        -I ${bam} \\
        -O ${prefix}.insert_size_metrics.txt \\
        -H ${prefix}.insert_size_histogram.pdf \\
        -M 0.5 \\
        --VALIDATION_STRINGENCY LENIENT
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.insert_size_metrics.txt
    touch ${prefix}.insert_size_histogram.pdf
    """
}

process PICARD_COLLECTWGSMETRICS {
    tag "$meta.id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'
    
    publishDir "${params.outdir}/stats", mode: params.publish_dir_mode
    
    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    
    output:
    tuple val(meta), path("*.wgs_metrics.txt"), emit: metrics
    path("*.wgs_metrics.txt")                 , emit: metrics_only
    
    script:
    def prefix = "${meta.id}"
    def mem_gb = task.memory ? task.memory.toGiga() : 8
    def min_cov = params.min_read_depth ?: 20
    
    """
    picard CollectWgsMetrics \\
        -Xmx${mem_gb}g \\
        -I ${bam} \\
        -O ${prefix}.wgs_metrics.txt \\
        -R ${fasta} \\
        --MINIMUM_MAPPING_QUALITY ${params.min_mapping_quality} \\
        --MINIMUM_BASE_QUALITY ${params.min_base_quality} \\
        --COVERAGE_CAP 10000 \\
        --VALIDATION_STRINGENCY LENIENT
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.wgs_metrics.txt
    """
}
