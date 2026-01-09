/*
 * SnpEff - Genetic variant annotation and effect prediction
 */

process SNPEFF_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode,
        pattern: "*.ann.vcf.gz*"
    publishDir "${params.outdir}/snpeff", mode: params.publish_dir_mode,
        pattern: "*.{csv,html,txt}"
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("*.ann.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.ann.vcf.gz.tbi"), emit: tbi
    path("*.snpeff_stats.csv")               , emit: stats
    path("*.snpeff_summary.html")            , emit: html
    path("snpEff_genes.txt")                 , emit: genes, optional: true
    
    script:
    def prefix = "${meta.id}"
    def cache_arg = params.snpeff_cache ? "-dataDir ${params.snpeff_cache}" : ""
    
    """
    # Run SnpEff annotation
    snpEff \\
        -Xmx${task.memory.toGiga()}g \\
        ${params.snpeff_db} \\
        ${cache_arg} \\
        -csvStats ${prefix}.snpeff_stats.csv \\
        -stats ${prefix}.snpeff_summary.html \\
        -nodownload \\
        -canon \\
        ${vcf} \\
        > ${prefix}.ann.vcf
    
    # Compress and index
    bgzip -c ${prefix}.ann.vcf > ${prefix}.ann.vcf.gz
    tabix -p vcf ${prefix}.ann.vcf.gz
    
    # Cleanup
    rm -f ${prefix}.ann.vcf
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.ann.vcf.gz
    touch ${prefix}.ann.vcf.gz.tbi
    touch ${prefix}.snpeff_stats.csv
    touch ${prefix}.snpeff_summary.html
    """
}

process SNPEFF_DOWNLOAD {
    tag "${db}"
    label 'process_single'
    
    publishDir "${params.outdir}/reference/snpeff", mode: params.publish_dir_mode
    
    input:
    val(db)
    
    output:
    path("data/${db}"), emit: cache
    
    script:
    """
    snpEff download -v ${db}
    """
    
    stub:
    """
    mkdir -p data/${db}
    touch data/${db}/snpEffectPredictor.bin
    """
}
