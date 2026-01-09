/*
 * bcftools - VCF/BCF manipulation utilities
 * Includes normalization, filtering, and annotation
 */

process BCFTOOLS_NORM {
    tag "$meta.id"
    label 'process_low'
    
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode,
        pattern: "*.norm.vcf.gz*"
    
    input:
    tuple val(meta), path(vcf)
    path(fasta)
    
    output:
    tuple val(meta), path("*.norm.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.norm.vcf.gz.tbi"), emit: tbi
    
    script:
    def prefix = "${meta.id}"
    
    """
    # Normalize variants (left-align and split multi-allelic)
    bcftools norm \\
        -f ${fasta} \\
        -m -any \\
        --check-ref w \\
        -o ${prefix}.norm.vcf.gz \\
        -O z \\
        ${vcf}
    
    # Index normalized VCF
    tabix -p vcf ${prefix}.norm.vcf.gz
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.norm.vcf.gz
    touch ${prefix}.norm.vcf.gz.tbi
    """
}

process BCFTOOLS_FILTER_GERMLINE {
    tag "$meta.id"
    label 'process_low'
    
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode,
        pattern: "*.somatic.vcf.gz*"
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path(gnomad_vcf)
    path(gnomad_tbi)
    
    output:
    tuple val(meta), path("*.somatic.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.somatic.vcf.gz.tbi"), emit: tbi
    path("*.filter_stats.txt")                   , emit: stats
    
    script:
    def prefix = "${meta.id}"
    def max_af = params.max_population_af ?: 0.01
    
    """
    # Annotate with gnomAD allele frequencies
    bcftools annotate \\
        -a ${gnomad_vcf} \\
        -c INFO/gnomAD_AF:=INFO/AF \\
        -O z \\
        -o ${prefix}.annotated.vcf.gz \\
        ${vcf}
    
    tabix -p vcf ${prefix}.annotated.vcf.gz
    
    # Filter out common germline variants
    bcftools filter \\
        -e 'INFO/gnomAD_AF > ${max_af}' \\
        -s 'GERMLINE' \\
        -m + \\
        -O z \\
        -o ${prefix}.somatic.vcf.gz \\
        ${prefix}.annotated.vcf.gz
    
    tabix -p vcf ${prefix}.somatic.vcf.gz
    
    # Generate filter statistics
    echo "=== Germline Filter Statistics ===" > ${prefix}.filter_stats.txt
    echo "Sample: ${meta.id}" >> ${prefix}.filter_stats.txt
    echo "Max population AF threshold: ${max_af}" >> ${prefix}.filter_stats.txt
    echo "" >> ${prefix}.filter_stats.txt
    echo "Variants before filtering:" >> ${prefix}.filter_stats.txt
    bcftools view -H ${vcf} | wc -l >> ${prefix}.filter_stats.txt
    echo "Variants after germline filtering (PASS):" >> ${prefix}.filter_stats.txt
    bcftools view -f PASS -H ${prefix}.somatic.vcf.gz | wc -l >> ${prefix}.filter_stats.txt
    echo "Variants filtered as GERMLINE:" >> ${prefix}.filter_stats.txt
    bcftools view -f GERMLINE -H ${prefix}.somatic.vcf.gz | wc -l >> ${prefix}.filter_stats.txt
    
    # Cleanup
    rm -f ${prefix}.annotated.vcf.gz ${prefix}.annotated.vcf.gz.tbi
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.somatic.vcf.gz
    touch ${prefix}.somatic.vcf.gz.tbi
    touch ${prefix}.filter_stats.txt
    """
}

process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'
    
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    
    publishDir "${params.outdir}/stats", mode: params.publish_dir_mode
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*.bcftools_stats.txt"), emit: stats
    path("*.bcftools_stats.txt")                 , emit: stats_only
    
    script:
    def prefix = "${meta.id}"
    
    """
    bcftools stats ${vcf} > ${prefix}.bcftools_stats.txt
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.bcftools_stats.txt
    """
}
