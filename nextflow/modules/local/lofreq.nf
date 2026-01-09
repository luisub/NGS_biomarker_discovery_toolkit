/*
 * LoFreq - Sensitive variant calling from sequencing data
 * Designed for low-frequency variant detection in ctDNA
 * 
 * Uses mulled container with LoFreq + samtools + htslib for all dependencies
 */

/*
 * LOFREQ_VITERBI - Base Quality Score Recalibration
 * Uses Viterbi algorithm to recalibrate base quality scores
 * Improves accuracy of low-frequency variant detection
 */
process LOFREQ_VITERBI {
    tag "$meta.id"
    label 'process_high'
    
    // Mulled container with lofreq, samtools, and htslib
    container 'quay.io/biocontainers/mulled-v2-4f8e2ed02c274b738b62a4e99cb6e02c84be1a6c:8d8db227fb4e07410c6e7fe07be0eae17b859ee4-0'
    
    publishDir "${params.outdir}/aligned", mode: params.publish_dir_mode,
        pattern: "*.viterbi.bam*", enabled: params.save_bqsr_bam
    
    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    
    output:
    tuple val(meta), path("*.viterbi.bam"), path("*.viterbi.bam.bai"), emit: bam
    path("*.viterbi_stats.txt"), emit: stats
    
    script:
    def prefix = "${meta.id}"
    """
    echo "=== LoFreq Viterbi BQSR ===" > ${prefix}.viterbi_stats.txt
    echo "Sample: ${meta.id}" >> ${prefix}.viterbi_stats.txt
    echo "Input BAM: ${bam}" >> ${prefix}.viterbi_stats.txt
    echo "" >> ${prefix}.viterbi_stats.txt
    
    # Get base quality distribution before BQSR
    echo "Base quality distribution BEFORE BQSR:" >> ${prefix}.viterbi_stats.txt
    samtools view ${bam} | head -10000 | cut -f11 | fold -w1 | sort | uniq -c | sort -rn | head -5 >> ${prefix}.viterbi_stats.txt
    echo "" >> ${prefix}.viterbi_stats.txt
    
    # Run LoFreq Viterbi for base quality recalibration
    # This uses a Hidden Markov Model to recalibrate base qualities
    lofreq viterbi \\
        -f ${fasta} \\
        -o ${prefix}.viterbi.bam \\
        ${bam}
    
    # Index recalibrated BAM
    samtools index ${prefix}.viterbi.bam
    
    # Get base quality distribution after BQSR
    echo "Base quality distribution AFTER BQSR:" >> ${prefix}.viterbi_stats.txt
    samtools view ${prefix}.viterbi.bam | head -10000 | cut -f11 | fold -w1 | sort | uniq -c | sort -rn | head -5 >> ${prefix}.viterbi_stats.txt
    echo "" >> ${prefix}.viterbi_stats.txt
    
    # Alignment statistics
    echo "Alignment statistics:" >> ${prefix}.viterbi_stats.txt
    samtools flagstat ${prefix}.viterbi.bam >> ${prefix}.viterbi_stats.txt
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.viterbi.bam
    touch ${prefix}.viterbi.bam.bai
    touch ${prefix}.viterbi_stats.txt
    """
}


process LOFREQ_CALL {
    tag "$meta.id"
    label 'process_high'
    
    // Mulled container with lofreq, samtools, and htslib (tabix/bgzip)
    container 'quay.io/biocontainers/mulled-v2-4f8e2ed02c274b738b62a4e99cb6e02c84be1a6c:8d8db227fb4e07410c6e7fe07be0eae17b859ee4-0'
    
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode
    
    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    
    output:
    tuple val(meta), path("*.lofreq.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.lofreq.vcf.gz.tbi"), emit: tbi
    path("*.lofreq_stats.txt")                  , emit: stats
    
    script:
    def prefix = "${meta.id}"
    def region_arg = params.target_bed ? "--bed ${params.target_bed}" : 
                     (params.whole_genome ? "" : "--region ${params.target_chromosome}:${params.target_start}-${params.target_end}")
    def sb_arg = params.max_strand_bias ? "--sb-thresh ${params.max_strand_bias}" : ""
    
    """
    # Add indel qualities (required for LoFreq)
    lofreq indelqual \\
        --dindel \\
        -f ${fasta} \\
        -o ${prefix}.indelqual.bam \\
        ${bam}
    
    # Index the new BAM
    samtools index ${prefix}.indelqual.bam
    
    # Call variants with LoFreq
    lofreq call-parallel \\
        --pp-threads ${task.cpus} \\
        --call-indels \\
        -f ${fasta} \\
        ${region_arg} \\
        ${sb_arg} \\
        -q ${params.min_base_quality} \\
        -Q ${params.min_base_quality} \\
        -m ${params.min_mapping_quality} \\
        -C ${params.min_read_depth} \\
        -o ${prefix}.lofreq.raw.vcf \\
        ${prefix}.indelqual.bam
    
    # Filter by allele frequency
    lofreq filter \\
        -i ${prefix}.lofreq.raw.vcf \\
        -o ${prefix}.lofreq.vcf \\
        --af-min ${params.min_allele_frequency}
    
    # Compress and index
    bgzip -c ${prefix}.lofreq.vcf > ${prefix}.lofreq.vcf.gz
    tabix -p vcf ${prefix}.lofreq.vcf.gz
    
    # Generate statistics
    echo "=== LoFreq Variant Calling Statistics ===" > ${prefix}.lofreq_stats.txt
    echo "Sample: ${meta.id}" >> ${prefix}.lofreq_stats.txt
    echo "Min depth: ${params.min_read_depth}" >> ${prefix}.lofreq_stats.txt
    echo "Min VAF: ${params.min_allele_frequency}" >> ${prefix}.lofreq_stats.txt
    echo "Min base quality: ${params.min_base_quality}" >> ${prefix}.lofreq_stats.txt
    echo "Min mapping quality: ${params.min_mapping_quality}" >> ${prefix}.lofreq_stats.txt
    echo "" >> ${prefix}.lofreq_stats.txt
    echo "Raw variants called: \$(grep -v '^#' ${prefix}.lofreq.raw.vcf | wc -l)" >> ${prefix}.lofreq_stats.txt
    echo "Variants after AF filter: \$(zcat ${prefix}.lofreq.vcf.gz | grep -v '^#' | wc -l)" >> ${prefix}.lofreq_stats.txt
    echo "" >> ${prefix}.lofreq_stats.txt
    echo "SNVs: \$(zcat ${prefix}.lofreq.vcf.gz | grep -v '^#' | awk 'length(\$4)==1 && length(\$5)==1' | wc -l)" >> ${prefix}.lofreq_stats.txt
    echo "Indels: \$(zcat ${prefix}.lofreq.vcf.gz | grep -v '^#' | awk 'length(\$4)!=length(\$5)' | wc -l)" >> ${prefix}.lofreq_stats.txt
    
    # Cleanup intermediate files
    rm -f ${prefix}.indelqual.bam ${prefix}.indelqual.bam.bai
    rm -f ${prefix}.lofreq.raw.vcf ${prefix}.lofreq.vcf
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.lofreq.vcf.gz
    touch ${prefix}.lofreq.vcf.gz.tbi
    touch ${prefix}.lofreq_stats.txt
    """
}
