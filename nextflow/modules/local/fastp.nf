/*
 * fastp - Fast All-in-One Preprocessing for FASTQ files
 * Performs adapter trimming, quality filtering, and polyG/polyX trimming
 */

process FASTP {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/fastp", mode: params.publish_dir_mode,
        pattern: "*.{json,html}"
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")             , emit: json
    tuple val(meta), path("*.html")             , emit: html
    path("*.json")                              , emit: json_only
    
    script:
    def prefix = "${meta.id}"
    def adapter_arg = params.adapter_fasta ? "--adapter_fasta ${params.adapter_fasta}" : ""
    
    if (meta.single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${prefix}_trimmed.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            --thread ${task.cpus} \\
            --qualified_quality_phred ${params.min_base_quality} \\
            --length_required ${params.min_read_length} \\
            --trim_front1 ${params.trim_front} \\
            --trim_tail1 ${params.trim_tail} \\
            --detect_adapter_for_pe \\
            ${adapter_arg}
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_trimmed_1.fastq.gz \\
            --out2 ${prefix}_trimmed_2.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            --thread ${task.cpus} \\
            --qualified_quality_phred ${params.min_base_quality} \\
            --length_required ${params.min_read_length} \\
            --trim_front1 ${params.trim_front} \\
            --trim_front2 ${params.trim_front} \\
            --trim_tail1 ${params.trim_tail} \\
            --trim_tail2 ${params.trim_tail} \\
            --detect_adapter_for_pe \\
            --correction \\
            ${adapter_arg}
        """
    }
    
    stub:
    def prefix = "${meta.id}"
    if (meta.single_end) {
        """
        touch ${prefix}_trimmed.fastq.gz
        touch ${prefix}.fastp.json
        touch ${prefix}.fastp.html
        """
    } else {
        """
        touch ${prefix}_trimmed_1.fastq.gz
        touch ${prefix}_trimmed_2.fastq.gz
        touch ${prefix}.fastp.json
        touch ${prefix}.fastp.html
        """
    }
}
