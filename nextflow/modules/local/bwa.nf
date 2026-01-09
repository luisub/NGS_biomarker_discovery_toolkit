/*
 * BWA - Burrows-Wheeler Aligner
 * Alignment of sequencing reads to reference genome
 */

process BWA_INDEX {
    tag "index"
    label 'process_high'
    
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode
    
    input:
    path(fasta)
    
    output:
    path("bwa_index")    , emit: index
    path("${fasta}.fai") , emit: fai
    
    script:
    """
    # Create index directory
    mkdir -p bwa_index
    
    # Create BWA index
    bwa index -p bwa_index/genome ${fasta}
    
    # Create FASTA index
    samtools faidx ${fasta}
    """
    
    stub:
    """
    mkdir -p bwa_index
    touch bwa_index/genome.amb
    touch bwa_index/genome.ann
    touch bwa_index/genome.bwt
    touch bwa_index/genome.pac
    touch bwa_index/genome.sa
    touch ${fasta}.fai
    """
}

process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    
    publishDir "${params.outdir}/aligned", mode: params.publish_dir_mode,
        pattern: "*.bam"
    
    input:
    tuple val(meta), path(reads)
    path(index)
    path(fasta)
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    
    script:
    def prefix = "${meta.id}"
    def index_prefix = "${index}/genome"
    
    // Read group information
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.patient_id}\\tPL:ILLUMINA\\tLB:${meta.id}"
    
    if (meta.single_end) {
        """
        bwa mem \\
            -t ${task.cpus} \\
            -R '${rg}' \\
            ${index_prefix} \\
            ${reads[0]} \\
            | samtools view -@ ${task.cpus} -bS - \\
            > ${prefix}.bam
        """
    } else {
        """
        bwa mem \\
            -t ${task.cpus} \\
            -R '${rg}' \\
            ${index_prefix} \\
            ${reads[0]} \\
            ${reads[1]} \\
            | samtools view -@ ${task.cpus} -bS - \\
            > ${prefix}.bam
        """
    }
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
