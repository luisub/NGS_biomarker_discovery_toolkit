/*
 * FastQC - Quality Control for High Throughput Sequencing Data
 */

process FASTQC {
    tag "$meta.id"
    label 'process_low'
    
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename -> filename.endsWith('.zip') ? "zips/${filename}" : filename }
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path("*.{html,zip}")           , emit: reports
    
    script:
    def prefix = "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ''
    """
    fastqc \\
        --threads ${task.cpus} \\
        ${memory} \\
        --outdir . \\
        ${reads}
    """
    
    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}_1_fastqc.html
    touch ${prefix}_1_fastqc.zip
    touch ${prefix}_2_fastqc.html
    touch ${prefix}_2_fastqc.zip
    """
}
