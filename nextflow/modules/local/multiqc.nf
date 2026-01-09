/*
 * MultiQC - Aggregate results from bioinformatics analyses
 */

process MULTIQC {
    tag "multiqc"
    label 'process_single'
    
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode
    
    input:
    path(multiqc_files)
    
    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data")       , emit: data
    
    script:
    """
    multiqc \\
        --force \\
        --title "VCA Pipeline Report" \\
        --comment "Variant Calling Analysis for Circulating Tumor DNA" \\
        --filename multiqc_report \\
        .
    """
    
    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html
    """
}
