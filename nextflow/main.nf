#!/usr/bin/env nextflow

/*
 * VCA Pipeline - Variant Calling Analysis for Circulating Tumor DNA
 * 
 * Author: Luis Aguilera
 * 
 * This pipeline implements:
 *   1. Quality control (FastQC)
 *   2. Read trimming (fastp)
 *   3. Alignment (BWA-MEM)
 *   4. Duplicate removal (samtools markdup)
 *   5. Variant calling (LoFreq)
 *   6. Variant annotation (SnpEff)
 *   7. Quality report aggregation (MultiQC)
 */

nextflow.enable.dsl = 2

// =============================================================================
// Help Message
// =============================================================================
def helpMessage() {
    log.info """
    ╔═══════════════════════════════════════════════════════════════════════════╗
    ║           VCA Pipeline - Variant Calling Analysis for ctDNA               ║
    ╚═══════════════════════════════════════════════════════════════════════════╝
    
    Usage:
        nextflow run main.nf --input samplesheet.csv --genome_fasta ref.fa -profile singularity
    
    Required arguments:
        --input             Path to samplesheet CSV file
        --genome_fasta      Path to reference genome FASTA file
    
    Optional arguments:
        --outdir            Output directory [default: ${params.outdir}]
        --bwa_index         Path to pre-built BWA index (skip indexing)
        --dbsnp             Path to dbSNP VCF for annotation
        --snpeff_db         SnpEff database name [default: ${params.snpeff_db}]
        --target_bed        BED file for target regions
        
    Quality filters:
        --min_read_depth        Min read depth for variant calling [default: ${params.min_read_depth}]
        --min_allele_frequency  Min VAF threshold [default: ${params.min_allele_frequency}]
        --min_base_quality      Min base quality score [default: ${params.min_base_quality}]
        --min_mapping_quality   Min mapping quality [default: ${params.min_mapping_quality}]
    
    Skip options:
        --skip_fastqc       Skip FastQC step
        --skip_trimming     Skip read trimming
        --skip_markdup      Skip duplicate marking
        --skip_annotation   Skip SnpEff annotation
        --skip_multiqc      Skip MultiQC report
    
    Execution profiles:
        -profile singularity    Run with Singularity containers (recommended for HPC)
        -profile docker         Run with Docker containers
        -profile slurm          Submit jobs to SLURM scheduler
        -profile test           Run with test dataset
    
    Examples:
        # Basic run with Singularity
        nextflow run main.nf \\
            --input samplesheet.csv \\
            --genome_fasta /path/to/GRCh38.fa \\
            -profile singularity
        
        # Run on SLURM cluster
        nextflow run main.nf \\
            --input samplesheet.csv \\
            --genome_fasta /path/to/GRCh38.fa \\
            -profile singularity,slurm
        
        # Resume failed run
        nextflow run main.nf -resume
    
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// =============================================================================
// Parameter Validation
// =============================================================================
if (!params.input) {
    error "ERROR: --input samplesheet.csv is required"
}

if (!params.genome_fasta && !params.bwa_index) {
    error "ERROR: --genome_fasta or --bwa_index is required"
}

// =============================================================================
// Log Pipeline Info
// =============================================================================
log.info """
╔═══════════════════════════════════════════════════════════════════════════╗
║                    VCA Pipeline - Execution Summary                        ║
╚═══════════════════════════════════════════════════════════════════════════╝

Input/Output:
    samplesheet     : ${params.input}
    output dir      : ${params.outdir}

Reference:
    genome FASTA    : ${params.genome_fasta ?: 'Using pre-built index'}
    BWA index       : ${params.bwa_index ?: 'Will be generated'}
    SnpEff DB       : ${params.snpeff_db}

Filters:
    min depth       : ${params.min_read_depth}
    min VAF         : ${params.min_allele_frequency}
    min base qual   : ${params.min_base_quality}
    min map qual    : ${params.min_mapping_quality}

Target region:
    gene            : ${params.target_gene}
    region          : ${params.target_chromosome}:${params.target_start}-${params.target_end}

Execution:
    profile         : ${workflow.profile}
    threads         : ${params.threads}

═══════════════════════════════════════════════════════════════════════════
"""

// =============================================================================
// Include Modules
// =============================================================================
include { FASTQC                        } from './modules/local/fastqc'
include { FASTP                         } from './modules/local/fastp'
include { BWA_INDEX                     } from './modules/local/bwa'
include { BWA_MEM                       } from './modules/local/bwa'
include { SAMTOOLS_SORT                 } from './modules/local/samtools'
include { SAMTOOLS_INDEX                } from './modules/local/samtools'
include { SAMTOOLS_MARKDUP              } from './modules/local/samtools'
include { SAMTOOLS_FLAGSTAT             } from './modules/local/samtools'
include { SAMTOOLS_FAIDX                } from './modules/local/samtools'
include { LOFREQ_CALL                   } from './modules/local/lofreq'
include { SNPEFF_ANNOTATE               } from './modules/local/snpeff'
include { BCFTOOLS_NORM                 } from './modules/local/bcftools'
include { BCFTOOLS_FILTER_GERMLINE      } from './modules/local/bcftools'
include { BCFTOOLS_STATS                } from './modules/local/bcftools'
include { MOSDEPTH                      } from './modules/local/mosdepth'
include { PICARD_COLLECTINSERTSIZEMETRICS } from './modules/local/picard'
include { AGGREGATE_VARIANTS            } from './modules/local/aggregate'
include { MULTIQC                       } from './modules/local/multiqc'

// =============================================================================
// Functions
// =============================================================================

// Parse samplesheet CSV
def parseSamplesheet(samplesheet_path) {
    Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [:]
            meta.id         = row.sample
            meta.patient_id = row.patient_id ?: row.sample
            meta.timepoint  = row.timepoint ?: 'unknown'
            meta.single_end = row.fastq_2 ? false : true
            
            def fastq_1 = file(row.fastq_1, checkIfExists: true)
            def fastq_2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : null
            
            if (meta.single_end) {
                return [meta, [fastq_1]]
            } else {
                return [meta, [fastq_1, fastq_2]]
            }
        }
}

// =============================================================================
// Main Workflow
// =============================================================================
workflow {
    
    // Create channels from input
    ch_input = parseSamplesheet(params.input)
    
    // Reference genome channel
    ch_fasta = params.genome_fasta ? Channel.fromPath(params.genome_fasta, checkIfExists: true) : Channel.empty()
    
    // ==========================================================================
    // Step 1a: Reference Genome Indexing - BWA (if needed)
    // ==========================================================================
    if (params.bwa_index) {
        // Use pre-built index
        ch_bwa_index = Channel.fromPath("${params.bwa_index}*", checkIfExists: true).collect()
    } else {
        // Build BWA index
        BWA_INDEX(ch_fasta)
        ch_bwa_index = BWA_INDEX.out.index
    }
    
    // ==========================================================================
    // Step 1b: Reference Genome Indexing - FAIDX (if needed)
    // ==========================================================================
    if (params.genome_fasta_fai) {
        ch_fai = Channel.fromPath(params.genome_fasta_fai, checkIfExists: true)
    } else {
        SAMTOOLS_FAIDX(ch_fasta)
        ch_fai = SAMTOOLS_FAIDX.out.fai
    }
    
    // ==========================================================================
    // Setup: gnomAD channel for germline filtering
    // ==========================================================================
    ch_gnomad_vcf = params.gnomad_vcf ? Channel.fromPath(params.gnomad_vcf, checkIfExists: true) : Channel.empty()
    ch_gnomad_tbi = params.gnomad_tbi ? Channel.fromPath(params.gnomad_tbi, checkIfExists: true) : Channel.empty()
    
    // Target BED file (optional)
    ch_target_bed = params.target_bed ? Channel.fromPath(params.target_bed, checkIfExists: true) : Channel.empty()
    
    // ==========================================================================
    // Step 2: Quality Control (FastQC)
    // ==========================================================================
    if (!params.skip_fastqc) {
        FASTQC(ch_input)
        ch_fastqc_reports = FASTQC.out.reports
    } else {
        ch_fastqc_reports = Channel.empty()
    }
    
    // ==========================================================================
    // Step 3: Read Trimming (fastp)
    // ==========================================================================
    if (!params.skip_trimming) {
        FASTP(ch_input)
        ch_trimmed_reads = FASTP.out.reads
        ch_fastp_reports = FASTP.out.json_only
    } else {
        ch_trimmed_reads = ch_input
        ch_fastp_reports = Channel.empty()
    }
    
    // ==========================================================================
    // Step 4: Alignment (BWA-MEM)
    // ==========================================================================
    BWA_MEM(
        ch_trimmed_reads,
        ch_bwa_index.collect(),
        ch_fasta.collect()
    )
    
    // ==========================================================================
    // Step 5: Sort BAM
    // ==========================================================================
    SAMTOOLS_SORT(BWA_MEM.out.bam)
    
    // ==========================================================================
    // Step 6: Mark Duplicates
    // ==========================================================================
    if (!params.skip_markdup) {
        SAMTOOLS_MARKDUP(SAMTOOLS_SORT.out.bam)
        ch_bam_for_calling = SAMTOOLS_MARKDUP.out.bam
    } else {
        ch_bam_for_calling = SAMTOOLS_SORT.out.bam
    }
    
    // ==========================================================================
    // Step 7: Index BAM
    // ==========================================================================
    SAMTOOLS_INDEX(ch_bam_for_calling)
    
    // Combine BAM with index
    ch_bam_indexed = ch_bam_for_calling
        .join(SAMTOOLS_INDEX.out.bai)
    
    // ==========================================================================
    // Step 8: Alignment Statistics
    // ==========================================================================
    SAMTOOLS_FLAGSTAT(ch_bam_indexed)
    
    // ==========================================================================
    // Step 8b: Coverage Analysis (mosdepth)
    // ==========================================================================
    ch_mosdepth_summary = Channel.empty()
    if (!params.skip_coverage) {
        def bed_for_mosdepth = params.target_bed ? ch_target_bed.first() : []
        MOSDEPTH(ch_bam_indexed, bed_for_mosdepth)
        ch_mosdepth_summary = MOSDEPTH.out.summary_only
    }
    
    // ==========================================================================
    // Step 8c: Insert Size Metrics (Picard)
    // ==========================================================================
    ch_insert_size_metrics = Channel.empty()
    if (!params.skip_coverage) {
        PICARD_COLLECTINSERTSIZEMETRICS(ch_bam_indexed)
        ch_insert_size_metrics = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics_only
    }
    
    // ==========================================================================
    // Step 9: Variant Calling (LoFreq)
    // ==========================================================================
    LOFREQ_CALL(
        ch_bam_indexed,
        ch_fasta.collect(),
        ch_fai.collect()
    )
    ch_lofreq_stats = LOFREQ_CALL.out.stats
    
    // ==========================================================================
    // Step 9b: Variant Normalization (bcftools)
    // ==========================================================================
    if (!params.skip_normalization) {
        BCFTOOLS_NORM(LOFREQ_CALL.out.vcf, ch_fasta.collect())
        ch_vcf_normalized = BCFTOOLS_NORM.out.vcf
        ch_vcf_normalized_tbi = BCFTOOLS_NORM.out.tbi
    } else {
        ch_vcf_normalized = LOFREQ_CALL.out.vcf
        ch_vcf_normalized_tbi = LOFREQ_CALL.out.tbi
    }
    
    // ==========================================================================
    // Step 10: Variant Annotation (SnpEff)
    // ==========================================================================
    if (!params.skip_annotation) {
        // Join VCF with its index
        ch_vcf_for_annotation = ch_vcf_normalized
            .join(ch_vcf_normalized_tbi)
        
        SNPEFF_ANNOTATE(ch_vcf_for_annotation)
        ch_vcf_annotated = SNPEFF_ANNOTATE.out.vcf
        ch_vcf_annotated_tbi = SNPEFF_ANNOTATE.out.tbi
        ch_snpeff_reports = SNPEFF_ANNOTATE.out.stats
    } else {
        ch_vcf_annotated = ch_vcf_normalized
        ch_vcf_annotated_tbi = ch_vcf_normalized_tbi
        ch_snpeff_reports = Channel.empty()
    }
    
    // ==========================================================================
    // Step 10b: Germline Filtering (bcftools with gnomAD)
    // ==========================================================================
    ch_germline_stats = Channel.empty()
    if (!params.skip_germline_filter && params.gnomad_vcf) {
        ch_vcf_for_germline = ch_vcf_annotated
            .join(ch_vcf_annotated_tbi)
        
        BCFTOOLS_FILTER_GERMLINE(
            ch_vcf_for_germline,
            ch_gnomad_vcf.first(),
            ch_gnomad_tbi.first()
        )
        ch_vcf_final = BCFTOOLS_FILTER_GERMLINE.out.vcf
        ch_germline_stats = BCFTOOLS_FILTER_GERMLINE.out.stats
    } else {
        ch_vcf_final = ch_vcf_annotated
    }
    
    // ==========================================================================
    // Step 10c: VCF Statistics
    // ==========================================================================
    BCFTOOLS_STATS(ch_vcf_final)
    ch_bcftools_stats = BCFTOOLS_STATS.out.stats_only
    
    // ==========================================================================
    // Step 11: Aggregate Variants
    // ==========================================================================
    // Extract just VCF paths (not metadata) for aggregation
    ch_vcf_paths = ch_vcf_final.map { meta, vcf -> vcf }
    AGGREGATE_VARIANTS(ch_vcf_paths.collect())
    
    // ==========================================================================
    // Step 12: MultiQC Report
    // ==========================================================================
    if (!params.skip_multiqc) {
        // Extract just paths (not metadata) for MultiQC
        ch_flagstat_paths = SAMTOOLS_FLAGSTAT.out.stats.map { meta, path -> path }
        
        ch_multiqc_files = ch_fastqc_reports
            .mix(ch_fastp_reports)
            .mix(ch_flagstat_paths)
            .mix(ch_snpeff_reports)
            .mix(ch_mosdepth_summary)
            .mix(ch_insert_size_metrics)
            .mix(ch_bcftools_stats)
            .collect()
        
        MULTIQC(ch_multiqc_files)
    }
}

// =============================================================================
// Workflow Completion
// =============================================================================
workflow.onComplete {
    log.info """
    ═══════════════════════════════════════════════════════════════════════════
                            Pipeline Completed!
    ═══════════════════════════════════════════════════════════════════════════
    
    Status      : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration    : ${workflow.duration}
    Work dir    : ${workflow.workDir}
    Output dir  : ${params.outdir}
    
    ═══════════════════════════════════════════════════════════════════════════
    """.stripIndent()
}

workflow.onError {
    log.error """
    ═══════════════════════════════════════════════════════════════════════════
                            Pipeline FAILED!
    ═══════════════════════════════════════════════════════════════════════════
    
    Error message: ${workflow.errorMessage}
    
    Check the log files in: ${params.tracedir}
    
    ═══════════════════════════════════════════════════════════════════════════
    """.stripIndent()
}
