# Nextflow Variant Calling Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/singularity-compatible-blue.svg)](https://sylabs.io/singularity/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Version 1.1.0** | Variant Calling Analysis for Circulating Tumor DNA

This directory contains a Nextflow implementation of the VCA pipeline optimized for low-frequency variant detection in ctDNA samples.

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start-after-implementation)
- [Pipeline Steps](#pipeline-steps)
- [Documentation](#documentation)
- [Development Status](#todo-list)
- [Directory Structure](#proposed-directory-structure)
- [Samplesheet Format](#samplesheet-format)
- [Future Improvements](#future-improvements)
- [References](#references)

---

## Overview

The Nextflow pipeline provides:

- **Reproducibility**: Containerized execution with Docker/Singularity
- **Scalability**: Parallel processing across samples and HPC/cloud deployment
- **Portability**: Run on local machines, clusters (SLURM, PBS), or cloud (AWS, GCP, Azure)
- **Resume capability**: Automatic checkpointing and restart from failures

## Pipeline Steps

```text
┌─────────────────────────────────────────────────────────────────────────────┐
│                         VCA Pipeline Workflow                                │
└─────────────────────────────────────────────────────────────────────────────┘

  FASTQ files ──► FASTQC ──► FASTP ──► BWA-MEM ──► Sort ──► MarkDup ──► Index
                    │           │                              │
                    │           │                              ▼
                    │           │                        ┌──────────┐
                    │           │                        │ MOSDEPTH │ Coverage
                    │           │                        │ PICARD   │ Insert Size
                    │           │                        └────┬─────┘
                    │           │                             │
                    ▼           ▼                             ▼
              ┌─────────────────────────────────────────────────────────┐
              │                      LOFREQ                              │
              │              (Low-frequency variant calling)             │
              └───────────────────────┬─────────────────────────────────┘
                                      │
                                      ▼
                              ┌───────────────┐
                              │ BCFTOOLS_NORM │ Normalize variants
                              └───────┬───────┘
                                      │
                                      ▼
                              ┌───────────────┐
                              │    SNPEFF     │ Annotate variants
                              └───────┬───────┘
                                      │
                                      ▼
                         ┌────────────────────────┐
                         │ BCFTOOLS_FILTER_GERMLINE│ (Optional)
                         │    gnomAD filtering     │
                         └───────────┬────────────┘
                                     │
                                     ▼
                           ┌─────────────────┐
                           │    AGGREGATE    │ Combine all samples
                           └────────┬────────┘
                                    │
                                    ▼
                              ┌──────────┐
                              │ MULTIQC  │ Final QC Report
                              └──────────┘
```

## Documentation

| Document | Description |
|----------|-------------|
| [PARAMETERS.md](docs/PARAMETERS.md) | Complete parameter reference |
| [OUTPUT.md](docs/OUTPUT.md) | Output files description |
| [troubleshooting.md](docs/troubleshooting.md) | Common issues and solutions |
| [docker/README.md](docker/README.md) | Container build instructions |
| [CHANGELOG.md](CHANGELOG.md) | Version history |

---

## TODO List

### Phase 1: Core Pipeline Structure ✅

- [x] **1.1** Create `main.nf` - Main Nextflow workflow file
  - [x] Define input channels for FASTQ files
  - [x] Define input channels for reference genome
  - [x] Set up configuration parameters

- [x] **1.2** Create `nextflow.config` - Configuration file
  - [x] Define default parameters (threads, memory, genome paths)
  - [x] Configure executor profiles (local, slurm, cloud)
  - [x] Set up container definitions (Docker/Singularity)

- [x] **1.3** Create `modules/` directory for process definitions
  - [x] `fastqc.nf` - Quality control module
  - [x] `fastp.nf` - Read trimming module
  - [x] `bwa.nf` - Alignment module
  - [x] `samtools.nf` - BAM processing module
  - [x] `lofreq.nf` - Variant calling module
  - [x] `snpeff.nf` - Variant annotation module

### Phase 2: Process Implementation ✅

- [x] **2.1** Implement `BWA_INDEX` process
  - [x] Download GRCh38 reference genome (handled externally)
  - [x] Create BWA index
  - [x] Create FASTA index (.fai)

- [x] **2.2** Implement `FASTQC` process
  - [x] Input: Raw FASTQ files
  - [x] Output: HTML quality reports

- [x] **2.3** Implement `FASTP` process
  - [x] Input: Raw FASTQ files
  - [x] Output: Trimmed FASTQ files + JSON report
  - [x] Parameters: quality threshold, adapter sequences

- [x] **2.4** Implement `BWA_MEM` process
  - [x] Input: Trimmed FASTQ + Reference index
  - [x] Output: BAM file
  - [x] Parameters: threads, read group info

- [x] **2.5** Implement `SAMTOOLS_MARKDUP` process
  - [x] Input: Sorted BAM
  - [x] Output: Deduplicated BAM + metrics
  - [x] Tool: samtools markdup

- [x] **2.6** Implement `LOFREQ_CALL` process
  - [x] Input: Deduplicated BAM + Reference
  - [x] Output: VCF file
  - [x] Parameters: min depth, min allele frequency

- [x] **2.7** Implement `SNPEFF_ANNOTATE` process
  - [x] Input: VCF file
  - [x] Output: Annotated VCF + summary HTML
  - [x] Database: GRCh38.p14

### Phase 3: Workflow Integration ✅

- [x] **3.1** Create main workflow in `main.nf`
  - [x] Chain all processes in correct order
  - [x] Handle sample metadata (patient ID, timepoint)
  - [x] Implement conditional execution (skip options)

- [x] **3.2** Create `subworkflows/` for reusable components
  - [x] `preprocessing.nf` - QC + trimming
  - [x] `alignment.nf` - BWA + BAM processing
  - [x] `variant_calling.nf` - Lofreq + annotation

- [x] **3.3** Implement input validation
  - [x] Samplesheet CSV parser
  - [x] Reference genome validation
  - [x] Parameter validation

### Phase 4: Configuration & Execution Profiles ✅

- [x] **4.1** Configure execution profiles
  - [x] `standard` - Local execution
  - [x] `docker` - Docker containers
  - [x] `singularity` - Singularity containers (for HPC)
  - [x] `slurm` - SLURM cluster submission
  - [x] `awsbatch` - AWS Batch execution

- [x] **4.2** Create container definitions
  - [x] `Dockerfile` for pipeline tools
  - [ ] Push to Docker Hub or GitHub Container Registry

- [x] **4.3** Set resource requirements per process
  - [x] CPU allocation (process labels)
  - [x] Memory limits (process labels)
  - [x] Time limits (process labels)

### Phase 5: Output & Reporting ✅

- [x] **5.1** Implement MultiQC integration
  - [x] Aggregate all QC reports
  - [x] Generate summary HTML report

- [x] **5.2** Create results aggregation process
  - [x] Merge variant calls across samples
  - [x] Generate `all_variants.csv`
  - [x] Calculate summary statistics

- [x] **5.3** Implement pipeline reporting
  - [x] Execution timeline (built-in)
  - [x] Resource usage (trace file)
  - [x] Error logging (onError handler)

### Phase 6: Testing & Documentation ✅

- [x] **6.1** Create test dataset
  - [ ] Subset of PRJNA714799 data (use real data when available)
  - [x] Small reference region for quick testing

- [x] **6.2** Write test workflow
  - [x] `test` profile configuration
  - [x] CI/CD integration (GitHub Actions)

- [x] **6.3** Complete documentation
  - [x] Installation instructions
  - [x] Usage examples
  - [x] Parameter reference (--help)
  - [x] Troubleshooting guide

### Phase 7: Advanced Analysis (v1.1.0) ✅

- [x] **7.1** Coverage Analysis
  - [x] `mosdepth.nf` - Depth of coverage statistics
  - [x] `picard.nf` - Insert size metrics
  - [x] Integration with MultiQC

- [x] **7.2** Variant Processing
  - [x] `bcftools.nf` - Variant normalization (left-align, split multi-allelic)
  - [x] `bcftools.nf` - VCF statistics
  - [x] `SAMTOOLS_FAIDX` - Dedicated reference indexing

- [x] **7.3** Germline Filtering
  - [x] gnomAD annotation support
  - [x] Population frequency filtering
  - [x] Filter statistics reporting

- [x] **7.4** Enhanced Documentation
  - [x] `CHANGELOG.md` - Version history
  - [x] `docs/PARAMETERS.md` - Complete parameter reference
  - [x] `docs/OUTPUT.md` - Output files documentation

---

## Proposed Directory Structure

```text
nextflow/
├── main.nf                     # Main workflow orchestrator
├── nextflow.config             # Global configuration
├── CHANGELOG.md                # Version history
├── README.md                   # This file
│
├── modules/local/              # Process definitions
│   ├── fastqc.nf               # Quality control
│   ├── fastp.nf                # Read trimming
│   ├── bwa.nf                  # BWA indexing & alignment
│   ├── samtools.nf             # BAM processing & FAIDX
│   ├── lofreq.nf               # Variant calling
│   ├── snpeff.nf               # Variant annotation
│   ├── bcftools.nf             # VCF normalization & filtering (v1.1.0)
│   ├── mosdepth.nf             # Coverage analysis (v1.1.0)
│   ├── picard.nf               # Insert size metrics (v1.1.0)
│   ├── aggregate.nf            # VCF aggregation
│   └── multiqc.nf              # QC report aggregation
│
├── subworkflows/               # Reusable subworkflows
│   ├── preprocessing.nf        # FASTQC + FASTP
│   ├── alignment.nf            # BWA + sorting + dedup
│   └── variant_calling.nf      # LoFreq + SnpEff
│
├── conf/                       # Configuration profiles
│   └── singularity.config      # Singularity/HPC settings
│
├── docs/                       # Documentation
│   ├── PARAMETERS.md           # Complete parameter reference
│   ├── OUTPUT.md               # Output files description
│   └── troubleshooting.md      # Common issues & solutions
│
├── docker/                     # Container definitions
│   ├── Dockerfile              # Custom container build
│   ├── environment.yml         # Tool versions
│   └── README.md               # Build instructions
│
├── assets/                     # Pipeline assets
│   └── samplesheet.csv         # Example samplesheet
│
├── test/                       # Test configuration
│   ├── test.config             # Test profile settings
│   ├── test_samplesheet.csv    # Test input
│   └── data/                   # Small test files
│
└── .gitignore                  # Ignore work/results
```

---

## Quick Start (After Implementation)

```bash
# Run with Docker (recommended)
nextflow run main.nf \
    --input samplesheet.csv \
    --genome GRCh38 \
    --outdir results \
    -profile docker

# Run on SLURM cluster
nextflow run main.nf \
    --input samplesheet.csv \
    --genome GRCh38 \
    --outdir results \
    -profile singularity,slurm

# Resume failed run
nextflow run main.nf -resume
```

---

## Samplesheet Format

```csv
sample,fastq_1,fastq_2,patient_id,timepoint
SRR13948001,/path/to/SRR13948001_1.fastq.gz,/path/to/SRR13948001_2.fastq.gz,CTC030,baseline
SRR13948002,/path/to/SRR13948002_1.fastq.gz,/path/to/SRR13948002_2.fastq.gz,CTC030,follow_up_1
```

---

## Future Improvements

This section documents recommended enhancements identified during expert bioinformatics review. These improvements would elevate the pipeline from a functional prototype to a clinical-grade ctDNA analysis tool.

### Priority 1: Critical for Clinical Use

#### 1.1 Add Base Quality Score Recalibration (BQSR)

**Current Gap**: The pipeline does not perform base quality score recalibration before variant calling.

**Scientific Impact**: Systematic sequencing errors (e.g., from specific sequence contexts or flow cell positions) can be miscalled as variants. This is especially critical for ctDNA analysis where we're detecting variants at 0.1-1% VAF—well within the range of sequencing error rates.

**Recommended Implementation**:

```text
Current:  FASTQ → Trim → Align → Markdup → LoFreq
Proposed: FASTQ → Trim → Align → Markdup → BQSR → LoFreq
```

**Options**:

- **GATK BaseRecalibrator**: Gold standard, requires known sites VCF (dbSNP, Mills indels)
- **LoFreq viterbi**: Realignment-based quality correction, integrated with LoFreq
- **abra2**: Assembly-based realignment for better indel calling

**New Process Required**:

```groovy
process GATK_BQSR {
    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(known_sites_vcf)
    
    output:
    tuple val(meta), path("*.recal.bam"), emit: bam
}
```

---

#### 1.2 Implement Germline Variant Filtering

**Current Gap**: The pipeline calls all variants but does not distinguish somatic mutations from germline polymorphisms.

**Scientific Impact**: Without germline filtering, common SNPs (present in ~1% or more of the population) will be reported as potential cancer biomarkers, leading to false positives. For clinical ctDNA monitoring, this is unacceptable.

**Recommended Approach**:

1. **Population Frequency Filtering**:
   - Annotate with gnomAD (Genome Aggregation Database)
   - Filter variants with population AF > 0.01 (1%)
   - Add 1000 Genomes, ExAC annotations

2. **Matched Normal Comparison** (if available):
   - Compare ctDNA calls against PBMC/buffy coat baseline
   - Remove shared variants (germline)
   - Requires samplesheet extension: `normal_bam` column

3. **Panel of Normals (PoN)**:
   - Build PoN from healthy control samples
   - Filter recurrent technical artifacts

**New Processes Required**:

```groovy
process ANNOTATE_GNOMAD {
    // bcftools annotate -a gnomad.vcf.gz -c AF
}

process FILTER_GERMLINE {
    // bcftools filter -e 'gnomAD_AF > 0.01 || PoN_COUNT > 2'
}
```

**New Parameters**:

```groovy
params.gnomad_vcf = null           // Path to gnomAD VCF
params.panel_of_normals = null     // Path to PoN VCF
params.max_population_af = 0.01    // Max gnomAD AF to keep
```

---

#### 1.3 Fix Container Dependencies in LoFreq Module

**Current Issue**: The `LOFREQ_CALL` process uses `samtools index` but the LoFreq container may not include samtools, causing runtime failures on some systems.

**Current Code** (`modules/local/lofreq.nf` line 39):

```bash
samtools index ${prefix}.indelqual.bam  # May fail!
```

**Solutions**:

1. **Use a mulled container** with both LoFreq and samtools
2. **Stage samtools separately** using Nextflow's container directive
3. **Use LoFreq's built-in indexing** (if available in version)

**Recommended Fix**:

```groovy
withName: 'LOFREQ_CALL' {
    container = 'quay.io/biocontainers/mulled-v2-lofreq-samtools:...'
}
```

---

#### 1.4 Add Input Validation with nf-schema

**Current Gap**: Samplesheet parsing has edge cases that can cause silent failures.

**Known Issues**:

- Empty string `""` vs null handling for `fastq_2`
- Missing required columns not caught early
- Invalid file paths discovered late in pipeline

**Recommended Implementation**:

```groovy
// Install: nf-core/nf-schema plugin
include { validateParameters } from 'plugin/nf-schema'

workflow {
    validateParameters()
    // ... rest of workflow
}
```

**Samplesheet Schema** (`assets/schema_input.json`):

```json
{
  "properties": {
    "sample": {"type": "string", "pattern": "^[a-zA-Z0-9_]+$"},
    "fastq_1": {"type": "string", "format": "file-path", "exists": true},
    "fastq_2": {"type": "string", "format": "file-path"},
    "patient_id": {"type": "string"},
    "timepoint": {"type": "string", "enum": ["baseline", "follow_up", "progression"]}
  },
  "required": ["sample", "fastq_1"]
}
```

---

### Priority 2: Important Enhancements

#### 2.1 Add Depth of Coverage Analysis

**Current Gap**: No coverage metrics are generated to assess sequencing depth across target regions.

**Scientific Impact**: Low coverage regions produce unreliable variant calls. For ctDNA, minimum 500-1000x coverage is typically required.

**Recommended Tool**: `mosdepth` (fast and resource-efficient)

**New Process**:

```groovy
process MOSDEPTH {
    container = 'quay.io/biocontainers/mosdepth:0.3.5--hd299d5a_0'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path(target_bed)
    
    output:
    tuple val(meta), path("*.mosdepth.summary.txt"), emit: summary
    tuple val(meta), path("*.per-base.bed.gz"), emit: per_base
}
```

**New Parameters**:

```groovy
params.min_target_coverage = 500   // Minimum required coverage
params.coverage_threshold = 0.9    // % of targets meeting min coverage
```

---

#### 2.2 Add Insert Size Metrics

**Current Gap**: Insert size distribution not captured, which can indicate sample quality issues.

**Scientific Impact**: Abnormal insert sizes indicate:

- DNA degradation (shift to smaller fragments)
- Library prep issues (bimodal distribution)
- ctDNA enrichment (tumor DNA often shorter than germline)

**Recommended Tool**: Picard CollectInsertSizeMetrics

**New Process**:

```groovy
process PICARD_COLLECTINSERTSIZEMETRICS {
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*.insert_size_metrics.txt"), emit: metrics
    tuple val(meta), path("*.insert_size_histogram.pdf"), emit: histogram
}
```

---

#### 2.3 Improve Reference Genome Handling

**Current Issue**: The FASTA index (`.fai`) is created inside processes but not cached, causing redundant indexing.

**Current Code**:

```bash
if [ ! -f ${fasta}.fai ]; then
    samtools faidx ${fasta}
fi
```

This check fails because `fasta` is staged fresh each time.

**Recommended Solution**: Create dedicated indexing process:

```groovy
process SAMTOOLS_FAIDX {
    input:
    path(fasta)
    
    output:
    path("*.fai"), emit: fai
    
    script:
    """
    samtools faidx ${fasta}
    """
}

// In workflow:
ch_fasta_fai = params.genome_fasta_fai 
    ? Channel.fromPath(params.genome_fasta_fai) 
    : SAMTOOLS_FAIDX(ch_fasta).fai
```

---

#### 2.4 Add Strand Bias Filtering

**Current Gap**: LoFreq detects strand bias but it's not used for filtering.

**Scientific Impact**: Variants called predominantly from one strand are often artifacts from:

- PCR amplification errors
- Sequencing chemistry bias
- DNA damage (e.g., oxidative damage in FFPE)

**Recommended Enhancement**:

```bash
# In LOFREQ_CALL process, add strand bias threshold:
lofreq call-parallel \
    --sb-thresh 0.001 \   # Add this parameter
    ...

# Or post-filter:
bcftools filter -e 'SB > 100' input.vcf > filtered.vcf
```

**New Parameter**:

```groovy
params.max_strand_bias = 100  // Maximum strand bias score
```

---

#### 2.5 Add Variant Normalization

**Current Gap**: VCF files are not normalized, leading to inconsistent variant representations.

**Scientific Impact**: The same variant can be represented differently:

- `chr12:25398284 C>CA` vs `chr12:25398285 ->A`
- This causes issues when merging samples or comparing to databases

**Recommended Process**:

```groovy
process BCFTOOLS_NORM {
    input:
    tuple val(meta), path(vcf), path(tbi)
    path(fasta)
    
    output:
    tuple val(meta), path("*.norm.vcf.gz"), emit: vcf
    
    script:
    """
    bcftools norm \
        -f ${fasta} \
        -m -any \
        --check-ref w \
        -o ${meta.id}.norm.vcf.gz \
        -O z \
        ${vcf}
    """
}
```

---

### Priority 3: Nice-to-Have Features

#### 3.1 Expand Gene Panel Support

**Current Limitation**: Pipeline is hard-coded for KRAS region only.

**Recommended Enhancement**:

- Support multiple genes via BED file
- Add common ctDNA panel genes: TP53, EGFR, BRAF, PIK3CA, NRAS, APC
- Allow whole-exome or whole-genome mode

**New Parameters**:

```groovy
params.analysis_mode = 'targeted'  // 'targeted', 'exome', 'genome'
params.gene_panel = null           // BED file for targeted panel
params.exome_bed = null            // Exome capture regions
```

---

#### 3.2 Add Oncogene/Tumor Suppressor Annotation

**Current Gap**: SnpEff provides functional annotation but not cancer-specific interpretation.

**Recommended Databases**:

- **COSMIC**: Catalogue of Somatic Mutations in Cancer
- **ClinVar**: Clinical significance of variants
- **OncoKB**: Oncogenic effect and drug sensitivity
- **CGI (Cancer Genome Interpreter)**: Actionability

**New Process**:

```groovy
process ANNOTATE_CANCER_DB {
    input:
    tuple val(meta), path(vcf)
    path(cosmic_vcf)
    path(clinvar_vcf)
    
    output:
    tuple val(meta), path("*.cancer_annotated.vcf.gz"), emit: vcf
}
```

---

#### 3.3 Generate IGV Session Files

**Current Gap**: No easy way to visualize variants in context.

**Recommended Enhancement**: Generate IGV session XML for each sample:

```groovy
process CREATE_IGV_SESSION {
    input:
    tuple val(meta), path(bam), path(bai), path(vcf)
    path(fasta)
    
    output:
    path("*.igv_session.xml"), emit: session
}
```

---

#### 3.4 Create Clinical Summary Report

**Current Gap**: Results are in multiple files; no unified clinical interpretation.

**Recommended Enhancement**: Generate PDF/HTML report with:

- Sample QC summary (pass/fail)
- Coverage across target regions
- Detected variants with clinical interpretation
- VAF timeline (for longitudinal samples)
- Treatment recommendations (if integrated with OncoKB)

---

### Implementation Roadmap

| Phase | Enhancement | Estimated Effort | Impact |
|-------|-------------|------------------|--------|
| **Immediate** | Fix LoFreq container | 1 hour | Critical |
| **Week 1** | Add germline filtering | 4-8 hours | Critical |
| **Week 1** | Add BQSR | 4-8 hours | Critical |
| **Week 2** | Add mosdepth coverage | 2-4 hours | High |
| **Week 2** | Add variant normalization | 2-4 hours | High |
| **Week 3** | Add input validation | 4-8 hours | Medium |
| **Week 4** | Expand gene panel | 4-8 hours | Medium |
| **Future** | Cancer database annotation | 8-16 hours | Medium |
| **Future** | Clinical report generator | 16-24 hours | High |

---

## Current Limitations

1. **Not FDA-approved**: This pipeline is for research use only
2. **No germline filtering**: Cannot distinguish somatic vs germline variants
3. **Single-gene focus**: Currently optimized for KRAS only
4. **No tumor-normal pairing**: Requires matched PBMC for somatic calling
5. **No UMI support**: Cannot remove PCR duplicates at molecular level

---

## References

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html)
- [nf-core Guidelines](https://nf-co.re/docs/contributing/guidelines)
- [nf-core/sarek](https://nf-co.re/sarek) - Reference variant calling pipeline
- [LoFreq Publication](https://doi.org/10.1093/nar/gks918) - Low-frequency variant calling
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711) - Germline/somatic workflows
- [gnomAD](https://gnomad.broadinstitute.org/) - Population frequency database
