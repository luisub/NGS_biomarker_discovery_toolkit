# Nextflow Variant Calling Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/singularity-compatible-blue.svg)](https://sylabs.io/singularity/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Version 1.1.0** | Variant Calling Analysis for Circulating Tumor DNA

A production-ready Nextflow pipeline for detecting low-frequency variants in circulating tumor DNA (ctDNA) samples. Optimized for cancer biomarker discovery with variant allele frequencies as low as 1%.

---

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Pipeline Workflow](#pipeline-workflow)
- [Installation](#installation)
- [Usage](#usage)
- [Samplesheet Format](#samplesheet-format)
- [Output Files](#output-files)
- [Configuration](#configuration)
- [Documentation](#documentation)
- [Known Limitations](#known-limitations)
- [References](#references)

---

## Features

- **Low-Frequency Variant Detection**: Optimized for ctDNA with VAF thresholds as low as 0.1%
- **Reproducibility**: Fully containerized with Docker and Singularity support
- **Scalability**: Parallel processing across samples; runs on laptops, HPC clusters, or cloud
- **Germline Filtering**: Optional gnomAD-based filtering to identify somatic variants
- **Comprehensive QC**: Coverage analysis, insert size metrics, and aggregated MultiQC reports
- **Resume Capability**: Automatic checkpointing—restart failed runs without reprocessing

---

## Quick Start

```bash
# Basic run with Singularity (recommended for HPC)
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    -profile singularity

# Run with Docker
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    -profile docker

# Resume a failed run
nextflow run main.nf -resume
```

---

## Pipeline Workflow

```text
┌─────────────────────────────────────────────────────────────────────────┐
│                      VCA Pipeline Workflow                               │
└─────────────────────────────────────────────────────────────────────────┘

  FASTQ ──► FastQC ──► fastp ──► BWA-MEM ──► Sort ──► MarkDup ──► Index
              │          │                               │
              ▼          ▼                               ▼
         QC Report   Trimmed                      ┌──────────────┐
                     Reads                        │  MOSDEPTH    │ Coverage
                                                  │  PICARD      │ Insert Size
                                                  └──────┬───────┘
                                                         │
                                                         ▼
                                    ┌────────────────────────────────┐
                                    │           LOFREQ               │
                                    │  (Low-frequency variant calls) │
                                    └───────────────┬────────────────┘
                                                    │
                                                    ▼
                                          ┌─────────────────┐
                                          │  BCFTOOLS NORM  │
                                          │  (Normalize)    │
                                          └────────┬────────┘
                                                   │
                                                   ▼
                                          ┌─────────────────┐
                                          │     SNPEFF      │
                                          │  (Annotate)     │
                                          └────────┬────────┘
                                                   │
                                                   ▼
                                    ┌──────────────────────────┐
                                    │  GERMLINE FILTER         │
                                    │  (Optional, gnomAD)      │
                                    └────────────┬─────────────┘
                                                 │
                                                 ▼
                                       ┌──────────────────┐
                                       │    AGGREGATE     │
                                       │  (Merge samples) │
                                       └─────────┬────────┘
                                                 │
                                                 ▼
                                          ┌───────────┐
                                          │  MULTIQC  │
                                          └───────────┘
```

---

## Installation

### Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04.0
- [Singularity](https://sylabs.io/singularity/) or [Docker](https://www.docker.com/)
- Reference genome (GRCh38 recommended)

### Install Nextflow

```bash
# Using conda (recommended)
conda install -c bioconda nextflow

# Or direct installation
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### Download Reference Genome

```bash
# Download GRCh38 from Ensembl
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Index the reference (or let the pipeline do it)
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

---

## Usage

### Basic Usage

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/reference.fa \
    --outdir results \
    -profile singularity
```

### With Pre-built Indexes

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/reference.fa \
    --bwa_index /path/to/bwa_index/ \
    -profile singularity
```

### With Germline Filtering

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/reference.fa \
    --gnomad_vcf /path/to/gnomad.vcf.gz \
    --gnomad_tbi /path/to/gnomad.vcf.gz.tbi \
    --skip_germline_filter false \
    -profile singularity
```

### On SLURM Cluster

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/reference.fa \
    --max_memory 128.GB \
    --max_cpus 32 \
    -profile singularity,slurm
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Samplesheet CSV file |
| `--genome_fasta` | required | Reference genome FASTA |
| `--outdir` | `results` | Output directory |
| `--min_allele_frequency` | `0.01` | Minimum VAF (1%) |
| `--min_read_depth` | `20` | Minimum read depth |

For complete parameter reference, see [docs/PARAMETERS.md](docs/PARAMETERS.md).

---

## Samplesheet Format

Create a CSV file with your sample information:

```csv
sample,fastq_1,fastq_2,patient_id,timepoint
Sample_01,/data/Sample_01_R1.fastq.gz,/data/Sample_01_R2.fastq.gz,Patient_A,baseline
Sample_02,/data/Sample_02_R1.fastq.gz,/data/Sample_02_R2.fastq.gz,Patient_A,progression
Sample_03,/data/Sample_03_R1.fastq.gz,,Patient_B,baseline
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `fastq_1` | Yes | Path to R1 FASTQ (gzipped) |
| `fastq_2` | No | Path to R2 FASTQ (leave empty for single-end) |
| `patient_id` | No | Patient identifier for grouping |
| `timepoint` | No | Sample timepoint (baseline, progression, etc.) |

---

## Output Files

```text
results/
├── aligned/              # BAM files
│   ├── *.markdup.bam     # Deduplicated alignments
│   └── *.bai             # BAM indexes
├── variants/             # VCF files
│   ├── *.lofreq.vcf.gz   # Raw variant calls
│   ├── *.norm.vcf.gz     # Normalized variants
│   └── *.ann.vcf.gz      # Annotated variants
├── coverage/             # Coverage statistics
│   └── *.mosdepth.*      # Depth metrics
├── stats/                # QC metrics
│   ├── *.flagstat.txt    # Alignment statistics
│   └── *.insert_size*    # Insert size metrics
├── results/              # Final summaries
│   ├── all_variants.csv  # Combined variant table
│   └── summary_stats.csv # Per-sample summary
└── multiqc/              # Aggregated QC report
    └── multiqc_report.html
```

For detailed output descriptions, see [docs/OUTPUT.md](docs/OUTPUT.md).

---

## Configuration

### Execution Profiles

| Profile | Description |
|---------|-------------|
| `singularity` | Use Singularity containers (recommended for HPC) |
| `docker` | Use Docker containers |
| `slurm` | Submit jobs to SLURM scheduler |
| `pbs` | Submit jobs to PBS scheduler |
| `awsbatch` | Run on AWS Batch |
| `test` | Quick test with minimal resources |

Combine profiles as needed: `-profile singularity,slurm`

### Resource Requirements

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| BWA_MEM | 8 | 32 GB | 4h |
| LOFREQ_CALL | 4 | 16 GB | 2h |
| SNPEFF | 2 | 8 GB | 30m |
| Other | 1-2 | 4-8 GB | <30m |

Adjust limits with `--max_memory` and `--max_cpus`.

---

## Documentation

| Document | Description |
|----------|-------------|
| [PARAMETERS.md](docs/PARAMETERS.md) | Complete parameter reference |
| [OUTPUT.md](docs/OUTPUT.md) | Output file descriptions |
| [troubleshooting.md](docs/troubleshooting.md) | Common issues and solutions |
| [docker/README.md](docker/README.md) | Container build instructions |
| [CHANGELOG.md](CHANGELOG.md) | Version history |

---

## Known Limitations

| Limitation | Description |
|------------|-------------|
| **Research Use Only** | This pipeline is not FDA-approved for clinical diagnostics |
| **No BQSR** | Base quality score recalibration is not yet implemented |
| **Target Region** | Currently optimized for KRAS (chr12:25200000-25250000) |
| **No UMI Support** | Molecular barcodes for PCR duplicate removal not supported |

---

## Directory Structure

```text
nextflow/
├── main.nf                  # Main workflow
├── nextflow.config          # Configuration
├── CHANGELOG.md             # Version history
│
├── modules/local/           # Process definitions
│   ├── fastqc.nf            # Quality control
│   ├── fastp.nf             # Read trimming
│   ├── bwa.nf               # Alignment
│   ├── samtools.nf          # BAM processing
│   ├── lofreq.nf            # Variant calling
│   ├── snpeff.nf            # Annotation
│   ├── bcftools.nf          # VCF manipulation
│   ├── mosdepth.nf          # Coverage
│   ├── picard.nf            # Metrics
│   └── multiqc.nf           # Report aggregation
│
├── docs/                    # Documentation
│   ├── PARAMETERS.md
│   ├── OUTPUT.md
│   └── troubleshooting.md
│
├── docker/                  # Container files
│   └── Dockerfile
│
└── test/                    # Test data
    └── test_samplesheet.csv
```

---

## References

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html)
- [nf-core Guidelines](https://nf-co.re/docs/contributing/guidelines)
- [nf-core/sarek](https://nf-co.re/sarek) - Variant calling pipeline reference
- [LoFreq Publication](https://doi.org/10.1093/nar/gks918) - Wilm et al., 2012
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711)
- [gnomAD](https://gnomad.broadinstitute.org/) - Population frequency database

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## Support

For issues and questions:

1. Check [troubleshooting.md](docs/troubleshooting.md)
2. Search existing issues on GitHub
3. Open a new issue with:
   - Full error message
   - `nextflow log` output
   - System information (`nextflow info`)
