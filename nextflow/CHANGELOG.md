# Changelog

All notable changes to the VCA Nextflow Pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.1.0] - 2026-01-09

### Added

- **Coverage Analysis**: `MOSDEPTH` process for depth of coverage statistics
- **Insert Size Metrics**: `PICARD_COLLECTINSERTSIZEMETRICS` for library QC
- **Variant Normalization**: `BCFTOOLS_NORM` for left-alignment and multi-allelic splitting
- **Germline Filtering**: `BCFTOOLS_FILTER_GERMLINE` with gnomAD integration
- **VCF Statistics**: `BCFTOOLS_STATS` for variant call quality metrics
- **FAIDX Process**: Dedicated `SAMTOOLS_FAIDX` for reference genome indexing

### Changed

- **LoFreq Container**: Now uses mulled container with samtools and htslib included
- **LoFreq Process**: Takes FAIDX output as explicit input (no redundant indexing)
- **MultiQC**: Now includes mosdepth, insert size, and bcftools stats

### New Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gnomad_vcf` | null | Path to gnomAD VCF for germline filtering |
| `--gnomad_tbi` | null | Path to gnomAD VCF index |
| `--max_population_af` | 0.01 | Maximum population AF to consider somatic |
| `--max_strand_bias` | null | Strand bias threshold for filtering |
| `--skip_coverage` | false | Skip mosdepth and Picard metrics |
| `--skip_germline_filter` | true | Skip germline filtering |
| `--skip_normalization` | false | Skip VCF normalization |
| `--whole_genome` | false | Analyze whole genome instead of target region |

---

## [1.0.0] - 2026-01-09

### Added

- Initial release of the Nextflow VCA Pipeline
- **Core Processes**:
  - `FASTQC` - Quality control of raw reads
  - `FASTP` - Read trimming and adapter removal
  - `BWA_INDEX` - Reference genome indexing
  - `BWA_MEM` - Read alignment
  - `SAMTOOLS_SORT` - BAM sorting
  - `SAMTOOLS_MARKDUP` - Duplicate marking
  - `SAMTOOLS_INDEX` - BAM indexing
  - `SAMTOOLS_FLAGSTAT` - Alignment statistics
  - `LOFREQ_CALL` - Low-frequency variant calling
  - `SNPEFF_ANNOTATE` - Variant annotation
  - `AGGREGATE_VARIANTS` - VCF merging and summary
  - `MULTIQC` - QC report aggregation

- **Subworkflows**:
  - `preprocessing.nf` - QC and trimming
  - `alignment.nf` - BWA and BAM processing
  - `variant_calling.nf` - LoFreq and annotation

- **Configuration**:
  - Multiple execution profiles (standard, docker, singularity, slurm, pbs, awsbatch)
  - Resource labels for process resource allocation
  - Comprehensive parameter validation

- **Documentation**:
  - Main README with quick start guide
  - Troubleshooting guide
  - Docker container documentation
  - Samplesheet format specification

### Execution Profiles

- `standard` - Local execution with Conda
- `docker` - Docker container execution
- `singularity` - Singularity for HPC (recommended)
- `slurm` - SLURM scheduler integration
- `pbs` - PBS/Torque scheduler integration
- `awsbatch` - AWS Batch cloud execution
- `test` - Minimal resources for testing
- `debug` - Verbose logging

---

## [Unreleased]

### Planned

- GATK BaseRecalibrator (BQSR) integration
- Panel of Normals (PoN) support
- UMI-based duplicate removal
- OncoKB/ClinVar annotation
- Clinical summary report generation
- IGV session file generation
