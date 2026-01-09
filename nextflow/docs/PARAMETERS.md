# Parameter Reference

Complete reference for all VCA Pipeline parameters.

---

## Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--input` | Path | Path to samplesheet CSV file (see [Samplesheet Format](#samplesheet-format)) |
| `--genome_fasta` | Path | Path to reference genome FASTA file |

---

## Reference Files

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genome_fasta_fai` | null | Pre-built FASTA index (`.fai`). If not provided, will be generated |
| `--bwa_index` | null | Path to pre-built BWA index directory. If not provided, will be generated |
| `--dbsnp` | null | Path to dbSNP VCF for annotation |
| `--dbsnp_tbi` | null | Path to dbSNP VCF index |
| `--snpeff_db` | `GRCh38.p14` | SnpEff database name |
| `--snpeff_cache` | null | Path to SnpEff cache directory |

---

## Target Region

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--target_bed` | null | BED file with target regions (overrides gene-based targeting) |
| `--target_gene` | `KRAS` | Target gene name (for logging purposes) |
| `--target_chromosome` | `chr12` | Chromosome of target region |
| `--target_start` | `25200000` | Start position of target region |
| `--target_end` | `25250000` | End position of target region |
| `--whole_genome` | `false` | If true, analyze whole genome (ignore target region) |

---

## Quality Filters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_base_quality` | `30` | Minimum base quality score (Phred) |
| `--min_mapping_quality` | `20` | Minimum mapping quality score |
| `--min_read_depth` | `20` | Minimum read depth at variant position |
| `--min_allele_frequency` | `0.01` | Minimum variant allele frequency (1% for ctDNA) |
| `--max_strand_bias` | null | Maximum strand bias score (null = no filter) |

---

## Germline Filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gnomad_vcf` | null | Path to gnomAD VCF for germline filtering |
| `--gnomad_tbi` | null | Path to gnomAD VCF index |
| `--panel_of_normals` | null | Path to Panel of Normals VCF |
| `--max_population_af` | `0.01` | Maximum population allele frequency to consider somatic |

---

## Read Trimming (fastp)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--adapter_fasta` | null | Custom adapter sequences FASTA |
| `--trim_front` | `0` | Trim N bases from front of reads |
| `--trim_tail` | `0` | Trim N bases from end of reads |
| `--min_read_length` | `50` | Minimum read length after trimming |

---

## Coverage Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_target_coverage` | `500` | Minimum required coverage for ctDNA analysis |

---

## Base Quality Score Recalibration (BQSR)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_bqsr` | `false` | Skip base quality score recalibration (LoFreq Viterbi) |
| `--save_bqsr_bam` | `false` | Save BQSR-recalibrated BAM files to output directory |

BQSR uses LoFreq's Viterbi algorithm to recalibrate base quality scores using a Hidden Markov Model. This improves accuracy of low-frequency variant detection by correcting systematic errors in base quality scores from the sequencer.

---

## Skip Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_fastqc` | `false` | Skip FastQC quality control |
| `--skip_trimming` | `false` | Skip fastp read trimming |
| `--skip_markdup` | `false` | Skip duplicate marking |
| `--skip_annotation` | `false` | Skip SnpEff variant annotation |
| `--skip_multiqc` | `false` | Skip MultiQC report generation |
| `--skip_coverage` | `false` | Skip mosdepth coverage analysis |
| `--skip_germline_filter` | `true` | Skip germline filtering (requires `--gnomad_vcf`) |
| `--skip_normalization` | `false` | Skip VCF normalization |

---

## Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_time` | `240.h` | Maximum time per process |
| `--threads` | `4` | Default thread count |

---

## Output Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Output directory path |
| `--publish_dir_mode` | `copy` | How to save outputs: `copy`, `link`, `symlink`, `move` |
| `--tracedir` | `${outdir}/pipeline_info` | Directory for execution reports |

---

## Utility Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--help` | `false` | Show help message and exit |
| `--version` | `false` | Show version and exit |
| `--monochrome_logs` | `false` | Disable colored log output |

---

## Samplesheet Format

The input samplesheet must be a CSV file with the following columns:

```csv
sample,fastq_1,fastq_2,patient_id,timepoint
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `fastq_1` | Yes | Path to R1 FASTQ file (gzipped) |
| `fastq_2` | No | Path to R2 FASTQ file (for paired-end) |
| `patient_id` | No | Patient identifier (defaults to sample) |
| `timepoint` | No | Sample timepoint (e.g., `baseline`, `follow_up`) |

### Example

```csv
sample,fastq_1,fastq_2,patient_id,timepoint
SRR13948001,/data/SRR13948001_1.fastq.gz,/data/SRR13948001_2.fastq.gz,CTC030,baseline
SRR13948002,/data/SRR13948002_1.fastq.gz,/data/SRR13948002_2.fastq.gz,CTC030,follow_up_1
SRR13948003,/data/SRR13948003_1.fastq.gz,,CTC031,baseline
```

Note: The third sample has no `fastq_2`, indicating single-end sequencing.

---

## Usage Examples

### Basic Run

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    -profile singularity
```

### With Pre-built Indexes

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    --genome_fasta_fai /path/to/GRCh38.fa.fai \
    --bwa_index /path/to/bwa_index/ \
    -profile singularity
```

### With Germline Filtering

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    --gnomad_vcf /path/to/gnomad.vcf.gz \
    --gnomad_tbi /path/to/gnomad.vcf.gz.tbi \
    --skip_germline_filter false \
    -profile singularity
```

### Low-stringency for Low-coverage Samples

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    --min_read_depth 10 \
    --min_allele_frequency 0.005 \
    -profile singularity
```

### SLURM Cluster with High Resources

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    --max_memory 256.GB \
    --max_cpus 32 \
    -profile singularity,slurm
```

### Skip Specific Steps

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    --skip_fastqc \
    --skip_coverage \
    -profile singularity
```

### With BQSR (Default - Recommended for ctDNA)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    --skip_bqsr false \
    --save_bqsr_bam true \
    -profile singularity
```

### Skip BQSR (Faster, Less Accurate)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    --skip_bqsr true \
    -profile singularity
```

### Full ctDNA Analysis (BQSR + Germline Filter)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta /path/to/GRCh38.fa \
    --gnomad_vcf /path/to/gnomad.exomes.v4.0.sites.chr12.vcf.bgz \
    --gnomad_tbi /path/to/gnomad.exomes.v4.0.sites.chr12.vcf.bgz.tbi \
    --skip_bqsr false \
    --skip_germline_filter false \
    --max_population_af 0.01 \
    --min_allele_frequency 0.01 \
    -profile singularity
```
