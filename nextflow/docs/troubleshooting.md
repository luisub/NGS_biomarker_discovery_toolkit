# Troubleshooting Guide

This guide covers common issues and solutions when running the VCA Nextflow pipeline.

---

## Table of Contents

1. [Installation Issues](#installation-issues)
2. [Singularity Container Issues](#singularity-container-issues)
3. [Resource Errors](#resource-errors)
4. [Input/Output Errors](#inputoutput-errors)
5. [Tool-Specific Errors](#tool-specific-errors)
6. [BQSR Issues](#bqsr-issues)
7. [Germline Filtering Issues](#germline-filtering-issues)
8. [HPC/Cluster Issues](#hpccluster-issues)

---

## Installation Issues

### Nextflow not found

**Error:**

```
command not found: nextflow
```

**Solution:**

```bash
# Install via conda
conda install -c bioconda nextflow

# Or install directly
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### Java version incompatible

**Error:**

```
ERROR: Cannot find Java or it's a wrong version -- please make sure that Java 11 (or later) is installed
```

**Solution:**

```bash
# Install Java 11+
conda install -c conda-forge openjdk=17

# Verify
java -version
```

---

## Singularity Container Issues

### Container pull timeout

**Error:**

```
WARN: Singularity image pull failed (exit 1) -- retrying
```

**Solution:**

```bash
# Increase pull timeout in nextflow.config
singularity {
    pullTimeout = '120 min'
}

# Or pre-pull the container manually
singularity pull docker://quay.io/biocontainers/bwa:0.7.18--he4a0461_0
```

### Bind mount errors

**Error:**

```
FATAL: container creation failed: mount /scratch: no such file or directory
```

**Solution:**
Edit `conf/singularity.config` to match your HPC file systems:

```groovy
singularity {
    runOptions = "--bind /your/path:/your/path"
}
```

### Permission denied in container

**Error:**

```
FATAL: could not open image: permission denied
```

**Solution:**

```bash
# Set proper cache directory permissions
export SINGULARITY_CACHEDIR=/path/to/writable/dir
chmod 755 $SINGULARITY_CACHEDIR

# Or use a different cache location
nextflow run main.nf -profile singularity \
    --singularity.cacheDir /tmp/singularity_cache
```

---

## Resource Errors

### Out of memory (OOM)

**Error:**

```
Process `BWA_MEM (sample1)` terminated with exit status 137
```

**Solution:**
Increase memory allocation in `nextflow.config`:

```groovy
process {
    withName: 'BWA_MEM' {
        memory = '32 GB'
    }
}
```

Or use command line:

```bash
nextflow run main.nf --max_memory '64.GB'
```

### Job killed by scheduler

**Error:**

```
slurmstepd: error: Exceeded job memory limit
```

**Solution:**

```bash
# Request more resources
nextflow run main.nf -profile slurm \
    --max_memory '128.GB' \
    --max_cpus 16
```

### Disk space full

**Error:**

```
No space left on device
```

**Solution:**

```bash
# Clean work directory
nextflow clean -f

# Use scratch space for work directory
nextflow run main.nf -w /scratch/$USER/nf_work
```

---

## Input/Output Errors

### Samplesheet parsing error

**Error:**

```
ERROR ~ Cannot find file: /path/to/sample.fastq.gz
```

**Solution:**

1. Check file paths in samplesheet are absolute paths
2. Verify files exist and are readable:

```bash
# Check each file listed in samplesheet
while IFS=, read -r sample fq1 fq2 patient timepoint; do
    [ -f "$fq1" ] || echo "Missing: $fq1"
    [ -f "$fq2" ] || echo "Missing: $fq2"
done < samplesheet.csv
```

### Invalid samplesheet format

**Error:**

```
ERROR ~ Samplesheet header does not contain required columns
```

**Solution:**
Ensure your samplesheet has the correct header:

```csv
sample,fastq_1,fastq_2,patient_id,timepoint
```

### Reference genome not indexed

**Error:**

```
[E::bwa_idx_load_from_disk] fail to locate the index files
```

**Solution:**

```bash
# Let the pipeline create the index
nextflow run main.nf --genome_fasta /path/to/reference.fa

# Or provide pre-built index
nextflow run main.nf --bwa_index /path/to/bwa_index/
```

---

## Tool-Specific Errors

### LoFreq: No variants called

**Symptoms:** Empty VCF file

**Possible causes:**

1. BAM file has no reads in target region
2. Quality thresholds too stringent

**Solution:**

```bash
# Check BAM coverage in target region
samtools depth -r chr12:25200000-25250000 sample.bam | head

# Lower thresholds
nextflow run main.nf \
    --min_read_depth 10 \
    --min_allele_frequency 0.005
```

### SnpEff: Database not found

**Error:**

```
ERROR: Database 'GRCh38.p14' not found
```

**Solution:**

```bash
# Download database first
snpEff download GRCh38.p14

# Or specify cache directory
nextflow run main.nf --snpeff_cache /path/to/snpeff/data
```

### FastQC: Memory error

**Error:**

```
Exception in thread "main" java.lang.OutOfMemoryError
```

**Solution:**

```groovy
// In nextflow.config
process {
    withName: 'FASTQC' {
        memory = '8 GB'
    }
}
```

---

## BQSR Issues

### LoFreq Viterbi: Reference genome mismatch

**Error:**

```text
[E::fai_retrieve] Region "chr12:25200000-25250000" not found in FASTA file
```

**Solution:**
Ensure your reference genome uses the same chromosome naming convention as your BAM files:

```bash
# Check chromosome names in FASTA
grep "^>" reference.fa | head

# Check chromosome names in BAM
samtools view -H sample.bam | grep "^@SQ" | head

# If mismatch, use correct reference or rename chromosomes
```

### BQSR: Memory issues with large BAMs

**Error:**

```text
Process `LOFREQ_VITERBI` terminated with exit status 137
```

**Solution:**

```bash
# Increase memory allocation
nextflow run main.nf --max_memory '64.GB'
```

Or in `nextflow.config`:

```groovy
process {
    withName: 'LOFREQ_VITERBI' {
        memory = '32 GB'
    }
}
```

### BQSR: Skipping when not desired

**Symptoms:** Pipeline skips BQSR step unexpectedly

**Solution:**
Check your parameters - BQSR runs by default unless disabled:

```bash
# Ensure BQSR is enabled (default)
nextflow run main.nf --skip_bqsr false

# To verify BQSR ran, check for .viterbi.bam files in output
ls results/aligned/*.viterbi.bam
```

### BQSR: No improvement in variant calls

**Symptoms:** Variant calls similar with/without BQSR

**Possible causes:**

1. High-quality input data (BQSR has less effect)
2. Very low coverage regions
3. Sequencing platform with accurate base qualities

**Solution:**
Compare base quality distributions before and after:

```bash
# Check BQSR stats file
cat results/aligned/sample.viterbi_stats.txt

# Compare quality scores
samtools view input.bam | head -1000 | cut -f11 | fold -w1 | sort | uniq -c
samtools view output.viterbi.bam | head -1000 | cut -f11 | fold -w1 | sort | uniq -c
```

---

## Germline Filtering Issues

### gnomAD file not found

**Error:**

```text
ERROR ~ Cannot find gnomAD VCF: /path/to/gnomad.vcf.gz
```

**Solution:**

```bash
# Download gnomAD for your region (chr12 for KRAS)
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr12.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr12.vcf.bgz.tbi

# Run with gnomAD files
nextflow run main.nf \
    --gnomad_vcf /path/to/gnomad.exomes.v4.0.sites.chr12.vcf.bgz \
    --gnomad_tbi /path/to/gnomad.exomes.v4.0.sites.chr12.vcf.bgz.tbi
```

### gnomAD index missing

**Error:**

```text
[E::hts_open_format] Failed to open file gnomad.vcf.gz.tbi
```

**Solution:**

```bash
# Create index if missing
tabix -p vcf gnomad.exomes.v4.0.sites.chr12.vcf.bgz

# Verify index exists
ls -la gnomad.exomes.v4.0.sites.chr12.vcf.bgz.tbi
```

### All variants filtered as germline

**Symptoms:** Empty or nearly empty somatic VCF

**Possible causes:**

1. AF threshold too strict
2. Tumor sample has significant germline contamination
3. Chromosome naming mismatch between VCF and gnomAD

**Solution:**

```bash
# Increase AF threshold (allow more variants through)
nextflow run main.nf --max_population_af 0.05

# Check chromosome naming
bcftools view -h your_variants.vcf.gz | grep "^##contig"
bcftools view -h gnomad.vcf.gz | grep "^##contig"

# Check filter statistics
cat results/variants/sample.filter_stats.txt
```

### No variants filtered (all kept)

**Symptoms:** Somatic VCF identical to input

**Possible causes:**

1. gnomAD file doesn't contain your variants' positions
2. Chromosome mismatch (chr12 vs 12)
3. gnomAD annotation failed silently

**Solution:**

```bash
# Check if gnomAD contains your region
tabix gnomad.vcf.gz chr12:25200000-25250000 | head

# Check annotation worked
bcftools query -f '%CHROM\t%POS\t%INFO/gnomAD_AF\n' annotated.vcf.gz | head

# Verify chromosome naming matches
bcftools view -H your.vcf.gz | cut -f1 | sort -u
bcftools view -H gnomad.vcf.gz | cut -f1 | sort -u
```

### Germline filter skipped unexpectedly

**Symptoms:** No `.somatic.vcf.gz` files created

**Solution:**
Germline filtering requires gnomAD files AND must be enabled:

```bash
# Enable germline filtering (disabled by default)
nextflow run main.nf \
    --skip_germline_filter false \
    --gnomad_vcf /path/to/gnomad.vcf.gz \
    --gnomad_tbi /path/to/gnomad.vcf.gz.tbi
```

### bcftools annotate: Region not found

**Error:**

```text
[W::bcf_sr_add_reader] No BGZF EOF marker; file may be truncated
```

**Solution:**
gnomAD file may be corrupted or incomplete:

```bash
# Verify file integrity
bgzip -t gnomad.exomes.v4.0.sites.chr12.vcf.bgz

# If corrupted, re-download
wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr12.vcf.bgz

# Recreate index
tabix -f -p vcf gnomad.exomes.v4.0.sites.chr12.vcf.bgz
```

---

## HPC/Cluster Issues

### SLURM job pending indefinitely

**Symptoms:** Jobs stay in `PD` state

**Solution:**

```bash
# Check queue status
squeue -u $USER

# Use different partition
nextflow run main.nf -profile slurm \
    --clusterOptions '--partition=short'
```

### Module load failures

**Error:**

```
module: command not found
```

**Solution:**
Add module initialization to process:

```groovy
process {
    beforeScript = 'source /etc/profile.d/modules.sh'
}
```

### Scratch directory not accessible

**Error:**

```
FATAL: Cannot access work directory
```

**Solution:**

```bash
# Use accessible scratch space
export NXF_WORK=/accessible/scratch/$USER/nf_work
nextflow run main.nf
```

---

## Debugging Tips

### Enable verbose logging

```bash
nextflow run main.nf -with-trace -with-report -with-timeline
```

### Check specific process logs

```bash
# Find work directory for failed process
cat .nextflow.log | grep "ERROR"

# Examine process directory
ls -la work/xx/xxxxxxxx/
cat work/xx/xxxxxxxx/.command.log
cat work/xx/xxxxxxxx/.command.err
```

### Resume failed pipeline

```bash
nextflow run main.nf -resume
```

### Run in debug mode

```bash
nextflow run main.nf -profile debug
```

---

## Getting Help

1. Check Nextflow documentation: <https://www.nextflow.io/docs/latest/>
2. Search nf-core Slack: <https://nf-co.re/join>
3. Open an issue on GitHub with:
   - Full error message
   - `nextflow log` output
   - System information (`nextflow info`)
