# Troubleshooting Guide

This guide covers common issues and solutions when running the VCA Nextflow pipeline.

---

## Table of Contents

1. [Installation Issues](#installation-issues)
2. [Singularity Container Issues](#singularity-container-issues)
3. [Resource Errors](#resource-errors)
4. [Input/Output Errors](#inputoutput-errors)
5. [Tool-Specific Errors](#tool-specific-errors)
6. [HPC/Cluster Issues](#hpccluster-issues)

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
