# VCA Pipeline Docker Container

This directory contains the Dockerfile and environment specification for building a complete container with all bioinformatics tools needed for the VCA pipeline.

## Included Tools

| Tool | Version | Purpose |
|------|---------|---------|
| **FastQC** | 0.12.1 | Quality control of raw reads |
| **fastp** | 0.23.4 | Read trimming and filtering |
| **BWA** | 0.7.18 | Read alignment to reference genome |
| **Samtools** | 1.19 | BAM file manipulation |
| **bcftools** | 1.19 | VCF file manipulation |
| **LoFreq** | 2.1.5 | Low-frequency variant calling |
| **SnpEff** | 5.2 | Variant annotation |
| **MultiQC** | 1.21 | Aggregate QC reports |

## Building the Container

### Local Build

```bash
cd nextflow/docker

# Build with default tag
docker build -t vca-pipeline:1.0.0 .

# Build with your Docker Hub username
docker build -t yourusername/vca-pipeline:1.0.0 .
```

### Verify the Build

```bash
# Test each tool
docker run --rm vca-pipeline:1.0.0 bwa 2>&1 | head -3
docker run --rm vca-pipeline:1.0.0 samtools --version
docker run --rm vca-pipeline:1.0.0 lofreq version
docker run --rm vca-pipeline:1.0.0 fastqc --version
docker run --rm vca-pipeline:1.0.0 snpEff -version
docker run --rm vca-pipeline:1.0.0 multiqc --version
```

## Pushing to Docker Hub

```bash
# Login to Docker Hub
docker login

# Push the image
docker push yourusername/vca-pipeline:1.0.0

# Also tag as latest
docker tag yourusername/vca-pipeline:1.0.0 yourusername/vca-pipeline:latest
docker push yourusername/vca-pipeline:latest
```

## Converting to Singularity

On your HPC cluster:

```bash
# Pull and convert Docker image to Singularity
singularity pull docker://yourusername/vca-pipeline:1.0.0

# This creates: vca-pipeline_1.0.0.sif

# Test the Singularity image
singularity exec vca-pipeline_1.0.0.sif bwa
singularity exec vca-pipeline_1.0.0.sif samtools --version
```

## Using with Nextflow

### Option 1: Use the Custom Container for All Processes

Edit `nextflow.config`:

```groovy
process {
    container = 'yourusername/vca-pipeline:1.0.0'
}
```

### Option 2: Use with Singularity Profile

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome_fasta reference.fa \
    -profile singularity \
    -with-singularity yourusername/vca-pipeline:1.0.0
```

### Option 3: Pre-download for Air-Gapped Clusters

```bash
# On a machine with internet access
singularity pull docker://yourusername/vca-pipeline:1.0.0

# Transfer the .sif file to your cluster
scp vca-pipeline_1.0.0.sif user@cluster:/path/to/containers/

# Use the local .sif file
nextflow run main.nf \
    -profile singularity \
    -with-singularity /path/to/containers/vca-pipeline_1.0.0.sif
```

## Container Size

Expected size: ~2-3 GB (includes all tools and dependencies)

## Rebuilding After Updates

If you need to update tool versions:

1. Edit `environment.yml` with new versions
2. Rebuild the container:

   ```bash
   docker build --no-cache -t vca-pipeline:1.0.1 .
   ```

3. Push the new version
4. Update `nextflow.config` with the new tag
