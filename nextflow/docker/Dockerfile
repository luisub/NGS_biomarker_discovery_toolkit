# =============================================================================
# VCA Pipeline - Complete Docker Container
# Variant Calling Analysis for Circulating Tumor DNA
#
# Author: Luis Aguilera
# 
# Contains all bioinformatics tools needed for the pipeline:
#   - FastQC (quality control)
#   - fastp (read trimming)
#   - BWA (alignment)
#   - Samtools (BAM processing)
#   - LoFreq (variant calling)
#   - SnpEff (variant annotation)
#   - MultiQC (report aggregation)
#   - bcftools, htslib, tabix (VCF processing)
#
# Build:
#   docker build -t vca-pipeline:1.0.0 .
#   docker build -t yourusername/vca-pipeline:1.0.0 .
#
# Push to Docker Hub:
#   docker login
#   docker push yourusername/vca-pipeline:1.0.0
#
# Convert to Singularity:
#   singularity pull docker://yourusername/vca-pipeline:1.0.0
#
# Test the container:
#   docker run --rm vca-pipeline:1.0.0 bwa
#   docker run --rm vca-pipeline:1.0.0 samtools --version
#   docker run --rm vca-pipeline:1.0.0 lofreq version
# =============================================================================

FROM mambaorg/micromamba:1.5.6-jammy

LABEL maintainer="Luis Aguilera"
LABEL description="VCA Pipeline - Variant Calling Analysis for Circulating Tumor DNA"
LABEL version="1.0.0"
LABEL org.opencontainers.image.source="https://github.com/luisub/NGS_biomarker_discovery_toolkit"

# Set environment variables
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# Switch to root for system package installation
USER root

# Install essential system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    wget \
    curl \
    gzip \
    bzip2 \
    unzip \
    ca-certificates \
    fontconfig \
    libfreetype6 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Switch back to micromamba user
USER $MAMBA_USER

# Create conda environment with all bioinformatics tools
# Using specific versions for reproducibility
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

# Alternative: Install tools directly (if not using environment.yml)
# RUN micromamba install -y -n base -c conda-forge -c bioconda \
#     python=3.12 \
#     bwa=0.7.18 \
#     samtools=1.19 \
#     bcftools=1.19 \
#     htslib=1.19 \
#     lofreq=2.1.5 \
#     snpeff=5.2 \
#     fastqc=0.12.1 \
#     fastp=0.23.4 \
#     multiqc=1.21 \
#     && micromamba clean --all --yes

# Set working directory
WORKDIR /data

# Verify installations
RUN bwa 2>&1 | head -3 || true && \
    samtools --version | head -1 && \
    bcftools --version | head -1 && \
    lofreq version && \
    fastqc --version && \
    fastp --version 2>&1 | head -1 && \
    snpEff -version 2>&1 | head -1 && \
    multiqc --version

# Default entrypoint
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
