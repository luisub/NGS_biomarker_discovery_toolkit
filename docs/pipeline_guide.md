# VCA Pipeline Setup and Usage Guide

## Quick Start

### Make Script Executable

```bash
cd pipelines
chmod +x run_vca_pipeline.py
```

### Configure Pipeline

Edit `config/pipeline_config.yml` to customize:
- Patient ID filter
- Number of threads
- Target gene coordinates
- Quality filters

### Run Pipeline

```bash
# From pipelines/ directory
python run_vca_pipeline.py ../config/pipeline_config.yml

# Or if made executable:
./run_vca_pipeline.py ../config/pipeline_config.yml
```

## Directory Structure

After running the pipeline, your structure will be:

```
NGS_dPCR_biomarker_toolkit/
├── pipelines/
│   ├── run_vca_pipeline.py          # Main pipeline script
│   └── run_pipeline.sh              # Launcher script
├── config/
│   └── pipeline_config.yml          # Configuration file
├── data_cluster/                     # All data stored here
│   ├── raw/                          # FASTQ files from SRA
│   ├── reference/                    # GRCh38 reference genome
│   ├── aligned/                      # SAM/BAM alignment files
│   ├── variants/                     # VCF variant files
│   ├── metadata/                     # SRA metadata tables
│   └── results/                      # Final outputs
│       ├── all_variants.csv          # Complete variant table
│       ├── candidate_variant.csv     # Selected variant
│       ├── summary_stats.txt         # Summary statistics
│       └── vaf_over_time.png         # VAF plot
```

## Pipeline Steps

The pipeline executes the following steps:

1. **Reference Genome Setup**
   - Downloads GRCh38 reference
   - Creates BWA index

2. **Metadata Acquisition**
   - Fetches SRA metadata
   - Extracts patient IDs and timepoints

3. **Data Download**
   - Downloads SRA data (prefetch)
   - Converts to FASTQ (fasterq-dump)

4. **Read Alignment**
   - Aligns reads with BWA-MEM
   - Converts SAM to BAM
   - Sorts and indexes BAM

5. **Quality Processing**
   - Removes PCR duplicates
   - Indexes deduplicated BAM

6. **Variant Calling**
   - Calls variants with bcftools
   - Filters by depth and allele frequency

7. **Analysis & Visualization**
   - Identifies candidate variants
   - Generates VAF plots
   - Saves results

## Configuration Options

### Key Parameters

| Parameter | Description | Default | Notes |
|-----------|-------------|---------|-------|
| `patient_id_filter` | Specific patient to analyze | "CTC030" | Set to `null` for all |
| `threads` | Number of CPU cores | 4 | Adjust based on system |
| `min_depth` | Minimum read depth | 20 | Higher = more confident |
| `min_allele_frequency` | Minimum VAF | 0.1 | 0.1 = 10% |

### Target Gene Configuration

To analyze a different gene, modify `variant_calling.target_gene`:

```yaml
target_gene:
  name: "TP53"
  chromosome: "chr17"
  start: 7668421
  end: 7687490
```


## Output Files

### all_variants.csv

All detected variants meeting filter criteria:
- `run_id`: SRA run accession
- `patient_id`: Patient identifier
- `timepoint`: Treatment stage
- `chromosome`: Chromosome location
- `position`: Genomic position
- `ref`: Reference allele
- `alt`: Alternate allele
- `depth`: Read depth
- `allele_frequency`: Variant allele frequency


### Batch Processing Multiple Patients

Remove patient filter in config:
```yaml
patient_id_filter: null
```

### Custom Quality Filters

```yaml
filters:
  min_depth: 30           # More stringent
  min_allele_frequency: 0.05  # More sensitive
```

### Running on HPC Cluster

```bash
#!/bin/bash
#SBATCH --job-name=vca_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

conda activate vca_env
cd pipelines
python run_vca_pipeline.py ../config/pipeline_config.yml
```