# TODO - Variant Calling Analysis Project

## Top Priority

### Reproduce original paper analysis
- [X] Task: Implement comprehensive quality control metrics for read alignment and variant calls.
    - [X] Task: Add base quality score filtering (minimum Q30) and mapping quality thresholds.
    - [X] Task: Integrate duplicate read removal using samtools markdup.
    - [ ] Task: Add read depth coverage analysis and visualization per sample.
    - [ ] Task: Implement strand bias detection to filter potential artifacts.

### Data Quality
- [ ] Task: Include normal tissue/germline samples to filter out non-ctDNA specific variants.
- [ ] Task: Add population frequency databases (gnomAD, 1000 Genomes) to filter common SNPs.
- [ ] Task: Implement variant allele frequency (VAF) threshold optimization.
- [ ] Task: Add dbSNP annotation to distinguish known polymorphisms from somatic mutations.
- [ ] Task: Create filtering pipeline to remove sequencing artifacts and technical noise.

### Code Organization
- [X] Task: Migrate notebook code to modular Python source files.
- [X] Task: Create command-line interface (CLI) for BASH execution.
- [X] Task: Add configuration file support YAML for pipeline parameters.
- [X] Task: Implement proper error handling and logging throughout pipeline.
- [X] Task: Restructure repository layout (pipelines, config, src).
- [ ] Task: Add Nextflow workflow for reproducibility.

### Scalability
- [ ] Task: Process data for all patients in the PRJNA714799 dataset.
    - [ ] Task: Create batch processing script to handle multiple patients simultaneously.
    - [ ] Task: Implement parallel processing for faster alignment and variant calling.
    - [ ] Task: Create summary statistics table for all processed patients.

## Digital PCR Assay Design
- [X] Task: Design allele-specific TaqMan probes (FAM for mutant, VIC for wildtype).
- [X] Task: Simulate droplet partitioning using Poisson statistics (20,000 droplets).
- [ ] Task: Create 2D fluorescence scatter plot showing four quadrants (WT+/MUT+, WT+/MUT-, WT-/MUT+, negative).
- [X] Task: Calculate limit of detection (LOD) for rare variants (target: 0.01% VAF sensitivity).
- [ ] Task: Simulate serial ctDNA monitoring during treatment (complete response, partial response, progression).

### Gene Panel Expansion
- [ ] Task: Expand analysis beyond KRAS to include full colorectal cancer gene panel.
    - [ ] Task: Add TP53, APC, PIK3CA, BRAF, and NRAS variant detection.
    - [X] Task: Create automated gene coordinate lookup for any gene of interest.
    - [X] Task: Implement variant effect prediction (VEP, SnpEff) for functional annotation.
    - [ ] Task: Generate comprehensive mutation timeline tables for all target genes.

## Second Priority

### Clinical Translation - PCR Diagnostic Development
- [X] Task: Design PCR primers targeting the most common KRAS mutation hotspots.
- [ ] Task: Develop multiplex PCR strategy to detect multiple variants simultaneously.
- [ ] Task: Optimize PCR conditions for low-input ctDNA samples.
- [ ] Task: Compare cost-effectiveness of PCR vs NGS for clinical monitoring.
- [X] Task: Create standard operating procedure (SOP) for PCR-based diagnostic.

## Long term goals

### Cross-Dataset Analysis
- [ ] Task: Identify and download similar circulating tumor DNA datasets from SRA/GEO.
- [ ] Task: Perform cross-dataset validation of detected KRAS variants.
- [ ] Task: Compare variant detection sensitivity across different sequencing platforms.
- [ ] Task: Validate detected variants against published results from original studies.
- [ ] Task: Implement statistical methods to assess variant call confidence.

### Regulatory and Clinical Implementation
- [ ] Task: Research FDA 510(k) requirements for ctDNA diagnostics.
    - [ ] Task: Prepare analytical validation plan with appropriate controls.
    - [ ] Task: Design clinical validation study protocol with patient cohorts.
    - [ ] Task: Document quality management system (QMS) requirements.
    - [ ] Task: Prepare pre-submission meeting materials for FDA consultation.
    - [ ] Task: Identify CLIA-certified laboratory partners for test development.
    - [ ] Task: Investigate reimbursement pathways and CPT coding.

### Documentation
- [X] Task: Create comprehensive documentation for pipeline usage.
    - [ ] Task: Write methods section suitable for publication.
    - [ ] Task: Generate example analysis reports for clinical interpretation.
    - [ ] Task: Create visualization dashboard for variant tracking over time.
    - [ ] Task: Document computational requirements and runtime benchmarks.

## Future Directions

### Advanced Analytics
- [ ] Task: Integrate machine learning for treatment response prediction.
- [ ] Task: Look deeper in the dataset to explore why some patients have fewer time points.
    - [ ] Task: Identify variants that may be more aggressive and resulted in early patient death. 


