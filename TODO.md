# TODO - Variant Calling Analysis Project

## Top Priority

### Reproduce original paper analysis
- [ ] Task: Implement comprehensive quality control metrics for read alignment and variant calls.
    - [ ] Task: Add base quality score filtering (minimum Q30) and mapping quality thresholds.
    - [ ] Task: Integrate duplicate read removal using samtools markdup.
    - [ ] Task: Add read depth coverage analysis and visualization per sample.
    - [ ] Task: Implement strand bias detection to filter potential artifacts.

### Data Quality
- [ ] Task: Include normal tissue/germline samples to filter out non-ctDNA specific variants.
- [ ] Task: Add population frequency databases (gnomAD, 1000 Genomes) to filter common SNPs.
- [ ] Task: Implement variant allele frequency (VAF) threshold optimization.
- [ ] Task: Add dbSNP annotation to distinguish known polymorphisms from somatic mutations.
- [ ] Task: Create filtering pipeline to remove sequencing artifacts and technical noise.

### Code Organization
- [ ] Task: Migrate notebook code to modular Python source files.
- [ ] Task: Create command-line interface (CLI) for BASH execution.
- [ ] Task: Add configuration file support YAML for pipeline parameters.
- [ ] Task: Implement proper error handling and logging throughout pipeline.
- [ ] Task: Add Nextflow workflow for reproducibility.

### Scalability
- [ ] Task: Process data for all patients in the PRJNA714799 dataset.
    - [ ] Task: Create batch processing script to handle multiple patients simultaneously.
    - [ ] Task: Implement parallel processing for faster alignment and variant calling.
    - [ ] Task: Create summary statistics table for all processed patients.

### Gene Panel Expansion
- [ ] Task: Expand analysis beyond KRAS to include full colorectal cancer gene panel.
    - [ ] Task: Add TP53, APC, PIK3CA, BRAF, and NRAS variant detection.
    - [ ] Task: Create automated gene coordinate lookup for any gene of interest.
    - [ ] Task: Implement variant effect prediction (VEP, SnpEff) for functional annotation.
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
- [ ] Task: Create comprehensive documentation for pipeline usage.
    - [ ] Task: Write methods section suitable for publication.
    - [ ] Task: Generate example analysis reports for clinical interpretation.
    - [ ] Task: Create visualization dashboard for variant tracking over time.
    - [ ] Task: Document computational requirements and runtime benchmarks.

## Future Directions

### Advanced Analytics
- [ ] Task: Integrate machine learning for treatment response prediction.
- [ ] Task: Look deeper in the dataset to explore why some patients have fewer time points.
    - [ ] Task: Identify variants that may be more aggressive and resulted in early patient death. 


