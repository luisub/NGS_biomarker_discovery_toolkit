#!/usr/bin/env python3
"""
Variant Calling Analysis Pipeline for Circulating Tumor DNA
Implements the pipeline described in vca_pipeline.ipynb
by: Luis Aguilera, December 2, 2025.
"""

import subprocess
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import re
import yaml
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysradb import SRAweb
import requests
import os
import logging
from typing import Dict, List, Tuple, Optional
from datetime import datetime
try:
    from plots_sequences import plot_gene_and_variants, plot_protein_mutations, get_protein_features
except ImportError:
    # Handle case where script is run from different directory
    sys.path.append(str(Path(__file__).parent))
    from plots_sequences import plot_gene_and_variants, plot_protein_mutations, get_protein_features

class VCAConfig:
    """Configuration manager for VCA pipeline."""
    def __init__(self, config_path: Path):
        """Load configuration from YAML file."""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        self._setup_paths()
    
    def _setup_paths(self):
        """Create all required directory paths."""
        base_dir = Path(self.config['paths']['base_dir'])
        self.data_dir = base_dir / 'data_cluster'
        self.raw_dir = self.data_dir / 'raw'
        self.reference_dir = self.data_dir / 'reference'
        self.aligned_dir = self.data_dir / 'aligned'
        self.variants_dir = self.data_dir / 'variants'
        self.metadata_dir = self.data_dir / 'metadata'
        self.results_dir = self.data_dir / 'results'
        self.qc_dir = self.data_dir / 'qc'
        for directory in [self.raw_dir, self.reference_dir, self.aligned_dir,
                         self.variants_dir, self.metadata_dir, self.results_dir, self.qc_dir]:
            directory.mkdir(parents=True, exist_ok=True)

class VCAPipeline:
    """Main pipeline for variant calling analysis.
    
    This class orchestrates the complete variant calling workflow including:
    - Reference genome download and indexing
    - Sample metadata fetching
    - Read QC, trimming, and alignment
    - Duplicate removal and BQSR
    - Coverage analysis
    - Variant calling, annotation, and filtering
    - Results aggregation and visualization
    
    Attributes:
        config: VCAConfig object with pipeline settings
        metadata_df: DataFrame with SRA metadata
        sample_info: Dictionary mapping run IDs to sample information
        logger: Logger instance for this pipeline run
    """
    
    def __init__(self, config: VCAConfig, log_file: Optional[Path] = None):
        """Initialize pipeline with configuration.
        
        Args:
            config: VCAConfig object with pipeline settings
            log_file: Optional path to log file. If None, creates timestamped log
                     in the QC directory.
        """
        self.config = config
        self.metadata_df = None
        self.sample_info = {}
        self._setup_logging(log_file)
    
    def _setup_logging(self, log_file: Optional[Path] = None) -> None:
        """Configure logging to both console and file.
        
        Args:
            log_file: Optional path to log file
        """
        # Create logger
        self.logger = logging.getLogger(f'VCAPipeline_{id(self)}')
        self.logger.setLevel(logging.DEBUG)
        
        # Clear any existing handlers
        self.logger.handlers = []
        
        # Console handler (INFO level)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_format = logging.Formatter('[%(levelname)s] %(message)s')
        console_handler.setFormatter(console_format)
        self.logger.addHandler(console_handler)
        
        # File handler (DEBUG level - more verbose)
        if log_file is None:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            log_file = self.config.qc_dir / f'pipeline_{timestamp}.log'
        
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_format = logging.Formatter(
            '%(asctime)s [%(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_format)
        self.logger.addHandler(file_handler)
        
        self.logger.info(f"Pipeline log: {log_file}")
    
    def run_command(self, cmd: List[str], step_name: str) -> bool:
        """Execute shell command with error handling and logging.
        
        Args:
            cmd: Command and arguments as list
            step_name: Human-readable name for logging
            
        Returns:
            True if command succeeded, False otherwise
        """
        self.logger.debug(f"Running: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info(f"{step_name}")
            print(f"[OK] {step_name}")
            return True
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr[:200] if e.stderr else str(e)
            self.logger.error(f"{step_name} failed: {error_msg}")
            print(f"[ERROR] {step_name} failed: {error_msg[:100]}")
            return False
        except FileNotFoundError:
            self.logger.error(f"{step_name} failed: Command not found - {cmd[0]}")
            print(f"[ERROR] {step_name} failed: Command not found")
            return False

    
    def download_reference_genome(self) -> bool:
        """Download and index GRCh38 reference genome."""
        ref_config = self.config.config['reference_genome']
        ref_url = ref_config['url']
        ref_path = self.config.reference_dir / ref_config['filename']
        if ref_path.exists():
            print(f"[SKIP] Reference genome already exists")
            return True
        print(f"Downloading reference genome...")
        try:
            response = requests.get(ref_url, stream=True)
            response.raise_for_status()
            with open(ref_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print("[OK] Reference genome downloaded")
        except Exception as e:
            print(f"[ERROR] Reference download failed: {str(e)[:100]}")
            return False
        index_files = [ref_path.with_suffix(ref_path.suffix + ext) 
                      for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        if all(f.exists() for f in index_files):
            print("[SKIP] Reference genome already indexed")
            return True
        cmd = ['bwa', 'index', str(ref_path)]
        return self.run_command(cmd, "BWA index creation")

    def download_dbsnp(self) -> Path:
        """Download common dbSNP VCF for GRCh38."""
        dbsnp_vcf = self.config.reference_dir / "common_all_20180418.vcf.gz"
        dbsnp_tbi = self.config.reference_dir / "common_all_20180418.vcf.gz.tbi"
        
        if dbsnp_vcf.exists() and dbsnp_tbi.exists():
            print(f"[SKIP] dbSNP files already exist")
            return dbsnp_vcf
            
        print(f"Downloading dbSNP common variants...")
        dbsnp_url = "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz"
        dbsnp_tbi_url = "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi"
        
        try:
            # Download VCF
            print(f"Downloading VCF from {dbsnp_url}...")
            with requests.get(dbsnp_url, stream=True) as r:
                r.raise_for_status()
                with open(dbsnp_vcf, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            
            # Download Index
            print(f"Downloading Index from {dbsnp_tbi_url}...")
            with requests.get(dbsnp_tbi_url, stream=True) as r:
                r.raise_for_status()
                with open(dbsnp_tbi, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                        
            print("[OK] dbSNP download completed")
            return dbsnp_vcf
        except Exception as e:
            print(f"[WARNING] dbSNP download failed: {e}")
            if dbsnp_vcf.exists(): dbsnp_vcf.unlink()
            if dbsnp_tbi.exists(): dbsnp_tbi.unlink()
            return None

    def download_gnomad(self) -> Path:
        """Download gnomAD VCF for germline filtering.
        
        Downloads population allele frequency data from gnomAD to filter
        common germline variants from somatic variant calls.
        
        Returns:
            Path to gnomAD VCF file, or None if download fails or disabled.
        """
        germline_config = self.config.config.get('germline_filtering', {})
        
        if not germline_config.get('enabled', False):
            print("[SKIP] Germline filtering is disabled")
            return None
        
        gnomad_filename = germline_config.get('gnomad_filename', 'gnomad.exomes.v4.0.sites.chr12.vcf.bgz')
        gnomad_url = germline_config.get('gnomad_url')
        
        gnomad_vcf = self.config.reference_dir / gnomad_filename
        gnomad_tbi = self.config.reference_dir / (gnomad_filename + '.tbi')
        
        if gnomad_vcf.exists() and gnomad_tbi.exists():
            print(f"[SKIP] gnomAD files already exist: {gnomad_filename}")
            return gnomad_vcf
        
        if not gnomad_url:
            print("[WARNING] gnomAD URL not configured, skipping download")
            return None
            
        print(f"Downloading gnomAD population frequencies...")
        print(f"  URL: {gnomad_url}")
        print(f"  This may take a while (~2GB for chr12)...")
        
        try:
            # Download VCF
            with requests.get(gnomad_url, stream=True) as r:
                r.raise_for_status()
                total_size = int(r.headers.get('content-length', 0))
                downloaded = 0
                with open(gnomad_vcf, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192*16):
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            pct = (downloaded / total_size) * 100
                            print(f"\r  Progress: {pct:.1f}%", end='', flush=True)
                print()  # New line after progress
            
            # Download Index
            tbi_url = gnomad_url + '.tbi'
            print(f"Downloading index from {tbi_url}...")
            with requests.get(tbi_url, stream=True) as r:
                r.raise_for_status()
                with open(gnomad_tbi, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                        
            print("[OK] gnomAD download completed")
            return gnomad_vcf
        except Exception as e:
            print(f"[WARNING] gnomAD download failed: {e}")
            if gnomad_vcf.exists(): gnomad_vcf.unlink()
            if gnomad_tbi.exists(): gnomad_tbi.unlink()
            return None

    
    def fetch_sra_metadata(self) -> bool:
        """Download SRA metadata for bioproject."""
        bioproject = self.config.config['data_source']['bioproject_id']
        metadata_path = self.config.metadata_dir / f"{bioproject}_metadata.csv"
        if metadata_path.exists():
            print(f"[SKIP] Metadata already exists")
            self.metadata_df = pd.read_csv(metadata_path)
            return True
        try:
            db = SRAweb()
            metadata = db.sra_metadata(bioproject, detailed=True)
            metadata.to_csv(metadata_path, index=False)
            self.metadata_df = metadata
            print(f"[OK] Metadata downloaded: {len(metadata)} samples")
            return True
        except Exception as e:
            print(f"[ERROR] Metadata download failed: {str(e)[:100]}")
            return False
    
    def extract_sample_info(self) -> bool:
        """Extract patient IDs and timepoints from metadata."""
        if self.metadata_df is None:
            print("[ERROR] No metadata available")
            return False
        try:
            pattern = r'patient:\s*(\w+).*?timepoint:\s*(\w+)'
            for _, row in self.metadata_df.iterrows():
                run_id = row['run_accession']
                description = str(row.get('sample_attribute', ''))
                match = re.search(pattern, description, re.IGNORECASE)
                if match:
                    patient_id = match.group(1)
                    timepoint = match.group(2)
                    self.sample_info[run_id] = {'patient_id': patient_id, 'timepoint': timepoint}
            print(f"[OK] Extracted info for {len(self.sample_info)} samples")
            return True
        except Exception as e:
            print(f"[ERROR] Sample info extraction failed: {str(e)[:100]}")
            return False
    
    def download_sra_data(self, run_id: str) -> bool:
        """Download SRA data using prefetch and fasterq-dump."""
        fastq_path = self.config.raw_dir / f"{run_id}_1.fastq"
        if fastq_path.exists():
            print(f"[SKIP] FASTQ already exists: {run_id}")
            return True
        prefetch_cmd = ['prefetch', run_id, '-O', str(self.config.raw_dir)]
        if not self.run_command(prefetch_cmd, f"Prefetch {run_id}"):
            return False
        dump_cmd = ['fasterq-dump', run_id, '-O', str(self.config.raw_dir), 
                   '-e', str(self.config.config['processing']['threads'])]
        return self.run_command(dump_cmd, f"FASTQ dump {run_id}")

    def run_fastqc(self, run_id: str) -> bool:
        """Run FastQC on raw FASTQ files."""
        fastq_r1 = self.config.raw_dir / f"{run_id}_1.fastq"
        fastq_r2 = self.config.raw_dir / f"{run_id}_2.fastq"
        
        if not fastq_r1.exists():
             print(f"[ERROR] FASTQ file not found: {fastq_r1}")
             return False
        
        cmd = ['fastqc', '-t', str(self.config.config['processing']['threads']), '-o', str(self.config.qc_dir), str(fastq_r1)]
        if fastq_r2.exists():
            cmd.append(str(fastq_r2))
            
        return self.run_command(cmd, f"FastQC {run_id}")

    def trim_reads(self, run_id: str) -> bool:
        """Trim reads using fastp."""
        fastq_r1 = self.config.raw_dir / f"{run_id}_1.fastq"
        fastq_r2 = self.config.raw_dir / f"{run_id}_2.fastq"
        
        trimmed_r1 = self.config.raw_dir / f"{run_id}_1.trimmed.fastq"
        trimmed_r2 = self.config.raw_dir / f"{run_id}_2.trimmed.fastq"
        
        html_report = self.config.qc_dir / f"{run_id}_fastp.html"
        json_report = self.config.qc_dir / f"{run_id}_fastp.json"
        
        if trimmed_r1.exists():
            print(f"[SKIP] Trimmed reads already exist: {run_id}")
            return True
            
        cmd = [
            "fastp",
            "-i", str(fastq_r1), "-I", str(fastq_r2),
            "-o", str(trimmed_r1), "-O", str(trimmed_r2),
            "-h", str(html_report), "-j", str(json_report),
            "--detect_adapter_for_pe",
            "--thread", str(self.config.config['processing']['threads'])
        ]
        
        return self.run_command(cmd, f"fastp trimming {run_id}")
    
    def align_reads(self, run_id: str) -> bool:
        """Align reads using BWA-MEM."""
        ref_path = self.config.reference_dir / self.config.config['reference_genome']['filename']
        # Use trimmed reads
        fastq_r1 = self.config.raw_dir / f"{run_id}_1.trimmed.fastq"
        fastq_r2 = self.config.raw_dir / f"{run_id}_2.trimmed.fastq"
        sam_path = self.config.aligned_dir / f"{run_id}.sam"
        
        if not fastq_r1.exists():
            print(f"[ERROR] Trimmed FASTQ file not found: {fastq_r1}")
            return False
        if sam_path.exists():
            print(f"[SKIP] Alignment already exists: {run_id}")
            return True
        
        threads = str(self.config.config['processing']['threads'])
        cmd = ['bwa', 'mem', '-t', threads, str(ref_path), str(fastq_r1), str(fastq_r2)]
        
        try:
            with open(sam_path, 'w') as out_file:
                subprocess.run(cmd, stdout=out_file, check=True, stderr=subprocess.PIPE)
            print(f"[OK] Alignment completed: {run_id}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Alignment failed: {e.stderr.decode()[:100]}")
            return False
    
    def convert_sort_index_bam(self, run_id: str) -> bool:
        """Convert SAM to BAM, sort, and index."""
        sam_path = self.config.aligned_dir / f"{run_id}.sam"
        bam_path = self.config.aligned_dir / f"{run_id}_sorted.bam"
        if bam_path.exists() and Path(str(bam_path) + '.bai').exists():
            print(f"[SKIP] Sorted BAM already exists: {run_id}")
            return True
        view_cmd = ['samtools', 'view', '-bS', str(sam_path)]
        sort_cmd = ['samtools', 'sort', '-o', str(bam_path), '-']
        try:
            view_proc = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            sort_proc = subprocess.Popen(sort_cmd, stdin=view_proc.stdout, stderr=subprocess.PIPE)
            view_proc.stdout.close()
            sort_proc.communicate()
            if sort_proc.returncode != 0:
                print(f"[ERROR] BAM sorting failed")
                return False
        except Exception as e:
            print(f"[ERROR] BAM conversion failed: {str(e)[:100]}")
            return False
        index_cmd = ['samtools', 'index', str(bam_path)]
        return self.run_command(index_cmd, f"BAM indexing {run_id}")
    
    def remove_duplicates(self, run_id: str) -> bool:
        """Remove PCR duplicates using samtools markdup."""
        bam_path = self.config.aligned_dir / f"{run_id}_sorted.bam"
        dedup_path = self.config.aligned_dir / f"{run_id}_dedup.bam"
        if dedup_path.exists():
            print(f"[SKIP] Deduplicated BAM already exists: {run_id}")
            return True
        markdup_cmd = ['samtools', 'markdup', '-r', str(bam_path), str(dedup_path)]
        if not self.run_command(markdup_cmd, f"Duplicate removal {run_id}"):
            return False
        index_cmd = ['samtools', 'index', str(dedup_path)]
        return self.run_command(index_cmd, f"Dedup BAM indexing {run_id}")

    def calculate_coverage_stats(self, run_id: str) -> bool:
        """Calculate depth of coverage statistics for the target region.
        
        Computes coverage metrics including mean depth, median depth,
        and fraction of bases covered at various thresholds (10x, 50x, 100x, 500x).
        This is critical for ctDNA analysis where high coverage is essential.
        
        Args:
            run_id: Sample run identifier
            
        Returns:
            True if coverage analysis succeeded, False otherwise
            
        Note:
            Uses samtools depth for coverage calculation. For more detailed
            analysis, consider installing mosdepth (faster, more features).
        """
        # Use BQSR BAM if available, otherwise use deduplicated BAM
        bqsr_bam = self.config.aligned_dir / f"{run_id}_dedup.bqsr.bam"
        dedup_bam = self.config.aligned_dir / f"{run_id}_dedup.bam"
        bam_path = bqsr_bam if bqsr_bam.exists() else dedup_bam
        
        coverage_file = self.config.qc_dir / f"{run_id}_coverage_stats.txt"
        
        if coverage_file.exists():
            print(f"[SKIP] Coverage stats already exist: {run_id}")
            return True
        
        if not bam_path.exists():
            print(f"[ERROR] BAM not found for coverage analysis: {bam_path}")
            return False
        
        gene_config = self.config.config['variant_calling']['target_gene']
        region = f"{gene_config['chromosome']}:{gene_config['start']}-{gene_config['end']}"
        
        try:
            print(f"  Calculating coverage for {region}...")
            
            # Get per-base depth using samtools depth
            depth_cmd = [
                "samtools", "depth",
                "-r", region,
                "-a",  # Output all positions including zero coverage
                str(bam_path)
            ]
            result = subprocess.run(depth_cmd, capture_output=True, text=True, check=True)
            
            # Parse depth values
            depths = []
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        depths.append(int(parts[2]))
            
            if not depths:
                print(f"[WARNING] No coverage data for region {region}")
                depths = [0]
            
            # Calculate statistics
            depths_array = np.array(depths)
            total_bases = len(depths_array)
            mean_depth = np.mean(depths_array)
            median_depth = np.median(depths_array)
            min_depth = np.min(depths_array)
            max_depth = np.max(depths_array)
            
            # Coverage thresholds (important for ctDNA)
            thresholds = [10, 50, 100, 500, 1000]
            coverage_at = {}
            for t in thresholds:
                coverage_at[t] = (depths_array >= t).sum() / total_bases * 100
            
            # Write stats file
            with open(coverage_file, 'w') as f:
                f.write(f"=== Coverage Statistics for {run_id} ===\n")
                f.write(f"Region: {region}\n")
                f.write(f"BAM: {bam_path.name}\n\n")
                f.write(f"Total bases in region: {total_bases:,}\n")
                f.write(f"Mean depth: {mean_depth:.1f}x\n")
                f.write(f"Median depth: {median_depth:.1f}x\n")
                f.write(f"Min depth: {min_depth}x\n")
                f.write(f"Max depth: {max_depth}x\n\n")
                f.write("Coverage thresholds:\n")
                for t in thresholds:
                    f.write(f"  >= {t}x: {coverage_at[t]:.1f}%\n")
                
                # Flag if coverage is insufficient for ctDNA
                min_ctdna_coverage = self.config.config.get('variant_calling', {}).get(
                    'filters', {}).get('min_depth', 100)
                if median_depth < min_ctdna_coverage:
                    f.write(f"\n⚠️  WARNING: Median coverage ({median_depth:.0f}x) below ")
                    f.write(f"recommended minimum ({min_ctdna_coverage}x) for ctDNA analysis\n")
            
            print(f"[OK] Coverage stats: mean={mean_depth:.0f}x, median={median_depth:.0f}x")
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Coverage analysis failed: {e.stderr[:200] if e.stderr else str(e)}")
            return False
        except Exception as e:
            print(f"[ERROR] Coverage analysis failed: {str(e)[:100]}")
            return False


    def recalibrate_base_qualities(self, run_id: str) -> bool:
        """Recalibrate base qualities using LoFreq Viterbi.
        
        Uses Hidden Markov Model to recalibrate base quality scores,
        improving accuracy of low-frequency variant detection.
        
        Args:
            run_id: Sample run identifier
            
        Returns:
            True if BQSR succeeded, False otherwise
        """
        bqsr_config = self.config.config.get('bqsr', {})
        
        if not bqsr_config.get('enabled', True):
            print(f"[SKIP] BQSR disabled for {run_id}")
            return True
        
        ref_path = self.config.reference_dir / self.config.config['reference_genome']['filename']
        input_bam = self.config.aligned_dir / f"{run_id}_dedup.bam"
        output_bam = self.config.aligned_dir / f"{run_id}_dedup.bqsr.bam"
        stats_file = self.config.qc_dir / f"{run_id}_bqsr_stats.txt"
        
        if output_bam.exists():
            print(f"[SKIP] BQSR BAM already exists: {run_id}")
            return True
        
        if not input_bam.exists():
            print(f"[ERROR] Deduplicated BAM not found: {input_bam}")
            return False
        
        try:
            # Run LoFreq Viterbi for base quality recalibration
            print(f"  Running LoFreq Viterbi BQSR...")
            viterbi_cmd = [
                "lofreq", "viterbi",
                "-f", str(ref_path),
                "-o", str(output_bam),
                str(input_bam)
            ]
            result = subprocess.run(viterbi_cmd, check=True, capture_output=True, text=True)
            
            # Index recalibrated BAM
            index_cmd = ["samtools", "index", str(output_bam)]
            subprocess.run(index_cmd, check=True, capture_output=True)
            
            # Generate statistics
            with open(stats_file, 'w') as f:
                f.write(f"=== LoFreq Viterbi BQSR Statistics ===\n")
                f.write(f"Sample: {run_id}\n")
                f.write(f"Input BAM: {input_bam}\n")
                f.write(f"Output BAM: {output_bam}\n\n")
                
                # Get flagstat for recalibrated BAM
                flagstat_result = subprocess.run(
                    ["samtools", "flagstat", str(output_bam)],
                    capture_output=True, text=True
                )
                f.write("Alignment statistics:\n")
                f.write(flagstat_result.stdout)
            
            print(f"[OK] BQSR completed: {run_id}")
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] BQSR failed: {e.stderr[:200] if e.stderr else str(e)}")
            return False
        except Exception as e:
            print(f"[ERROR] BQSR failed: {str(e)[:100]}")
            return False
    
    def call_variants_lofreq(self, run_id: str, dbsnp_path: Path = None) -> bool:
        """Call variants using Lofreq.
        
        Uses BQSR-recalibrated BAM if available, otherwise falls back to
        deduplicated BAM.
        """
        ref_path = self.config.reference_dir / self.config.config['reference_genome']['filename']
        
        # Prefer BQSR BAM if available
        bqsr_bam = self.config.aligned_dir / f"{run_id}_dedup.bqsr.bam"
        dedup_bam = self.config.aligned_dir / f"{run_id}_dedup.bam"
        
        if bqsr_bam.exists():
            bam_path = bqsr_bam
            print(f"  Using BQSR-recalibrated BAM for variant calling")
        else:
            bam_path = dedup_bam
            print(f"  Using deduplicated BAM (no BQSR)")
        
        vcf_path = self.config.variants_dir / f"{run_id}.lofreq.vcf"

        
        if vcf_path.exists():
            print(f"[SKIP] VCF already exists: {run_id}")
            return True
            
        # Indel qualities
        bam_indel_path = bam_path.with_suffix(".indel.bam")
        if not bam_indel_path.exists():
            cmd_indel = ["lofreq", "indelqual", "--dindel", "-f", str(ref_path), "-o", str(bam_indel_path), str(bam_path)]
            if not self.run_command(cmd_indel, f"Lofreq indelqual {run_id}"):
                return False
            self.run_command(["samtools", "index", str(bam_indel_path)], f"Index indel BAM {run_id}")
            
        gene_config = self.config.config['variant_calling']['target_gene']
        region = f"{gene_config['chromosome']}:{gene_config['start']}-{gene_config['end']}"
        
        cmd_call = [
            "lofreq", "call",
            "-f", str(ref_path),
            "-r", region,
            "-o", str(vcf_path),
            "--call-indels",
            str(bam_indel_path)
        ]
        
        if dbsnp_path:
            cmd_call.extend(["-d", str(dbsnp_path)])

        
        return self.run_command(cmd_call, f"Lofreq call {run_id}")

    def annotate_variants(self, run_id: str) -> bool:
        """Annotate variants using SnpEff."""
        input_vcf = self.config.variants_dir / f"{run_id}.lofreq.vcf"
        output_vcf = self.config.variants_dir / f"{run_id}.lofreq.ann.vcf"
        output_vcf_gz = self.config.variants_dir / f"{run_id}.lofreq.ann.vcf.gz"
        
        if output_vcf_gz.exists():
            print(f"[SKIP] Annotated VCF already exists: {run_id}")
            return True
            
        # Set Java heap
        os.environ["_JAVA_OPTIONS"] = "-Xmx4g"
        
        try:
            with open(output_vcf, "w") as f:
                subprocess.run(["snpEff", "-v", "GRCh38.86", str(input_vcf)], stdout=f, check=True, stderr=subprocess.PIPE)
            
            pysam.tabix_index(str(output_vcf), preset="vcf", force=True)
            if output_vcf.exists():
                output_vcf.unlink() # Remove uncompressed
            print(f"[OK] Annotation completed: {run_id}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Annotation failed: {e.stderr.decode()[:100]}")
            return False
        except Exception as e:
            print(f"[ERROR] Annotation failed: {str(e)[:100]}")
            return False

    def filter_germline_variants(self, run_id: str, gnomad_path: Path) -> bool:
        """Filter germline variants using gnomAD population frequencies.
        
        Uses bcftools to annotate variants with gnomAD allele frequencies
        and filter out common germline variants (AF > threshold).
        
        Args:
            run_id: Sample run identifier
            gnomad_path: Path to gnomAD VCF file
            
        Returns:
            True if filtering succeeded, False otherwise
        """
        germline_config = self.config.config.get('germline_filtering', {})
        
        if not germline_config.get('enabled', False):
            print(f"[SKIP] Germline filtering disabled for {run_id}")
            return True
            
        if gnomad_path is None or not gnomad_path.exists():
            print(f"[WARNING] gnomAD file not available, skipping germline filter for {run_id}")
            return True
        
        input_vcf = self.config.variants_dir / f"{run_id}.lofreq.ann.vcf.gz"
        annotated_vcf = self.config.variants_dir / f"{run_id}.annotated.vcf.gz"
        output_vcf = self.config.variants_dir / f"{run_id}.lofreq.ann.somatic.vcf.gz"
        stats_file = self.config.variants_dir / f"{run_id}.germline_filter_stats.txt"
        
        if output_vcf.exists():
            print(f"[SKIP] Somatic VCF already exists: {run_id}")
            return True
        
        if not input_vcf.exists():
            print(f"[ERROR] Annotated VCF not found: {input_vcf}")
            return False
        
        max_af = germline_config.get('max_population_af', 0.01)
        
        try:
            # Step 1: Annotate with gnomAD allele frequencies
            print(f"  Annotating with gnomAD AF...")
            annotate_cmd = [
                "bcftools", "annotate",
                "-a", str(gnomad_path),
                "-c", "INFO/gnomAD_AF:=INFO/AF",
                "-O", "z",
                "-o", str(annotated_vcf),
                str(input_vcf)
            ]
            result = subprocess.run(annotate_cmd, check=True, capture_output=True, text=True)
            
            # Index annotated VCF
            subprocess.run(["tabix", "-p", "vcf", str(annotated_vcf)], check=True, capture_output=True)
            
            # Step 2: Filter out common germline variants
            print(f"  Filtering germline variants (AF > {max_af})...")
            filter_cmd = [
                "bcftools", "filter",
                "-e", f"INFO/gnomAD_AF > {max_af}",
                "-s", "GERMLINE",
                "-m", "+",
                "-O", "z",
                "-o", str(output_vcf),
                str(annotated_vcf)
            ]
            result = subprocess.run(filter_cmd, check=True, capture_output=True, text=True)
            
            # Index filtered VCF
            subprocess.run(["tabix", "-p", "vcf", str(output_vcf)], check=True, capture_output=True)
            
            # Generate statistics
            before_count = subprocess.run(
                ["bcftools", "view", "-H", str(input_vcf)],
                capture_output=True, text=True
            ).stdout.count('\n')
            
            pass_count = subprocess.run(
                ["bcftools", "view", "-f", "PASS", "-H", str(output_vcf)],
                capture_output=True, text=True
            ).stdout.count('\n')
            
            germline_count = subprocess.run(
                ["bcftools", "view", "-f", "GERMLINE", "-H", str(output_vcf)],
                capture_output=True, text=True
            ).stdout.count('\n')
            
            with open(stats_file, 'w') as f:
                f.write(f"=== Germline Filter Statistics ===\n")
                f.write(f"Sample: {run_id}\n")
                f.write(f"Max population AF threshold: {max_af}\n\n")
                f.write(f"Variants before filtering: {before_count}\n")
                f.write(f"Variants after filtering (PASS): {pass_count}\n")
                f.write(f"Variants filtered as GERMLINE: {germline_count}\n")
            
            # Cleanup intermediate file
            if annotated_vcf.exists():
                annotated_vcf.unlink()
                Path(str(annotated_vcf) + '.tbi').unlink(missing_ok=True)
            
            print(f"[OK] Germline filtering completed: {run_id}")
            print(f"     Filtered {germline_count}/{before_count} variants as germline")
            print(f"     Remaining somatic candidates: {pass_count}")
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Germline filtering failed: {e.stderr[:200] if e.stderr else str(e)}")
            return False
        except Exception as e:
            print(f"[ERROR] Germline filtering failed: {str(e)[:100]}")
            return False
    
    def parse_variants(self) -> pd.DataFrame:
        """Parse VCF files and extract variant information.
        
        Prefers somatic (germline-filtered) VCFs if available, otherwise
        falls back to annotated VCFs. Skips variants marked as GERMLINE.
        """
        variants_list = []
        filters = self.config.config['variant_calling']['filters']
        min_depth = filters['min_depth']
        min_freq = filters['min_allele_frequency']
        
        for run_id, info in self.sample_info.items():
            # Prefer somatic VCF (germline-filtered) if available
            somatic_vcf_path = self.config.variants_dir / f"{run_id}.lofreq.ann.somatic.vcf.gz"
            annotated_vcf_path = self.config.variants_dir / f"{run_id}.lofreq.ann.vcf.gz"
            
            if somatic_vcf_path.exists():
                vcf_path = somatic_vcf_path
                print(f"  Using somatic VCF for {run_id}")
            elif annotated_vcf_path.exists():
                vcf_path = annotated_vcf_path
                print(f"  Using annotated VCF for {run_id} (no germline filtering)")
            else:
                print(f"  [WARNING] No VCF found for {run_id}")
                continue
                
            try:
                vcf = pysam.VariantFile(str(vcf_path))
                for record in vcf:
                    # Skip GERMLINE-filtered variants
                    if 'GERMLINE' in record.filter.keys():
                        continue
                        
                    depth = record.info.get('DP', 0)
                    if depth < min_depth:
                        continue
                    ref_allele = record.ref
                    alt_alleles = record.alts
                    if not alt_alleles:
                        continue
                    
                    # Lofreq AF is in AF info field
                    af = record.info.get('AF', [0.0])[0]
                    
                    # Get gnomAD AF if available (for reporting)
                    gnomad_af = record.info.get('gnomAD_AF', None)
                    
                    if af >= min_freq:
                        variant_entry = {
                            'run_id': run_id,
                            'patient_id': info['patient_id'],
                            'timepoint': info['timepoint'],
                            'chromosome': record.chrom,
                            'position': record.pos,
                            'ref': ref_allele,
                            'alt': alt_alleles[0],
                            'depth': depth,
                            'allele_frequency': af
                        }
                        if gnomad_af is not None:
                            variant_entry['gnomad_af'] = gnomad_af
                        variants_list.append(variant_entry)
                vcf.close()
            except Exception as e:
                print(f"[WARNING] Failed to parse VCF for {run_id}: {str(e)[:100]}")
                continue
        if not variants_list:
            print("[WARNING] No variants passed filters")
            return pd.DataFrame()
        variants_df = pd.DataFrame(variants_list)
        print(f"[OK] Parsed {len(variants_df)} variants from {len(set(variants_df['run_id']))} samples")
        return variants_df

    
    def analyze_variants(self, variants_df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict]:
        """Identify candidate variants and generate summary statistics."""
        if variants_df.empty:
            print("[ERROR] No variants to analyze")
            return pd.DataFrame(), {}
        position_counts = variants_df.groupby(['chromosome', 'position']).size()
        max_timepoints = len(variants_df['timepoint'].unique())
        persistent_positions = position_counts[position_counts == max_timepoints]
        if persistent_positions.empty:
            print("[WARNING] No persistent variants across all timepoints")
            candidate_position = variants_df.groupby(['chromosome', 'position']).size().idxmax()
        else:
            candidate_position = persistent_positions.idxmax()
        candidate_chrom, candidate_pos = candidate_position
        candidate_variants = variants_df[
            (variants_df['chromosome'] == candidate_chrom) & 
            (variants_df['position'] == candidate_pos)
        ].copy()
        summary_stats = {
            'total_variants': len(variants_df),
            'unique_positions': len(variants_df.groupby(['chromosome', 'position'])),
            'candidate_chromosome': candidate_chrom,
            'candidate_position': candidate_pos,
            'candidate_ref': candidate_variants.iloc[0]['ref'],
            'candidate_alt': candidate_variants.iloc[0]['alt'],
            'mean_depth': candidate_variants['depth'].mean(),
            'mean_allele_freq': candidate_variants['allele_frequency'].mean()
        }
        print(f"[OK] Identified candidate at {candidate_chrom}:{candidate_pos}")
        return candidate_variants, summary_stats
    
    def create_visualizations(self, candidate_variants: pd.DataFrame):
        """Generate analysis plots."""
        if candidate_variants.empty:
            print("[WARNING] No candidate variants to visualize")
            return
        timepoint_order = ['pre_treatment', 'during_treatment', 'post_treatment']
        # Map numerical timepoints to labels if needed, or rely on existing logic
        # For now, assuming timepoint column matches
        
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(candidate_variants['timepoint'], candidate_variants['allele_frequency'],
               marker='o', linewidth=2, markersize=10, color='#2E86AB')
        ax.set_title('Variant Allele Frequency Over Treatment', fontsize=14, fontweight='bold')
        ax.set_xlabel('Treatment Stage', fontsize=12)
        ax.set_ylabel('Allele Frequency', fontsize=12)
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plot_path = self.config.results_dir / 'vaf_over_time.png'
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[OK] Visualization saved: vaf_over_time.png")

        # Generate Gene Plot
        try:
            gene_config = self.config.config['variant_calling']['target_gene']
            gene_name = gene_config.get('name', 'Target_Gene')
            region = f"{gene_config['chromosome']}:{gene_config['start']}-{gene_config['end']}"
            
            vcf_files = {}
            for run_id, info in self.sample_info.items():
                vcf_path = self.config.variants_dir / f"{run_id}.lofreq.ann.vcf.gz"
                if vcf_path.exists():
                    label = f"{info['patient_id']} ({info['timepoint']})"
                    vcf_files[label] = str(vcf_path)
            
            if vcf_files:
                # 1. Gene Structure Plot
                output_plot = self.config.results_dir / 'gene_variants_plot.png'
                plot_gene_and_variants(
                    gene_name=gene_name,
                    vcf_files=vcf_files,
                    genomic_region=region,
                    output_file=str(output_plot)
                )
                print(f"[OK] Visualization saved: gene_variants_plot.png")

                # 2. Protein Mutation Plot
                output_protein_plot = self.config.results_dir / 'protein_mutations.png'
                features, length = get_protein_features(gene_name)
                plot_protein_mutations(
                    gene_name=gene_name,
                    vcf_files=vcf_files,
                    protein_features=features,
                    protein_length=length,
                    output_file=str(output_protein_plot)
                )
                print(f"[OK] Visualization saved: protein_mutations.png")

        except Exception as e:
            print(f"[WARNING] Failed to generate plots: {e}")
    
    def save_results(self, variants_df: pd.DataFrame, candidate_variants: pd.DataFrame,
                    summary_stats: Dict):
        """Save analysis results to CSV files."""
        all_variants_path = self.config.results_dir / 'all_variants.csv'
        candidate_path = self.config.results_dir / 'candidate_variant.csv'
        summary_path = self.config.results_dir / 'summary_stats.txt'
        variants_df.to_csv(all_variants_path, index=False)
        candidate_variants.to_csv(candidate_path, index=False)
        with open(summary_path, 'w') as f:
            f.write("VCA Pipeline Summary\n")
            f.write("=" * 40 + "\n\n")
            for key, value in summary_stats.items():
                f.write(f"{key}: {value}\n")
        print(f"[OK] Results saved")
    
    def run_full_pipeline(self):
        """Execute complete VCA pipeline."""
        print("\nVCA PIPELINE - Variant Calling Analysis")
        print("=" * 60)
        if not self.download_reference_genome():
            return False
        if not self.fetch_sra_metadata():
            return False
        if not self.extract_sample_info():
            return False
            
        # Download dbSNP once
        dbsnp_path = self.download_dbsnp()
        
        # Download gnomAD for germline filtering
        gnomad_path = self.download_gnomad()
        
        patient_filter = self.config.config['data_source'].get('patient_id_filter', None)

        if patient_filter:
            self.sample_info = {k: v for k, v in self.sample_info.items() 
                               if v['patient_id'] == patient_filter}
            print(f"[INFO] Filtered to patient: {patient_filter}")
        for run_id in self.sample_info.keys():
            print(f"\nProcessing: {run_id}")
            if not self.download_sra_data(run_id):
                continue
            if not self.run_fastqc(run_id):
                continue
            if not self.trim_reads(run_id):
                continue
            if not self.align_reads(run_id):
                continue
            if not self.convert_sort_index_bam(run_id):
                continue
            if not self.remove_duplicates(run_id):
                continue
            # Base Quality Score Recalibration (if enabled)
            if not self.recalibrate_base_qualities(run_id):
                continue
            # Coverage analysis (QC)
            if not self.calculate_coverage_stats(run_id):
                print(f"[WARNING] Coverage analysis failed, continuing anyway")
            if not self.call_variants_lofreq(run_id, dbsnp_path):
                continue
            if not self.annotate_variants(run_id):
                continue
            # Germline filtering (if enabled and gnomAD available)
            if not self.filter_germline_variants(run_id, gnomad_path):
                continue
                
        print("\nVARIANT ANALYSIS")
        print("=" * 60)
        variants_df = self.parse_variants()
        if variants_df.empty:
            return False
        candidate_variants, summary_stats = self.analyze_variants(variants_df)
        self.create_visualizations(candidate_variants)
        self.save_results(variants_df, candidate_variants, summary_stats)
        print("\nPIPELINE COMPLETED")
        print("=" * 60 + "\n")
        return True


def main():
    """Main entry point for VCA pipeline."""
    if len(sys.argv) < 2:
        print("Usage: python run_vca_pipeline.py <config.yml>")
        sys.exit(1)
    config_path = Path(sys.argv[1])
    if not config_path.exists():
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)
    try:
        config = VCAConfig(config_path)
        pipeline = VCAPipeline(config)
        success = pipeline.run_full_pipeline()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"[FATAL] Pipeline failed: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
