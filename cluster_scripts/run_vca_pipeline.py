#!/usr/bin/env python3
"""
Variant Calling Analysis Pipeline for Circulating Tumor DNA
by: Luis Aguilera, November 21, 2025.
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
from typing import Dict, List, Tuple

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
        for directory in [self.raw_dir, self.reference_dir, self.aligned_dir,
                         self.variants_dir, self.metadata_dir, self.results_dir]:
            directory.mkdir(parents=True, exist_ok=True)

class VCAPipeline:
    """Main pipeline for variant calling analysis."""
    
    def __init__(self, config: VCAConfig):
        """Initialize pipeline with configuration."""
        self.config = config
        self.metadata_df = None
        self.sample_info = {}
    
    def run_command(self, cmd: List[str], step_name: str) -> bool:
        """Execute shell command with error handling."""
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"[OK] {step_name}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] {step_name} failed: {e.stderr[:100]}")
            return False
        except FileNotFoundError:
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
    
    def align_reads(self, run_id: str) -> bool:
        """Align reads using BWA-MEM."""
        ref_path = self.config.reference_dir / self.config.config['reference_genome']['filename']
        fastq_r1 = self.config.raw_dir / f"{run_id}_1.fastq"
        fastq_r2 = self.config.raw_dir / f"{run_id}_2.fastq"
        sam_path = self.config.aligned_dir / f"{run_id}.sam"
        if not fastq_r1.exists():
            print(f"[ERROR] FASTQ file not found: {fastq_r1}")
            return False
        if sam_path.exists():
            print(f"[SKIP] Alignment already exists: {run_id}")
            return True
        threads = str(self.config.config['processing']['threads'])
        if fastq_r2.exists():
            cmd = ['bwa', 'mem', '-t', threads, str(ref_path), str(fastq_r1), str(fastq_r2)]
        else:
            cmd = ['bwa', 'mem', '-t', threads, str(ref_path), str(fastq_r1)]
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
    
    def call_variants(self, run_id: str) -> bool:
        """Call variants using bcftools mpileup and call."""
        ref_path = self.config.reference_dir / self.config.config['reference_genome']['filename']
        bam_path = self.config.aligned_dir / f"{run_id}_dedup.bam"
        vcf_path = self.config.variants_dir / f"{run_id}.vcf"
        if vcf_path.exists():
            print(f"[SKIP] VCF already exists: {run_id}")
            return True
        gene_config = self.config.config['variant_calling']['target_gene']
        region = f"{gene_config['chromosome']}:{gene_config['start']}-{gene_config['end']}"
        mpileup_cmd = ['bcftools', 'mpileup', '-f', str(ref_path), '-r', region, str(bam_path)]
        call_cmd = ['bcftools', 'call', '-mv', '-Ov', '-o', str(vcf_path)]
        try:
            mpileup_proc = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            call_proc = subprocess.Popen(call_cmd, stdin=mpileup_proc.stdout, stderr=subprocess.PIPE)
            mpileup_proc.stdout.close()
            call_proc.communicate()
            if call_proc.returncode != 0:
                print(f"[ERROR] Variant calling failed")
                return False
            print(f"[OK] Variant calling completed: {run_id}")
            return True
        except Exception as e:
            print(f"[ERROR] Variant calling failed: {str(e)[:100]}")
            return False
    
    def parse_variants(self) -> pd.DataFrame:
        """Parse VCF files and extract variant information."""
        variants_list = []
        filters = self.config.config['variant_calling']['filters']
        min_depth = filters['min_depth']
        min_freq = filters['min_allele_frequency']
        for run_id, info in self.sample_info.items():
            vcf_path = self.config.variants_dir / f"{run_id}.vcf"
            if not vcf_path.exists():
                continue
            try:
                vcf = pysam.VariantFile(str(vcf_path))
                for record in vcf:
                    depth = record.info.get('DP', 0)
                    if depth < min_depth:
                        continue
                    ref_allele = record.ref
                    alt_alleles = record.alts
                    if not alt_alleles:
                        continue
                    for sample in record.samples.values():
                        ad = sample.get('AD', None)
                        if ad and len(ad) > 1:
                            ref_count = ad[0]
                            alt_count = ad[1]
                            total_count = ref_count + alt_count
                            if total_count > 0:
                                allele_freq = alt_count / total_count
                                if allele_freq >= min_freq:
                                    variants_list.append({
                                        'run_id': run_id,
                                        'patient_id': info['patient_id'],
                                        'timepoint': info['timepoint'],
                                        'chromosome': record.chrom,
                                        'position': record.pos,
                                        'ref': ref_allele,
                                        'alt': alt_alleles[0],
                                        'depth': depth,
                                        'ref_count': ref_count,
                                        'alt_count': alt_count,
                                        'allele_frequency': allele_freq
                                    })
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
        candidate_variants['timepoint_cat'] = pd.Categorical(
            candidate_variants['timepoint'], categories=timepoint_order, ordered=True)
        candidate_sorted = candidate_variants.sort_values('timepoint_cat')
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(candidate_sorted['timepoint_cat'], candidate_sorted['allele_frequency'],
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
        print(f"[OK] Visualization saved")
    
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
        patient_filter = self.config.config['data_source'].get('patient_id_filter', None)
        if patient_filter:
            self.sample_info = {k: v for k, v in self.sample_info.items() 
                               if v['patient_id'] == patient_filter}
            print(f"[INFO] Filtered to patient: {patient_filter}")
        for run_id in self.sample_info.keys():
            print(f"\nProcessing: {run_id}")
            if not self.download_sra_data(run_id):
                continue
            if not self.align_reads(run_id):
                continue
            if not self.convert_sort_index_bam(run_id):
                continue
            if not self.remove_duplicates(run_id):
                continue
            if not self.call_variants(run_id):
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
