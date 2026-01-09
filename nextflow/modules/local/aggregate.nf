/*
 * Aggregate Variants - Combine VCF files and generate summary CSV
 */

process AGGREGATE_VARIANTS {
    tag "aggregate"
    label 'process_low'
    
    publishDir "${params.outdir}/results", mode: params.publish_dir_mode
    
    input:
    path(vcf_files)
    
    output:
    path("all_variants.csv")       , emit: csv
    path("candidate_variants.csv") , emit: candidates
    path("summary_stats.txt")      , emit: stats
    
    script:
    """
    #!/usr/bin/env python3
    
    import gzip
    import csv
    import re
    from pathlib import Path
    from collections import defaultdict
    
    def parse_vcf(vcf_path):
        \"\"\"Parse VCF file and extract variant information.\"\"\"
        variants = []
        
        # Handle both gzipped and plain VCF
        if str(vcf_path).endswith('.gz'):
            fh = gzip.open(vcf_path, 'rt')
        else:
            fh = open(vcf_path, 'r')
        
        # Extract sample name from filename
        sample_id = Path(vcf_path).stem.replace('.ann', '').replace('.lofreq', '').replace('.vcf', '')
        
        for line in fh:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\\t')
            if len(fields) < 8:
                continue
            
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
                else:
                    info_dict[item] = True
            
            # Extract allele frequency
            af = float(info_dict.get('AF', 0))
            dp = int(info_dict.get('DP', 0))
            
            # Extract SnpEff annotation if present
            ann_field = info_dict.get('ANN', '')
            gene = ''
            effect = ''
            aa_change = ''
            
            if ann_field:
                ann_parts = ann_field.split('|')
                if len(ann_parts) >= 11:
                    alt_allele = ann_parts[0]
                    effect = ann_parts[1]
                    gene = ann_parts[3]
                    aa_change = ann_parts[10] if len(ann_parts) > 10 else ''
            
            variants.append({
                'sample_id': sample_id,
                'chromosome': chrom,
                'position': int(pos),
                'ref': ref,
                'alt': alt,
                'quality': float(qual) if qual != '.' else 0,
                'filter': filt,
                'depth': dp,
                'allele_frequency': af,
                'gene': gene,
                'effect': effect,
                'aa_change': aa_change
            })
        
        fh.close()
        return variants
    
    # Collect all VCF files
    vcf_files = list(Path('.').glob('*.vcf.gz'))
    
    # Parse all variants
    all_variants = []
    for vcf_file in vcf_files:
        # Skip index files
        if str(vcf_file).endswith('.tbi'):
            continue
        variants = parse_vcf(vcf_file)
        all_variants.extend(variants)
    
    # Write all variants to CSV
    with open('all_variants.csv', 'w', newline='') as f:
        if all_variants:
            writer = csv.DictWriter(f, fieldnames=all_variants[0].keys())
            writer.writeheader()
            writer.writerows(all_variants)
        else:
            f.write('No variants found\\n')
    
    # Identify candidate variants (recurrent, high-impact)
    variant_counts = defaultdict(list)
    for v in all_variants:
        key = (v['chromosome'], v['position'], v['ref'], v['alt'])
        variant_counts[key].append(v)
    
    # Candidates: variants seen in multiple samples or high-impact
    high_impact_effects = {'missense_variant', 'stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant'}
    candidates = []
    for key, var_list in variant_counts.items():
        # Check if high-impact or recurrent
        is_high_impact = any(v['effect'] in high_impact_effects for v in var_list)
        is_recurrent = len(var_list) > 1
        
        if is_high_impact or is_recurrent:
            # Use the variant with highest AF
            best = max(var_list, key=lambda x: x['allele_frequency'])
            best['n_samples'] = len(var_list)
            candidates.append(best)
    
    # Sort by allele frequency
    candidates.sort(key=lambda x: x['allele_frequency'], reverse=True)
    
    # Write candidate variants
    with open('candidate_variants.csv', 'w', newline='') as f:
        if candidates:
            fieldnames = list(candidates[0].keys())
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(candidates)
        else:
            f.write('No candidate variants identified\\n')
    
    # Write summary statistics
    with open('summary_stats.txt', 'w') as f:
        f.write('VCA Pipeline - Variant Summary Statistics\\n')
        f.write('=' * 50 + '\\n\\n')
        f.write(f'Total VCF files processed: {len(vcf_files)}\\n')
        f.write(f'Total variants detected: {len(all_variants)}\\n')
        f.write(f'Candidate variants (high-impact/recurrent): {len(candidates)}\\n\\n')
        
        # Per-sample summary
        samples = defaultdict(int)
        for v in all_variants:
            samples[v['sample_id']] += 1
        
        f.write('Variants per sample:\\n')
        for sample, count in sorted(samples.items()):
            f.write(f'  {sample}: {count}\\n')
        
        # Gene summary
        genes = defaultdict(int)
        for v in all_variants:
            if v['gene']:
                genes[v['gene']] += 1
        
        if genes:
            f.write('\\nVariants per gene:\\n')
            for gene, count in sorted(genes.items(), key=lambda x: -x[1])[:20]:
                f.write(f'  {gene}: {count}\\n')
    
    print(f'Processed {len(vcf_files)} VCF files')
    print(f'Found {len(all_variants)} total variants')
    print(f'Identified {len(candidates)} candidate variants')
    """
    
    stub:
    """
    echo "sample_id,chromosome,position,ref,alt,quality,filter,depth,allele_frequency,gene,effect,aa_change" > all_variants.csv
    echo "sample_id,chromosome,position,ref,alt,quality,filter,depth,allele_frequency,gene,effect,aa_change,n_samples" > candidate_variants.csv
    echo "VCA Pipeline - Summary Stats" > summary_stats.txt
    """
}
