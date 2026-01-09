#!/usr/bin/env python
"""
Quick validation test for BQSR and Germline Filtering implementations.
Tests the logic without requiring full pipeline execution.
"""
import sys
import subprocess
from pathlib import Path

# Add parent dir to path
sys.path.insert(0, str(Path(__file__).parent))

def test_lofreq_viterbi_available():
    """Test that LoFreq Viterbi command is available."""
    print("=" * 60)
    print("TEST 1: LoFreq Viterbi Availability")
    print("=" * 60)
    try:
        result = subprocess.run(
            ["lofreq", "viterbi"],
            capture_output=True, text=True
        )
        # lofreq viterbi without args returns usage (non-zero exit)
        if "Usage: lofreq viterbi" in result.stderr or "Usage: lofreq viterbi" in result.stdout:
            print("✅ PASS: lofreq viterbi is available")
            return True
        else:
            print("❌ FAIL: lofreq viterbi not found")
            print(f"   stderr: {result.stderr[:200]}")
            return False
    except Exception as e:
        print(f"❌ FAIL: Error running lofreq viterbi: {e}")
        return False


def test_bcftools_available():
    """Test that bcftools is available."""
    print("\n" + "=" * 60)
    print("TEST 2: bcftools Availability")
    print("=" * 60)
    try:
        result = subprocess.run(
            ["bcftools", "--version"],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            version = result.stdout.split('\n')[0]
            print(f"✅ PASS: {version}")
            return True
        else:
            print("❌ FAIL: bcftools not working")
            return False
    except Exception as e:
        print(f"❌ FAIL: Error running bcftools: {e}")
        return False


def test_germline_filter_logic():
    """Test germline filtering logic with test VCF files."""
    print("\n" + "=" * 60)
    print("TEST 3: Germline Filter Logic")
    print("=" * 60)
    
    test_data_dir = Path(__file__).parent.parent / "nextflow" / "test" / "data"
    test_vcf = test_data_dir / "test_sample.ann.vcf.gz"
    gnomad_vcf = test_data_dir / "gnomad_test_kras.vcf.gz"
    
    if not test_vcf.exists():
        print(f"⚠️  SKIP: Test VCF not found: {test_vcf}")
        return None
    
    if not gnomad_vcf.exists():
        print(f"⚠️  SKIP: gnomAD VCF not found: {gnomad_vcf}")
        return None
    
    output_annotated = test_data_dir / "validation_test.annotated.vcf.gz"
    output_filtered = test_data_dir / "validation_test.somatic.vcf.gz"
    
    try:
        # Step 1: Annotate
        print("  Annotating with gnomAD...")
        subprocess.run([
            "bcftools", "annotate",
            "-a", str(gnomad_vcf),
            "-c", "INFO/gnomAD_AF:=INFO/AF",
            "-O", "z",
            "-o", str(output_annotated),
            str(test_vcf)
        ], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "vcf", str(output_annotated)], check=True, capture_output=True)
        
        # Step 2: Filter
        print("  Filtering germline variants (AF > 0.01)...")
        subprocess.run([
            "bcftools", "filter",
            "-e", "INFO/gnomAD_AF > 0.01",
            "-s", "GERMLINE",
            "-m", "+",
            "-O", "z",
            "-o", str(output_filtered),
            str(output_annotated)
        ], check=True, capture_output=True)
        subprocess.run(["tabix", "-p", "vcf", str(output_filtered)], check=True, capture_output=True)
        
        # Count results
        before_result = subprocess.run(
            ["bcftools", "view", "-H", str(test_vcf)],
            capture_output=True, text=True
        )
        before_count = before_result.stdout.count('\n')
        
        pass_result = subprocess.run(
            ["bcftools", "view", "-f", "PASS", "-H", str(output_filtered)],
            capture_output=True, text=True
        )
        pass_count = pass_result.stdout.count('\n')
        
        germline_result = subprocess.run(
            ["bcftools", "view", "-f", "GERMLINE", "-H", str(output_filtered)],
            capture_output=True, text=True
        )
        germline_count = germline_result.stdout.count('\n')
        
        # Cleanup
        output_annotated.unlink(missing_ok=True)
        Path(str(output_annotated) + '.tbi').unlink(missing_ok=True)
        output_filtered.unlink(missing_ok=True)
        Path(str(output_filtered) + '.tbi').unlink(missing_ok=True)
        
        # Validate
        print(f"\n  Results:")
        print(f"    Variants before filtering: {before_count}")
        print(f"    Variants after (PASS):     {pass_count}")
        print(f"    Filtered as GERMLINE:      {germline_count}")
        
        if before_count == 7 and pass_count == 4 and germline_count == 3:
            print("\n✅ PASS: Germline filtering works correctly")
            return True
        else:
            print(f"\n❌ FAIL: Expected 7 before, 4 pass, 3 germline")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"❌ FAIL: Command failed: {e.stderr[:200] if e.stderr else str(e)}")
        return False
    except Exception as e:
        print(f"❌ FAIL: Error: {str(e)[:200]}")
        return False


def test_python_module_imports():
    """Test that Python module imports work."""
    print("\n" + "=" * 60)
    print("TEST 4: Python Module Imports")  
    print("=" * 60)
    try:
        from run_vca_pipeline import VCAPipeline, VCAConfig
        print("✅ PASS: VCAPipeline and VCAConfig imported successfully")
        
        # Check for new methods
        if hasattr(VCAPipeline, 'recalibrate_base_qualities'):
            print("✅ PASS: recalibrate_base_qualities method exists")
        else:
            print("❌ FAIL: recalibrate_base_qualities method not found")
            return False
            
        if hasattr(VCAPipeline, 'filter_germline_variants'):
            print("✅ PASS: filter_germline_variants method exists")
        else:
            print("❌ FAIL: filter_germline_variants method not found")
            return False
            
        if hasattr(VCAPipeline, 'download_gnomad'):
            print("✅ PASS: download_gnomad method exists")
        else:
            print("❌ FAIL: download_gnomad method not found")
            return False
            
        return True
    except Exception as e:
        print(f"❌ FAIL: Import error: {e}")
        return False


def test_config_yaml():
    """Test that config.yaml has the new sections."""
    print("\n" + "=" * 60)
    print("TEST 5: Config YAML Structure")
    print("=" * 60)
    try:
        import yaml
        config_path = Path(__file__).parent / "config.yaml"
        
        with open(config_path) as f:
            config = yaml.safe_load(f)
        
        has_germline = 'germline_filtering' in config
        has_bqsr = 'bqsr' in config
        
        if has_germline:
            print("✅ PASS: germline_filtering section exists")
            print(f"         enabled: {config['germline_filtering'].get('enabled', 'N/A')}")
            print(f"         max_population_af: {config['germline_filtering'].get('max_population_af', 'N/A')}")
        else:
            print("❌ FAIL: germline_filtering section not found")
            
        if has_bqsr:
            print("✅ PASS: bqsr section exists")
            print(f"         enabled: {config['bqsr'].get('enabled', 'N/A')}")
            print(f"         method: {config['bqsr'].get('method', 'N/A')}")
        else:
            print("❌ FAIL: bqsr section not found")
            
        return has_germline and has_bqsr
    except Exception as e:
        print(f"❌ FAIL: Error reading config: {e}")
        return False


def main():
    """Run all validation tests."""
    print("\n" + "=" * 60)
    print("   BQSR & GERMLINE FILTER VALIDATION TESTS")
    print("=" * 60 + "\n")
    
    results = []
    
    results.append(("LoFreq Viterbi", test_lofreq_viterbi_available()))
    results.append(("bcftools", test_bcftools_available()))
    results.append(("Germline Filter Logic", test_germline_filter_logic()))
    results.append(("Python Module Imports", test_python_module_imports()))
    results.append(("Config YAML", test_config_yaml()))
    
    # Summary
    print("\n" + "=" * 60)
    print("   VALIDATION SUMMARY")
    print("=" * 60)
    
    passed = 0
    failed = 0
    skipped = 0
    
    for name, result in results:
        if result is True:
            status = "✅ PASS"
            passed += 1
        elif result is False:
            status = "❌ FAIL"
            failed += 1
        else:
            status = "⚠️  SKIP"
            skipped += 1
        print(f"  {status}: {name}")
    
    print("\n" + "-" * 60)
    print(f"  Total: {passed} passed, {failed} failed, {skipped} skipped")
    print("=" * 60 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
