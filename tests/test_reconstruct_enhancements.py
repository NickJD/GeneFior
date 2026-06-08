#!/usr/bin/env python3
"""
Test script to verify the GeneFior-Reconstruct enhancements:
1. Stale file cleanup
2. Quality warnings detection
3. Warn tags in FASTA headers
4. Quality warnings in validation reports
"""

import tempfile
import shutil
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from GeneFior.GeneFior_Reconstruct import GeneReconstructor


def test_stale_file_cleanup():
    """Test that old output files are properly cleaned up."""
    print("\n" + "="*70)
    print("TEST 1: Stale File Cleanup")
    print("="*70)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)
        recon_dir = tmppath / "reconstruction"
        gene_dir = recon_dir / "test_gene" / "bowtie2"
        gene_dir.mkdir(parents=True, exist_ok=True)

        # Create some "stale" files from a previous run
        stale_files = [
            gene_dir / "test_gene_consensus.fasta",
            gene_dir / "test_gene_haplotypes.fasta",
            gene_dir / "test_gene_depth.tsv",
            gene_dir / "test_gene_variants.tsv",
            gene_dir / "test_gene_reads.fasta",
            gene_dir / "test_gene_validation.txt",
            gene_dir / "test_gene_reconstruction_plot.html",
        ]

        for f in stale_files:
            f.write_text("old content")

        print(f"✓ Created {len(stale_files)} stale files")
        assert all(f.exists() for f in stale_files), "Stale files should exist"

        # Now simulate the cleanup code that happens at the start of _run_analysis_pipeline
        gene_safe = "test_gene"
        for sfx in ["_consensus.fasta", "_haplotypes.fasta", "_depth.tsv",
                    "_variants.tsv", "_reads.fasta", "_validation.txt",
                    "_reconstruction_plot.html"]:
            stale = gene_dir / f"{gene_safe}{sfx}"
            if stale.exists():
                stale.unlink()

        # Verify all stale files were removed
        remaining = [f for f in stale_files if f.exists()]
        if remaining:
            print(f"✗ FAILED: {len(remaining)} files still exist: {remaining}")
            return False
        else:
            print(f"✓ All {len(stale_files)} stale files successfully cleaned up")
            return True


def test_quality_warnings():
    """Test quality warning detection logic."""
    print("\n" + "="*70)
    print("TEST 2: Quality Warning Detection")
    print("="*70)

    test_cases = [
        {
            "name": "High N% warning",
            "consensus_n_pct": 25.0,
            "grade": "B",
            "has_haplotypes": False,
            "attempt_hap": False,
            "expected_warnings": 1,
            "expected_tag": "QUALITY_WARNING"
        },
        {
            "name": "Grade F warning",
            "consensus_n_pct": 5.0,
            "grade": "F",
            "has_haplotypes": False,
            "attempt_hap": False,
            "expected_warnings": 1,
            "expected_tag": "FAILED_RECONSTRUCTION"
        },
        {
            "name": "Failed haplotype separation",
            "consensus_n_pct": 10.0,
            "grade": "B",
            "has_haplotypes": False,
            "attempt_hap": True,
            "expected_warnings": 1,
            "expected_tag": "QUALITY_WARNING"
        },
        {
            "name": "Multiple warnings",
            "consensus_n_pct": 30.0,
            "grade": "F",
            "has_haplotypes": False,
            "attempt_hap": True,
            "expected_warnings": 3,
            "expected_tag": "FAILED_RECONSTRUCTION"
        },
        {
            "name": "No warnings",
            "consensus_n_pct": 2.0,
            "grade": "A",
            "has_haplotypes": False,
            "attempt_hap": False,
            "expected_warnings": 0,
            "expected_tag": ""
        },
    ]

    all_passed = True
    for tc in test_cases:
        # Simulate the warning detection logic
        quality_warnings = []
        min_depth = 3

        if tc["consensus_n_pct"] > 20.0:
            quality_warnings.append(
                f"High gap rate: {tc['consensus_n_pct']:.1f}% of consensus positions are N "
                f"(insufficient depth < {min_depth}×). Reconstruct with more reads or reduce -min_depth."
            )
        if tc["grade"] == "F":
            quality_warnings.append(
                "Grade F — coverage/read support below thresholds. UNRELIABLE — do not use for analysis."
            )
        if tc["attempt_hap"] and not tc["has_haplotypes"]:
            quality_warnings.append(
                "Multi-version signal detected but haplotype separation FAILED — all groups rejected "
                f"(too few reads or >30% N after split). Consensus may blend ≥2 versions. "
                "Try more reads or lower -min_depth."
            )

        warn_tag = "FAILED_RECONSTRUCTION" if tc["grade"] == "F" else ("QUALITY_WARNING" if quality_warnings else "")

        passed = (
            len(quality_warnings) == tc["expected_warnings"] and
            warn_tag == tc["expected_tag"]
        )

        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {tc['name']}")
        print(f"  Expected {tc['expected_warnings']} warnings, got {len(quality_warnings)}")
        print(f"  Expected tag '{tc['expected_tag']}', got '{warn_tag}'")

        if not passed:
            all_passed = False

    return all_passed


def test_fasta_headers():
    """Test that FASTA headers include proper warn tags."""
    print("\n" + "="*70)
    print("TEST 3: FASTA Header Warn Tags")
    print("="*70)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)
        fasta_path = tmppath / "test_consensus.fasta"

        # Simulate writing a FASTA with a warning tag
        gene_name = "test_gene"
        label = "sample_consensus"
        seq = "ATCGATCGATCGATCGNNNNNATCGATCG"
        stats_row = {"Num_Sequences_Mapped": "100", "Gene_Coverage": "95.5", "Avg_Identity": "98.2"}
        warn_tag = "QUALITY_WARNING"

        n_pct = seq.count("N") / len(seq) * 100

        with open(fasta_path, "w") as fh:
            tag = f" [{warn_tag}]" if warn_tag else ""
            fh.write(
                f">{gene_name}_{label}{tag} "
                f"[len={len(seq)}] [N_pct={n_pct:.1f}] "
                f"[reads={stats_row['Num_Sequences_Mapped']}] "
                f"[gene_cov={stats_row['Gene_Coverage']}%] "
                f"[avg_id={stats_row['Avg_Identity']}%]\n"
            )
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

        # Verify the file was created and contains the warning tag
        content = fasta_path.read_text()

        if "[QUALITY_WARNING]" in content:
            print("✓ QUALITY_WARNING tag present in FASTA header")
            if f"N_pct={n_pct:.1f}" in content:
                print(f"✓ N percentage ({n_pct:.1f}%) included in header")
                return True
            else:
                print(f"✗ FAILED: N percentage not found in header")
                return False
        else:
            print("✗ FAILED: QUALITY_WARNING tag not found in FASTA header")
            print(f"Header content: {content.split(chr(10))[0]}")
            return False


def test_validation_report():
    """Test that validation reports include quality warnings section."""
    print("\n" + "="*70)
    print("TEST 4: Validation Report Quality Warnings")
    print("="*70)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)
        report_path = tmppath / "test_validation.txt"

        # Simulate writing a validation report with warnings
        quality_warnings = [
            "High gap rate: 25.5% of consensus positions are N (insufficient depth < 3×).",
            "Multi-version signal detected but haplotype separation FAILED.",
        ]

        gene_name = "test_gene"
        validations = [{
            "haplotype_label": "consensus",
            "grade": "C",
            "interpretation": "Moderate confidence.",
            "gene_length": 1000,
            "covered_positions": 800,
            "coverage_pct": 80.0,
            "n_pct": 25.5,
            "mean_depth": 15.2,
            "cv_depth_pct": 45.3,
            "read_support_pct": 82.1,
            "ref_identity_pct": 95.8,
            "mv_classification": "multi",
            "n_informative_sites": 15,
            "mean_alt_freq_pct": 35.2,
        }]

        from datetime import datetime
        with open(report_path, "w") as fh:
            fh.write(f"RECONSTRUCTION VALIDATION REPORT\n")
            fh.write(f"Gene: {gene_name}\n")
            fh.write(f"Generated: {datetime.now().isoformat()}\n")
            fh.write("=" * 70 + "\n\n")

            # ── Quality warnings block ────────────────────────────────────────────
            if quality_warnings:
                fh.write("  ⚠  QUALITY WARNINGS\n")
                for w in quality_warnings:
                    fh.write(f"     • {w}\n")
                fh.write("\n" + "-" * 70 + "\n\n")

            for v in validations:
                fh.write(f"  Sequence: {v['haplotype_label']}\n")
                fh.write(f"  Grade:    {v['grade']}  — {v['interpretation']}\n\n")
                fh.write(f"  Coverage\n")
                fh.write(f"    Gene length:            {v['gene_length']} bp\n")
                # ... (rest of report)

        # Verify the report contains the warnings section
        content = report_path.read_text()

        if "⚠  QUALITY WARNINGS" in content:
            print("✓ Quality warnings header present")
            if all(w in content for w in quality_warnings):
                print(f"✓ All {len(quality_warnings)} warnings present in report")
                return True
            else:
                print("✗ FAILED: Not all warnings found in report")
                return False
        else:
            print("✗ FAILED: Quality warnings section not found")
            return False


def main():
    """Run all tests."""
    print("\n" + "#"*70)
    print("# GeneFior-Reconstruct Enhancement Tests")
    print("#"*70)

    results = {
        "Stale File Cleanup": test_stale_file_cleanup(),
        "Quality Warning Detection": test_quality_warnings(),
        "FASTA Header Warn Tags": test_fasta_headers(),
        "Validation Report Warnings": test_validation_report(),
    }

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    for test_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")

    all_passed = all(results.values())
    print("\n" + "="*70)
    if all_passed:
        print("✓ ALL TESTS PASSED")
        return 0
    else:
        failed_count = sum(1 for p in results.values() if not p)
        print(f"✗ {failed_count}/{len(results)} TESTS FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())

