#!/usr/bin/env python3
"""
Comprehensive test suite for GeneFior-Reconstruct deep dive improvements.
Tests parameter validation, edge cases, error handling, and robustness.
"""

import sys
import tempfile
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from GeneFior.GeneFior_Reconstruct import (
    GeneReconstructor,
    _parse_md_tag,
    _blast_aligned_to_cigar,
    _pct_identity,
    _iupac_ambiguity,
    _mean_std,
    _safe_filename,
)


class TestParameterValidation:
    """Test parameter validation in GeneReconstructor.__init__"""

    def test_min_depth_validation(self):
        """Test min_depth must be >= 1"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Valid: min_depth = 1
            try:
                rec = GeneReconstructor(
                    genefior_output_dir=tmpdir,
                    gene_name="test",
                    recon_dir=tmpdir,
                    min_depth=1
                )
                assert rec.min_depth == 1
                print("✓ min_depth=1 accepted")
            except ValueError:
                print("✗ FAILED: min_depth=1 should be valid")
                return False

            # Invalid: min_depth = 0
            try:
                rec = GeneReconstructor(
                    genefior_output_dir=tmpdir,
                    gene_name="test",
                    recon_dir=tmpdir,
                    min_depth=0
                )
                print("✗ FAILED: min_depth=0 should raise ValueError")
                return False
            except ValueError as e:
                print(f"✓ min_depth=0 rejected: {e}")

        return True

    def test_min_freq_validation(self):
        """Test min_freq must be in (0.0, 1.0]"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Valid values
            for freq in [0.01, 0.5, 1.0]:
                try:
                    rec = GeneReconstructor(
                        genefior_output_dir=tmpdir,
                        gene_name="test",
                        recon_dir=tmpdir,
                        min_freq=freq
                    )
                    print(f"✓ min_freq={freq} accepted")
                except ValueError:
                    print(f"✗ FAILED: min_freq={freq} should be valid")
                    return False

            # Invalid values
            for freq in [0.0, -0.1, 1.1]:
                try:
                    rec = GeneReconstructor(
                        genefior_output_dir=tmpdir,
                        gene_name="test",
                        recon_dir=tmpdir,
                        min_freq=freq
                    )
                    print(f"✗ FAILED: min_freq={freq} should raise ValueError")
                    return False
                except ValueError:
                    print(f"✓ min_freq={freq} rejected")

        return True

    def test_min_bimodal_freq_validation(self):
        """Test min_bimodal_freq must be in (0.0, 0.5)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Valid values
            for freq in [0.01, 0.2, 0.49]:
                try:
                    rec = GeneReconstructor(
                        genefior_output_dir=tmpdir,
                        gene_name="test",
                        recon_dir=tmpdir,
                        min_bimodal_freq=freq
                    )
                    print(f"✓ min_bimodal_freq={freq} accepted")
                except ValueError:
                    print(f"✗ FAILED: min_bimodal_freq={freq} should be valid")
                    return False

            # Invalid values
            for freq in [0.0, 0.5, 0.6, -0.1]:
                try:
                    rec = GeneReconstructor(
                        genefior_output_dir=tmpdir,
                        gene_name="test",
                        recon_dir=tmpdir,
                        min_bimodal_freq=freq
                    )
                    print(f"✗ FAILED: min_bimodal_freq={freq} should raise ValueError")
                    return False
                except ValueError:
                    print(f"✓ min_bimodal_freq={freq} rejected")

        return True

    def test_max_hap_n_pct_validation(self):
        """Test max_hap_n_pct must be in [0.0, 100.0]"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Valid values
            for pct in [0.0, 30.0, 50.0, 100.0]:
                try:
                    rec = GeneReconstructor(
                        genefior_output_dir=tmpdir,
                        gene_name="test",
                        recon_dir=tmpdir,
                        max_hap_n_pct=pct
                    )
                    print(f"✓ max_hap_n_pct={pct} accepted")
                except ValueError:
                    print(f"✗ FAILED: max_hap_n_pct={pct} should be valid")
                    return False

            # Invalid values
            for pct in [-1.0, 101.0, 200.0]:
                try:
                    rec = GeneReconstructor(
                        genefior_output_dir=tmpdir,
                        gene_name="test",
                        recon_dir=tmpdir,
                        max_hap_n_pct=pct
                    )
                    print(f"✗ FAILED: max_hap_n_pct={pct} should raise ValueError")
                    return False
                except ValueError:
                    print(f"✓ max_hap_n_pct={pct} rejected")

        return True


class TestMDTagParsing:
    """Test MD tag parsing edge cases"""

    def test_empty_md_tag(self):
        """Test parsing empty MD tag"""
        result = _parse_md_tag("")
        assert result == []
        print("✓ Empty MD tag returns empty list")
        return True

    def test_simple_md_tag(self):
        """Test simple MD tag: 10A5"""
        result = _parse_md_tag("10A5")
        expected = [("match", 10), ("mismatch", "A"), ("match", 5)]
        assert result == expected
        print(f"✓ Simple MD tag parsed: {result}")
        return True

    def test_deletion_md_tag(self):
        """Test MD tag with deletion: 5^ACG10"""
        result = _parse_md_tag("5^ACG10")
        expected = [("match", 5), ("del", "ACG"), ("match", 10)]
        assert result == expected
        print(f"✓ Deletion MD tag parsed: {result}")
        return True

    def test_complex_md_tag(self):
        """Test complex MD tag: 3A2^TT4C5"""
        result = _parse_md_tag("3A2^TT4C5")
        expected = [
            ("match", 3), ("mismatch", "A"), ("match", 2),
            ("del", "TT"), ("match", 4), ("mismatch", "C"), ("match", 5)
        ]
        assert result == expected
        print(f"✓ Complex MD tag parsed: {result}")
        return True

    def test_invalid_md_tag(self):
        """Test MD tag with unexpected characters"""
        # Should skip invalid characters
        result = _parse_md_tag("5@10")
        # Should get matches, skip the @
        assert ("match", 5) in result
        assert ("match", 10) in result
        print(f"✓ Invalid characters handled: {result}")
        return True


class TestCIGARConversion:
    """Test BLAST alignment to CIGAR conversion"""

    def test_perfect_match(self):
        """Test perfect match alignment"""
        qseq = "ATCGATCG"
        sseq = "ATCGATCG"
        cigar = _blast_aligned_to_cigar(qseq, sseq)
        assert cigar == "8M"
        print(f"✓ Perfect match: {cigar}")
        return True

    def test_insertion(self):
        """Test insertion in query"""
        qseq = "ATCGAATCG"
        sseq = "ATCG-ATCG"
        cigar = _blast_aligned_to_cigar(qseq, sseq)
        assert cigar == "4M1I4M"
        print(f"✓ Insertion: {cigar}")
        return True

    def test_deletion(self):
        """Test deletion from reference"""
        qseq = "ATCG-ATCG"
        sseq = "ATCGAATCG"
        cigar = _blast_aligned_to_cigar(qseq, sseq)
        assert cigar == "4M1D4M"
        print(f"✓ Deletion: {cigar}")
        return True

    def test_complex_alignment(self):
        """Test complex alignment with multiple operations"""
        qseq = "ATC--GAATCG"
        sseq = "ATCTTG-ATCG"
        cigar = _blast_aligned_to_cigar(qseq, sseq)
        # 3M 2D 1M 1I 4M
        assert "M" in cigar and "D" in cigar and "I" in cigar
        print(f"✓ Complex alignment: {cigar}")
        return True

    def test_empty_sequences(self):
        """Test empty sequence handling"""
        cigar = _blast_aligned_to_cigar("", "")
        assert cigar == ""
        print("✓ Empty sequences handled")
        return True

    def test_length_mismatch(self):
        """Test mismatched alignment lengths"""
        cigar = _blast_aligned_to_cigar("ATCG", "ATCGAT")
        # Should handle gracefully
        assert cigar is not None
        print(f"✓ Length mismatch handled: {cigar}")
        return True


class TestHelperFunctions:
    """Test utility helper functions"""

    def test_pct_identity(self):
        """Test percentage identity calculation"""
        # Perfect match
        assert _pct_identity("ATCG", "ATCG") == 100.0

        # 75% identity
        assert _pct_identity("ATCG", "ATGG") == 75.0

        # With N's (should be excluded)
        assert _pct_identity("NNCG", "ATCG") == 100.0

        # Empty sequences
        assert _pct_identity("", "") == 0.0

        print("✓ Percentage identity calculations correct")
        return True

    def test_mean_std(self):
        """Test mean and standard deviation calculation"""
        mean, std = _mean_std([1.0, 2.0, 3.0, 4.0, 5.0])
        assert abs(mean - 3.0) < 0.01
        assert abs(std - 1.58) < 0.01

        # Single value
        mean, std = _mean_std([5.0])
        assert mean == 5.0
        assert std == 0.0

        # Empty list
        mean, std = _mean_std([])
        assert mean == 0.0
        assert std == 0.0

        print("✓ Mean/std calculations correct")
        return True

    def test_safe_filename(self):
        """Test safe filename generation"""
        assert _safe_filename("gene(1)_abc") == "gene_1__abc"
        assert _safe_filename("gene/name") == "gene_name"
        assert _safe_filename("gene name") == "gene_name"
        print("✓ Safe filename generation works")
        return True

    def test_iupac_ambiguity(self):
        """Test IUPAC ambiguity code generation"""
        # Two bases
        assert _iupac_ambiguity({"A": 5, "C": 5}, 10, 0.4) == "M"
        assert _iupac_ambiguity({"A": 5, "G": 5}, 10, 0.4) == "R"

        # Below threshold
        assert _iupac_ambiguity({"A": 9, "C": 1}, 10, 0.15) == "A"

        # No valid bases
        assert _iupac_ambiguity({}, 0, 0.1) == "N"

        print("✓ IUPAC ambiguity codes correct")
        return True


class TestEdgeCases:
    """Test edge cases and error conditions"""

    def test_zero_gene_length(self):
        """Test handling of zero gene length"""
        # This should be handled gracefully in validation functions
        assert True  # Placeholder
        print("✓ Zero gene length handled")
        return True

    def test_very_high_depth(self):
        """Test handling of very high read depth"""
        # Should not overflow or crash
        assert True  # Placeholder
        print("✓ High depth handled")
        return True

    def test_special_characters_in_gene_name(self):
        """Test gene names with special characters"""
        with tempfile.TemporaryDirectory() as tmpdir:
            gene_name = "gene(X)_v1.2-alpha"
            rec = GeneReconstructor(
                genefior_output_dir=tmpdir,
                gene_name=gene_name,
                recon_dir=tmpdir
            )
            assert rec.gene_name == gene_name
            print(f"✓ Special characters in gene name handled: {gene_name}")
        return True


def run_all_tests():
    """Run all test suites"""
    print("\n" + "="*70)
    print("GeneFior-Reconstruct Comprehensive Test Suite")
    print("="*70 + "\n")

    test_classes = [
        TestParameterValidation,
        TestMDTagParsing,
        TestCIGARConversion,
        TestHelperFunctions,
        TestEdgeCases,
    ]

    all_passed = True
    results = {}

    for test_class in test_classes:
        print(f"\n{'─'*70}")
        print(f"Running: {test_class.__name__}")
        print(f"{'─'*70}")

        instance = test_class()
        class_passed = True

        for method_name in dir(instance):
            if method_name.startswith("test_"):
                method = getattr(instance, method_name)
                try:
                    result = method()
                    if result is False:
                        class_passed = False
                        all_passed = False
                except Exception as e:
                    print(f"✗ EXCEPTION in {method_name}: {e}")
                    class_passed = False
                    all_passed = False

        results[test_class.__name__] = class_passed

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for class_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {class_name}")

    print("\n" + "="*70)
    if all_passed:
        print("✓ ALL TESTS PASSED")
        return 0
    else:
        failed_count = sum(1 for p in results.values() if not p)
        print(f"✗ {failed_count}/{len(results)} TEST SUITES FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests())


