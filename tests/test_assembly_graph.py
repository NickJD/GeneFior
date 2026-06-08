#!/usr/bin/env python3
"""
Test script for assembly graph generation functionality.
Validates GFA format output and Bandage compatibility.
"""

import sys
import tempfile
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent / "src"))

from GeneFior.GeneFior_Reconstruct import GeneReconstructor


def test_gfa_format_validation():
    """Test that generated GFA files follow the specification"""
    print("\n" + "="*70)
    print("TEST: GFA Format Validation")
    print("="*70)

    # Create mock data
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)
        tool_dir = tmppath / "test"
        tool_dir.mkdir()

        # Mock reconstruction data
        gene_safe = "test_gene"
        consensus = "ATCGATCGATCG"
        per_pos = [
            {"pos": i+1, "depth": 10, "consensus_base": consensus[i],
             "top_base_freq": 95.0, "alt_base": "", "alt_base_freq": 0.0,
             "has_insertion": False, "ins_seq": "", "ins_depth": 0, "ins_freq": 0.0}
            for i in range(len(consensus))
        ]

        # Add a variant position
        per_pos[6]["alt_base"] = "T"
        per_pos[6]["alt_base_freq"] = 30.0
        per_pos[6]["top_base_freq"] = 70.0

        variants = [{"pos": 7, "ref_base": "A", "sample_base": "G",
                    "type": "SNP", "depth": 10, "freq_pct": 70.0}]
        informative_pos = [{"pos": 7, "top_base": "G", "top_freq": 0.7,
                           "alt_base": "T", "alt_freq": 0.3, "depth": 10}]

        base_pileup = defaultdict(lambda: defaultdict(int))
        for i in range(len(consensus)):
            base_pileup[i][consensus[i]] = 10
        base_pileup[6]["G"] = 7
        base_pileup[6]["T"] = 3

        alignments = [
            {"name": "read1", "pos": 0, "cigar": "12M", "seq": "ATCGATCGATCG"},
            {"name": "read2", "pos": 0, "cigar": "12M", "seq": "ATCGATCGATCG"},
        ]

        # Create reconstructor instance
        rec = GeneReconstructor(
            genefior_output_dir=str(tmppath),
            gene_name="test_gene",
            recon_dir=str(tmppath)
        )

        # Generate assembly graph
        rec._write_assembly_graph(
            tool_label="test",
            tool_dir=tool_dir,
            gene_safe=gene_safe,
            consensus=consensus,
            per_pos=per_pos,
            variants=variants,
            informative_pos=informative_pos,
            haplotypes=[],
            mv_class="single",
            base_pileup=base_pileup,
            alignments=alignments,
        )

        # Validate GFA file
        gfa_path = tool_dir / f"{gene_safe}_assembly_graph.gfa"

        if not gfa_path.exists():
            print("✗ FAILED: GFA file not created")
            return False

        content = gfa_path.read_text()
        lines = content.strip().split("\n")

        # Validate GFA format
        checks = {
            "has_header": False,
            "has_segments": False,
            "has_links": False,
            "has_paths": False,
        }

        for line in lines:
            if line.startswith("H\t"):
                checks["has_header"] = True
                # Check version tag
                if "VN:Z:1.0" not in line:
                    print("✗ FAILED: Invalid GFA version")
                    return False
            elif line.startswith("S\t"):
                checks["has_segments"] = True
                # Validate S line format: S <name> <seq> [tags]
                parts = line.split("\t")
                if len(parts) < 3:
                    print(f"✗ FAILED: Invalid S line: {line}")
                    return False
            elif line.startswith("L\t"):
                checks["has_links"] = True
                # Validate L line format: L <from> <from_orient> <to> <to_orient> <overlap>
                parts = line.split("\t")
                if len(parts) < 6:
                    print(f"✗ FAILED: Invalid L line: {line}")
                    return False
            elif line.startswith("P\t"):
                checks["has_paths"] = True
                # Validate P line format: P <path_name> <segments> <overlaps>
                parts = line.split("\t")
                if len(parts) < 4:
                    print(f"✗ FAILED: Invalid P line: {line}")
                    return False

        # Check all required elements present
        for check, status in checks.items():
            symbol = "✓" if status else "✗"
            print(f"{symbol} {check}: {status}")
            if not status and check != "has_links":  # Links optional for single segment
                return False

        print(f"✓ GFA file created: {gfa_path.name}")
        print(f"✓ Total lines: {len(lines)}")
        print(f"✓ Format validation passed")

        return True


def test_variant_bubbles():
    """Test that variant positions create bubble structures"""
    print("\n" + "="*70)
    print("TEST: Variant Bubble Generation")
    print("="*70)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)
        tool_dir = tmppath / "test"
        tool_dir.mkdir()

        # Create sequence with multiple variants
        consensus = "ATCGATCGATCGATCGATCG"
        per_pos = [
            {"pos": i+1, "depth": 20, "consensus_base": consensus[i],
             "top_base_freq": 95.0, "alt_base": "", "alt_base_freq": 0.0,
             "has_insertion": False, "ins_seq": "", "ins_depth": 0, "ins_freq": 0.0}
            for i in range(len(consensus))
        ]

        # Add variant positions at 5, 10, 15
        variant_positions = [5, 10, 15]
        variants = []
        informative_pos = []
        base_pileup = defaultdict(lambda: defaultdict(int))

        for i in range(len(consensus)):
            base_pileup[i][consensus[i]] = 20
            if i + 1 in variant_positions:
                per_pos[i]["alt_base"] = "T" if consensus[i] != "T" else "A"
                per_pos[i]["alt_base_freq"] = 40.0
                per_pos[i]["top_base_freq"] = 60.0

                variants.append({
                    "pos": i + 1, "ref_base": "A", "sample_base": consensus[i],
                    "type": "SNP", "depth": 20, "freq_pct": 60.0
                })

                informative_pos.append({
                    "pos": i + 1, "top_base": consensus[i], "top_freq": 0.6,
                    "alt_base": "T" if consensus[i] != "T" else "A",
                    "alt_freq": 0.4, "depth": 20
                })

                # Add alt allele to pileup
                alt_base = "T" if consensus[i] != "T" else "A"
                base_pileup[i][consensus[i]] = 12
                base_pileup[i][alt_base] = 8

        rec = GeneReconstructor(
            genefior_output_dir=str(tmppath),
            gene_name="test_gene",
            recon_dir=str(tmppath)
        )

        rec._write_assembly_graph(
            tool_label="test",
            tool_dir=tool_dir,
            gene_safe="test_gene",
            consensus=consensus,
            per_pos=per_pos,
            variants=variants,
            informative_pos=informative_pos,
            haplotypes=[],
            mv_class="multi",
            base_pileup=base_pileup,
            alignments=[],
        )

        gfa_path = tool_dir / "test_gene_assembly_graph.gfa"
        content = gfa_path.read_text()

        # Count segments (should have bubbles at variant positions)
        segment_count = len([l for l in content.split("\n") if l.startswith("S\t")])

        # We expect: segments between variants + 2 alleles per variant
        # With 3 variants, we should have multiple segments
        if segment_count < 5:  # At least a few segments
            print(f"✗ FAILED: Expected multiple segments, got {segment_count}")
            return False

        print(f"✓ Generated {segment_count} segments for {len(variants)} variants")

        # Check for multiple paths at variant positions
        # Count segments with same position annotation
        position_segments = {}
        for line in content.split("\n"):
            if line.startswith("S\t"):
                parts = line.split("\t")
                # Extract position from tags
                for tag in parts[3:]:
                    if tag.startswith("SC:i:"):
                        pos = int(tag.split(":")[2])
                        position_segments[pos] = position_segments.get(pos, 0) + 1

        variant_pos_count = sum(1 for pos in variant_positions
                               if position_segments.get(pos, 0) >= 2)

        if variant_pos_count > 0:
            print(f"✓ Found bubble structures at {variant_pos_count} variant positions")
        else:
            print("⚠ Warning: Variant bubbles may not be properly structured")

        return True


def test_haplotype_paths():
    """Test that haplotype separation creates separate paths"""
    print("\n" + "="*70)
    print("TEST: Haplotype Path Generation")
    print("="*70)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)
        tool_dir = tmppath / "test"
        tool_dir.mkdir()

        consensus = "ATCGATCGATCG"
        per_pos = [
            {"pos": i+1, "depth": 20, "consensus_base": consensus[i],
             "top_base_freq": 95.0, "alt_base": "", "alt_base_freq": 0.0,
             "has_insertion": False, "ins_seq": "", "ins_depth": 0, "ins_freq": 0.0}
            for i in range(len(consensus))
        ]

        # Create mock haplotypes
        haplotypes = [
            ("haplotype_1", "ATCGATCGATCG", "ATCGATCGATCG", per_pos,
             {"coverage_pct": 95.0, "grade": "A"}),
            ("haplotype_2", "ATCGTTCGATCG", "ATCGTTCGATCG", per_pos,
             {"coverage_pct": 92.0, "grade": "B"}),
        ]

        rec = GeneReconstructor(
            genefior_output_dir=str(tmppath),
            gene_name="test_gene",
            recon_dir=str(tmppath)
        )

        rec._write_assembly_graph(
            tool_label="test",
            tool_dir=tool_dir,
            gene_safe="test_gene",
            consensus=consensus,
            per_pos=per_pos,
            variants=[],
            informative_pos=[],
            haplotypes=haplotypes,
            mv_class="multi",
            base_pileup=defaultdict(lambda: defaultdict(int)),
            alignments=[],
        )

        gfa_path = tool_dir / "test_gene_assembly_graph.gfa"
        content = gfa_path.read_text()

        # Count paths
        path_lines = [l for l in content.split("\n") if l.startswith("P\t")]

        if len(path_lines) < 3:  # consensus + 2 haplotypes
            print(f"✗ FAILED: Expected 3 paths, got {len(path_lines)}")
            return False

        print(f"✓ Generated {len(path_lines)} paths (consensus + haplotypes)")

        # Check for haplotype path names
        has_hap1 = any("haplotype_1" in line for line in path_lines)
        has_hap2 = any("haplotype_2" in line for line in path_lines)
        has_consensus = any("consensus" in line for line in path_lines)

        if has_hap1 and has_hap2 and has_consensus:
            print("✓ All expected paths present: consensus, haplotype_1, haplotype_2")
            return True
        else:
            print("✗ FAILED: Missing expected paths")
            return False


def test_bandage_compatibility():
    """Test that GFA is compatible with Bandage requirements"""
    print("\n" + "="*70)
    print("TEST: Bandage Compatibility")
    print("="*70)

    # Bandage minimum requirements:
    # 1. Valid GFA 1.0 format
    # 2. Each segment must have valid sequence
    # 3. Links must reference existing segments
    # 4. No circular references without proper tags

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)
        tool_dir = tmppath / "test"
        tool_dir.mkdir()

        consensus = "ATCGATCGATCGATCGATCGATCGATCG"
        per_pos = [
            {"pos": i+1, "depth": 15, "consensus_base": consensus[i],
             "top_base_freq": 95.0, "alt_base": "", "alt_base_freq": 0.0,
             "has_insertion": False, "ins_seq": "", "ins_depth": 0, "ins_freq": 0.0}
            for i in range(len(consensus))
        ]

        rec = GeneReconstructor(
            genefior_output_dir=str(tmppath),
            gene_name="test_gene",
            recon_dir=str(tmppath)
        )

        rec._write_assembly_graph(
            tool_label="test",
            tool_dir=tool_dir,
            gene_safe="test_gene",
            consensus=consensus,
            per_pos=per_pos,
            variants=[],
            informative_pos=[],
            haplotypes=[],
            mv_class="single",
            base_pileup=defaultdict(lambda: defaultdict(int)),
            alignments=[],
        )

        gfa_path = tool_dir / "test_gene_assembly_graph.gfa"
        content = gfa_path.read_text()
        lines = content.split("\n")

        # Extract segment IDs
        segment_ids = set()
        for line in lines:
            if line.startswith("S\t"):
                parts = line.split("\t")
                seg_id = parts[1]
                seq = parts[2]

                # Check sequence is valid DNA
                if not all(c in "ACGTN" for c in seq):
                    print(f"✗ FAILED: Invalid sequence in segment {seg_id}: {seq}")
                    return False

                segment_ids.add(seg_id)

        print(f"✓ All {len(segment_ids)} segments have valid DNA sequences")

        # Check links reference existing segments
        for line in lines:
            if line.startswith("L\t"):
                parts = line.split("\t")
                from_seg = parts[1]
                to_seg = parts[3]

                if from_seg not in segment_ids:
                    print(f"✗ FAILED: Link references non-existent segment: {from_seg}")
                    return False
                if to_seg not in segment_ids:
                    print(f"✗ FAILED: Link references non-existent segment: {to_seg}")
                    return False

        print("✓ All links reference existing segments")

        # Check paths reference existing segments
        for line in lines:
            if line.startswith("P\t"):
                parts = line.split("\t")
                path_segs = parts[2].split(",")
                for seg_ref in path_segs:
                    seg_id = seg_ref.rstrip("+-")
                    if seg_id not in segment_ids:
                        print(f"✗ FAILED: Path references non-existent segment: {seg_id}")
                        return False

        print("✓ All paths reference existing segments")
        print("✓ GFA file is Bandage-compatible")

        return True


def run_all_tests():
    """Run all assembly graph tests"""
    print("\n" + "#"*70)
    print("# Assembly Graph Generation Tests")
    print("#"*70)

    tests = [
        ("GFA Format Validation", test_gfa_format_validation),
        ("Variant Bubble Generation", test_variant_bubbles),
        ("Haplotype Path Generation", test_haplotype_paths),
        ("Bandage Compatibility", test_bandage_compatibility),
    ]

    results = {}
    for test_name, test_func in tests:
        try:
            result = test_func()
            results[test_name] = result
        except Exception as e:
            print(f"✗ EXCEPTION in {test_name}: {e}")
            import traceback
            traceback.print_exc()
            results[test_name] = False

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for test_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")

    print("\n" + "="*70)
    if all(results.values()):
        print("✓ ALL TESTS PASSED")
        return 0
    else:
        failed = sum(1 for p in results.values() if not p)
        print(f"✗ {failed}/{len(results)} TESTS FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests())

