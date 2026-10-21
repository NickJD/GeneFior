from pathlib import Path
import random

import pytest

from GeneFior.GeneFior_Architect import align_to_contig, evaluate_database_file, make_aligner
from GeneFior.GeneFior_Reconstruct import calculate_boundary_extension, trim_extension_unknowns
from GeneFior.GeneFior_Reconstruction_Visualise import count_reference_loci, make_aligner as visual_aligner


def unique_reference():
    rng = random.Random(481516)
    return "".join(rng.choice("ACGT") for _ in range(240))


def test_boundary_extension_uses_minimum_fraction_and_codon_rounding():
    assert calculate_boundary_extension(100, 0.20, True) == 21
    assert calculate_boundary_extension(100, 0.20, False) == 20
    assert calculate_boundary_extension(10, 0.01, True) == 3
    assert calculate_boundary_extension(10, 0.01, False) == 3


def test_extension_trimming_removes_terminal_unknowns_only():
    assert trim_extension_unknowns("NNACN", from_left=True) == "ACN"
    assert trim_extension_unknowns("NACNN", from_left=False) == "NAC"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_architect_finds_full_reference_on_either_strand(strand):
    query = unique_reference()
    oriented = query if strand == "+" else query.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    hit = align_to_contig(query, "GGGG" + oriented + "CCCC", make_aligner())
    assert hit is not None
    assert hit.identity_pct == pytest.approx(100)
    assert hit.query_coverage_pct == pytest.approx(100)
    assert hit.strand == strand
    assert (hit.target_start, hit.target_end) == (4, 4 + len(query))


def test_architect_reports_truncation_as_partial_query_coverage():
    query = unique_reference()
    hit = align_to_contig(query, query[:150], make_aligner())
    assert hit is not None
    assert hit.identity_pct == pytest.approx(100)
    assert hit.query_coverage_pct < 80


def test_architect_database_mode_compares_without_sample_assemblies(tmp_path: Path):
    reference = unique_reference()
    reference_fasta = tmp_path / "reference.fasta"
    reference_fasta.write_text(f">gene1\n{reference}\n")
    reconstruction = tmp_path / "gene1_consensus.fasta"
    reconstruction.write_text(f">gene1_sample_consensus\n{'G' * 40}{reference}{'C' * 40}\n")
    rows = evaluate_database_file(reconstruction, "gene1", "sample1",
                                  reference_fasta, make_aligner())
    assert len(rows) == 1
    assert rows[0]["comparison_mode"] == "database"
    assert rows[0]["status"] == "PASS"
    assert rows[0]["identity_pct"] == pytest.approx(100)
    assert rows[0]["database_reference_coverage_pct"] == pytest.approx(100)
    assert rows[0]["reconstruction_coverage_pct"] < 90


def test_assembly_locus_counter_finds_multiple_separate_copies(tmp_path: Path):
    reference = unique_reference()
    fasta = tmp_path / "assembly.fasta"
    fasta.write_text(f">contig1\n{'G' * 60}{reference}{'C' * 60}{reference}\n")
    loci = count_reference_loci(reference, fasta, visual_aligner())
    assert len(loci) == 2
    assert loci[0]["identity_pct"] == pytest.approx(100)
    assert loci[1]["identity_pct"] == pytest.approx(100)


def test_architect_finds_divergent_gene_on_chromosome_sized_contig():
    query = unique_reference()
    # A deterministic 20% substitution divergence defeats the original 19-mer
    # sparse seed strategy but should be recovered by dense short seeds.
    rng = random.Random(19)
    mutated = list(query)
    for i in range(0, len(mutated), 5):
        mutated[i] = next(base for base in "ACGT" if base != mutated[i])
    target = "G" * 100_000 + "".join(mutated) + "C" * 100_000
    hit = align_to_contig(query, target, make_aligner())
    assert hit is not None
    assert hit.query_coverage_pct >= 95
    assert hit.identity_pct >= 75
    assert 99_900 <= hit.target_start < 100_100


def test_reconstruction_pileup_counts_overlapping_paired_fragment_once():
    from GeneFior.GeneFior_Reconstruct import _build_pileup
    alignments = [
        {"name": "pair/1", "seq": "ACGT", "qual": "IIII", "cigar": "4M", "pos": 0, "mapq": 60, "flag": 0x41},
        {"name": "pair/2", "seq": "ACGT", "qual": "IIII", "cigar": "4M", "pos": 0, "mapq": 60, "flag": 0x81},
    ]
    bases, _ = _build_pileup(alignments, 4)
    assert all(sum(bases[pos].values()) == 1 for pos in range(4))


def test_reconstruction_pileup_excludes_low_quality_mapping():
    from GeneFior.GeneFior_Reconstruct import _build_pileup
    alignments = [{"name": "low", "seq": "ACGT", "qual": "!!!!", "cigar": "4M",
                   "pos": 0, "mapq": 0, "flag": 0}]
    bases, _ = _build_pileup(alignments, 4)
    assert not bases
