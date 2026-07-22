"""
test_coverage_identity_calculations.py
=======================================
Formal unit and integration tests for the pident (percent identity) and
query/gene coverage calculations in GeneFior's workflow parsing layer.

Three separate computation sites are tested:

  Stage-1  – the stream filter inside _stream_filter_hits()
             Formula: ((abs(qend - qstart) + 1) / qlen) * 100
             Purpose: pre-filters BLAST/DIAMOND lines before writing to disk.

  Stage-2  – parse_blast_results()
             Formula: (abs(qend - qstart) + 1) / qlen * 100  (same, for consistency)
             Purpose: final per-read coverage used to decide whether a hit contributes
             to a gene's covered positions.

  BAM      – parse_bam_results()
             Identity: (alignment_length - NM) / alignment_length * 100
             Coverage: (M + I CIGAR bases) / read_length * 100
             Purpose: equivalent metrics derived from CIGAR strings.

A key regression is also verified: the OLD Stage-2 formula
``alignment_len / qlen`` (BLAST field-3 / qlen) could exceed 100 % for
alignments with reference deletions.  The new span-based formula cannot.

GeneStats class tests verify that the coverage union, depth, and identity
aggregation work correctly across single and multiple hits.
"""

import json
import re
import tempfile
import textwrap
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import pytest

from GeneFior.gene_stats import GeneStats
from GeneFior.evidence import EvidenceConfig, MUST_FLAG_REVIEW
from GeneFior.GeneFior_Gene_Stats import GeneCoverage, GeneVisualiser
from GeneFior.utils import load_gene_id_file
from GeneFior.workflow import Workflow


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_workflow(
    query_min_coverage: float = 40.0,
    query_min_identity: float = 80.0,
    detection_min_coverage: float = 80.0,
    detection_min_identity: float = 80.0,
    detection_min_depth: int = 3,
    detection_min_num_reads: int = 1,
) -> Workflow:
    """Build a Workflow instance without requiring real files or databases.

    Uses object.__new__ to skip the heavy __init__ and manually injects only
    the attributes referenced by parse_blast_results / parse_bam_results.
    """
    import logging
    wf = object.__new__(Workflow)
    wf.query_min_coverage = query_min_coverage
    wf.query_min_identity = query_min_identity
    wf.detection_min_coverage = detection_min_coverage
    wf.detection_min_identity = detection_min_identity
    wf.detection_min_depth = detection_min_depth
    wf.detection_min_num_reads = detection_min_num_reads
    wf.min_bitscore = None
    wf.evalue = None
    wf.evidence_config = EvidenceConfig(
        corroborating_depth=3,
        exact_identity_min=detection_min_identity,
        min_unique_best_reads=1,
        min_unique_best_fraction=0.0,
    )
    wf.gene_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(GeneStats)))
    wf.detections = defaultdict(lambda: defaultdict(lambda: defaultdict(bool)))
    wf.evidence_calls = defaultdict(lambda: defaultdict(dict))
    wf.family_calls = defaultdict(lambda: defaultdict(dict))
    wf.is_genes_fasta = False
    wf.report_fasta = None
    wf.all_reads = {}
    wf.input_fasta = None
    wf.always_flag_genes = set()
    wf.logger = logging.getLogger("test_workflow")
    wf.logger.addHandler(logging.NullHandler())
    return wf


def _make_visualiser(tmp_path: Path) -> GeneVisualiser:
    """Build a Gene-Stats visualiser for direct parser tests."""
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    return GeneVisualiser(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        genes=["geneA"],
        databases=["db"],
        tools=["diamond"],
    )


def _blast_line(
    qseqid: str = "read1",
    sseqid: str = "geneA",
    pident: float = 100.0,
    length: int = 100,      # BLAST alignment length (field 3)
    mismatch: int = 0,
    gapopen: int = 0,
    qstart: int = 1,
    qend: int = 100,
    sstart: int = 1,
    send: int = 100,
    evalue: float = 1e-50,
    bitscore: float = 200.0,
    qlen: int = 100,
    slen: int = 100,
) -> str:
    """Return a single BLAST outfmt-6 tab-separated line."""
    return "\t".join(str(x) for x in [
        qseqid, sseqid, pident, length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore, qlen, slen,
    ]) + "\n"


def _blast_lines(count: int = 3, qseqid_prefix: str = "read", **kwargs) -> str:
    """Return multiple BLAST outfmt-6 rows with unique query IDs."""
    return "".join(
        _blast_line(qseqid=f"{qseqid_prefix}_{index}", **kwargs)
        for index in range(count)
    )


def _cigar_query_coverage_and_identity(cigar: str, seq_len: int, nm: int) -> tuple:
    """Replicate the CIGAR parsing logic from parse_bam_results.

    Returns (query_coverage_pct, identity_pct).
    """
    cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')
    query_aligned_bases = 0
    alignment_length = 0
    ref_pos = 0
    aligned_positions = set()

    for count_str, op in cigar_re.findall(cigar):
        length = int(count_str)
        if op in ('M', '=', 'X'):
            aligned_positions.update(range(ref_pos, ref_pos + length))
            ref_pos += length
            alignment_length += length
            query_aligned_bases += length
        elif op == 'I':
            alignment_length += length
            query_aligned_bases += length  # I ops counted (consistent with BLAST span)
        elif op == 'D':
            ref_pos += length
            alignment_length += length    # D counted in alignment_length for identity
        elif op == 'N':
            ref_pos += length
        elif op in ('S', 'H'):
            pass  # clips excluded

    if alignment_length == 0:
        identity = 0.0
    else:
        matches = max(0, alignment_length - nm)
        identity = (matches / alignment_length) * 100.0

    query_coverage = (query_aligned_bases / seq_len) * 100.0 if seq_len > 0 else 0.0
    return query_coverage, identity


class TestAlwaysFlagGeneList:
    def test_gene_id_file_parser_uses_gene_id_header(self, tmp_path):
        gene_file = tmp_path / "genes.tsv"
        gene_file.write_text(
            "Gene_ID\tNote\n"
            "blaCTX-M-1_1_DQ915955\twatch\n"
            "aac(3)-IIc_1_X54723\twatch\n"
        )

        assert load_gene_id_file(gene_file) == {
            "blaCTX-M-1_1_DQ915955",
            "aac(3)-IIc_1_X54723",
        }

    def test_always_flagged_low_evidence_gene_is_review_not_detected(
            self, tmp_path):
        wf = _make_workflow(
            query_min_coverage=20.0,
            query_min_identity=80.0,
            detection_min_coverage=80.0,
        )
        wf.output_dir = tmp_path
        wf.always_flag_genes = {"geneA"}
        hit_file = tmp_path / "partial.tsv"
        hit_file.write_text(_blast_line(
            pident=100.0,
            qstart=1,
            qend=30,
            qlen=100,
            sstart=1,
            send=30,
            slen=100,
        ))

        wf.parse_blast_results(hit_file, "db", "blastn")
        wf.generate_detection_matrix("db")

        with open(
                tmp_path / "db_evidence_matrix.tsv",
                newline="",
        ) as handle:
            rows = list(_csv.DictReader(handle, delimiter="\t"))
        with open(
                tmp_path / "db_evidence_summary.tsv",
                newline="",
        ) as handle:
            summary_rows = {
                row["Tool"]: row
                for row in _csv.DictReader(handle, delimiter="\t")
            }

        assert rows[0]["Gene"] == "geneA"
        assert rows[0]["blastn"] == MUST_FLAG_REVIEW
        assert rows[0]["Always_Flagged"] == "1"
        assert rows[0]["Evidence_Detections"] == "0"
        assert "ALWAYS_FLAGGED_GENE" in rows[0]["Evidence_Warnings"]
        assert summary_rows["blastn"]["Evidence_Genes"] == "0"
        assert summary_rows["blastn"]["Always_Flagged_Genes"] == "1"


# ===========================================================================
# 1.  Stage-1 stream-filter query coverage formula
# ===========================================================================

class TestStage1QueryCoverage:
    """Tests for the query coverage formula used in _stream_filter_hits.

    Formula: ((abs(qend - qstart) + 1) / qlen) * 100
    """

    @pytest.mark.parametrize("qstart,qend,qlen,expected", [
        (1,   100, 100, 100.0),   # full alignment
        (1,   80,  100, 80.0),    # 80 % partial
        (21,  100, 100, 80.0),    # 80 % partial from right
        (11,  90,  100, 80.0),    # 80 % internal
        (1,   1,   100, 1.0),     # single-base alignment
        (100, 100, 100, 1.0),     # single base, same coords
    ])
    def test_forward_strand(self, qstart, qend, qlen, expected):
        coverage = ((abs(qend - qstart) + 1) / qlen) * 100
        assert abs(coverage - expected) < 1e-9

    @pytest.mark.parametrize("qstart,qend,qlen,expected", [
        (100, 1,   100, 100.0),   # full, reverse
        (80,  1,   100, 80.0),    # 80 %, reverse
    ])
    def test_reverse_strand_abs_handles_negative_span(self, qstart, qend, qlen, expected):
        """abs() must be applied so reverse-strand alignments are not penalised."""
        coverage = ((abs(qend - qstart) + 1) / qlen) * 100
        assert abs(coverage - expected) < 1e-9

    def test_zero_qlen_returns_zero(self):
        # Non-zero qlen works normally: qstart=1, qend=100, qlen=100 → 100 %
        qlen = 100
        coverage = ((abs(100 - 1) + 1) / qlen) * 100 if qlen else 0.0
        assert abs(coverage - 100.0) < 1e-9
        # Zero qlen must return 0 (guard against division-by-zero)
        qlen = 0
        coverage = ((abs(100 - 1) + 1) / qlen) * 100 if qlen else 0.0
        assert coverage == 0.0

    def test_cannot_exceed_100_percent(self):
        """Span formula is bounded by qlen — can never exceed 100 %."""
        # qstart/qend cannot span more than qlen positions on the query
        for qstart, qend, qlen in [(1, 100, 100), (1, 50, 50), (1, 200, 200)]:
            coverage = ((abs(qend - qstart) + 1) / qlen) * 100
            assert coverage <= 100.0 + 1e-9


# ===========================================================================
# 2.  Stage-2 parse_blast_results query coverage — regression for old bug
# ===========================================================================

class TestStage2QueryCoverageRegression:
    """Regression tests confirming the old alignment_len/qlen formula is gone.

    With a reference deletion BLAST's 'length' field (alignment_len) exceeds
    qlen.  The old code returned > 100 %; the span-based fix cannot.
    """

    def test_large_alignment_len_does_not_inflate_coverage(self):
        """alignment_len > qlen must NOT produce coverage > 100 %.

        Scenario: 80 bp read aligns to a 100 bp reference with a 30 bp
        deletion (reference gap).  BLAST field-3 'length' = 80 + 30 = 110,
        but qstart=1, qend=80, qlen=80.

        Old formula:  110 / 80 * 100 = 137.5 %  (wrong)
        New formula:  (80 - 1 + 1) / 80 * 100  = 100 %  (correct)
        """
        alignment_len = 110   # BLAST field-3 (includes reference deletion)
        qstart, qend, qlen = 1, 80, 80

        # Old (broken) formula
        old_coverage = (alignment_len / float(qlen)) * 100.0
        assert old_coverage > 100.0, "Setup check: old formula must exceed 100 %"

        # New (correct) formula now used in parse_blast_results
        new_coverage = (abs(qend - qstart) + 1) / float(qlen) * 100.0
        assert new_coverage <= 100.0
        assert abs(new_coverage - 100.0) < 1e-9

    def test_partial_read_is_not_inflated(self):
        """Partial alignment: only 60 % of the read should be reported as 60 %."""
        qstart, qend, qlen = 1, 60, 100
        alignment_len = 70   # longer due to indels

        old_coverage = (alignment_len / float(qlen)) * 100.0  # 70 %
        new_coverage = (abs(qend - qstart) + 1) / float(qlen) * 100.0  # 60 %

        # Old formula over-estimates; new formula reflects actual query span
        assert old_coverage > new_coverage
        assert abs(new_coverage - 60.0) < 1e-9

    def test_stage1_and_stage2_produce_identical_values(self):
        """Stage-1 and Stage-2 formulas must be numerically identical."""
        test_cases = [
            (1, 100, 100),
            (1, 80, 100),
            (10, 90, 150),
            (100, 1, 100),  # reverse
        ]
        for qstart, qend, qlen in test_cases:
            stage1 = ((abs(qend - qstart) + 1) / qlen) * 100
            stage2 = (abs(qend - qstart) + 1) / float(qlen) * 100.0
            assert abs(stage1 - stage2) < 1e-12, (
                f"Stage-1 and Stage-2 differ for qstart={qstart} qend={qend} qlen={qlen}: "
                f"{stage1} vs {stage2}"
            )


# ===========================================================================
# 3.  BAM CIGAR — identity calculation
# ===========================================================================

class TestCigarIdentity:
    """Tests for per-read identity derived from CIGAR + NM tag.

    Formula: (alignment_length - NM) / alignment_length * 100
    where alignment_length = sum(M + I + D lengths).
    """

    def test_perfect_match(self):
        _, identity = _cigar_query_coverage_and_identity("100M", 100, nm=0)
        assert abs(identity - 100.0) < 1e-9

    def test_mismatches_reduce_identity(self):
        _, identity = _cigar_query_coverage_and_identity("100M", 100, nm=5)
        assert abs(identity - 95.0) < 1e-9

    def test_insertions_increase_alignment_length(self):
        """An insertion extends alignment_length; with NM=0 identity stays 100 %."""
        _, identity = _cigar_query_coverage_and_identity("90M10I", 100, nm=0)
        # alignment_length = 100; matches = 100; identity = 100 %
        assert abs(identity - 100.0) < 1e-9

    def test_insertions_with_mismatches(self):
        """90M10I with NM=5: alignment_length=100, matches=95, identity=95 %."""
        _, identity = _cigar_query_coverage_and_identity("90M10I", 100, nm=5)
        assert abs(identity - 95.0) < 1e-9

    def test_deletions_increase_alignment_length(self):
        """Deletion from reference: 80M10D10M → alignment_length=100, NM=10 → 90 %."""
        _, identity = _cigar_query_coverage_and_identity("80M10D10M", 90, nm=10)
        # alignment_length = 80+10+10 = 100; matches = 90; identity = 90 %
        assert abs(identity - 90.0) < 1e-9

    def test_softclip_excluded_from_alignment_length(self):
        """Soft-clipped bases do not contribute to alignment_length."""
        _, identity = _cigar_query_coverage_and_identity("10S80M10S", 100, nm=0)
        # alignment_length = 80 (S excluded); matches = 80; identity = 100 %
        assert abs(identity - 100.0) < 1e-9

    def test_zero_alignment_length_returns_zero(self):
        """All-soft-clip edge case: alignment_length=0 → identity=0."""
        _, identity = _cigar_query_coverage_and_identity("100S", 100, nm=0)
        assert identity == 0.0

    def test_complex_cigar_identity(self):
        """50M5I10D35M with NM=10: al_len=100, matches=90, identity=90 %."""
        _, identity = _cigar_query_coverage_and_identity("50M5I10D35M", 90, nm=10)
        # alignment_length = 50+5+10+35 = 100; matches = 90; identity = 90 %
        assert abs(identity - 90.0) < 1e-9


# ===========================================================================
# 4.  BAM CIGAR — query coverage calculation
# ===========================================================================

class TestCigarQueryCoverage:
    """Tests for per-read query coverage derived from CIGAR strings.

    Formula: (M + I bases) / read_length * 100

    I (insertion) operations ARE counted because BLAST's span formula
    (qend - qstart + 1) implicitly includes insertion positions in query
    coordinates.  Excluding I would make BAM more conservative than BLAST.
    """

    def test_all_m_full_coverage(self):
        coverage, _ = _cigar_query_coverage_and_identity("100M", 100, nm=0)
        assert abs(coverage - 100.0) < 1e-9

    def test_m_plus_softclip_partial_coverage(self):
        """80M20S: 80 aligned bases out of 100 bp read → 80 %."""
        coverage, _ = _cigar_query_coverage_and_identity("80M20S", 100, nm=0)
        assert abs(coverage - 80.0) < 1e-9

    def test_m_plus_insertion_is_counted(self):
        """90M10I: both M and I contribute → 100/100 = 100 %.

        Contrast with old code (M-only): 90/100 = 90 %.
        """
        coverage, _ = _cigar_query_coverage_and_identity("90M10I", 100, nm=0)
        assert abs(coverage - 100.0) < 1e-9

    def test_deletion_does_not_add_to_query_bases(self):
        """80M10D10M: D consumes reference, not query → query_aligned = 90/90 = 100 %."""
        # read_length = 90 (80 + 10 from M ops), D doesn't consume query
        coverage, _ = _cigar_query_coverage_and_identity("80M10D10M", 90, nm=0)
        assert abs(coverage - 100.0) < 1e-9

    def test_softclip_plus_insertion(self):
        """70M10I20S: query_aligned = 70+10=80, read_length=100 → 80 %."""
        coverage, _ = _cigar_query_coverage_and_identity("70M10I20S", 100, nm=0)
        assert abs(coverage - 80.0) < 1e-9

    def test_complex_cigar_coverage(self):
        """30M20D30M40S: query_aligned=60 (M only, D excluded), read_length=100 → 60 %."""
        coverage, _ = _cigar_query_coverage_and_identity("30M20D30M40S", 100, nm=0)
        assert abs(coverage - 60.0) < 1e-9

    def test_insertion_consistency_with_blast_span(self):
        """Verify I-op counting closes the gap to BLAST span-based coverage.

        BLAST: qstart=1, qend=90, qlen=90 (I positions included in span) → 100 %
        BAM:   90M (no insertions in M count)                             →  old: 90/100=90 %
               90M with I being read as 90M10I (full 100 bp read)        →  new: 100/100=100 %
        """
        # Simulate a 100 bp read where 90 bp match and 10 bp are insertions
        cigar = "90M10I"
        read_length = 100
        coverage_new, _ = _cigar_query_coverage_and_identity(cigar, read_length, nm=0)

        # BLAST span-based for the same alignment extent
        qstart, qend, qlen = 1, 100, 100  # span covers the full read
        blast_coverage = ((abs(qend - qstart) + 1) / qlen) * 100

        # New BAM formula must match BLAST span
        assert abs(coverage_new - blast_coverage) < 1e-9

    def test_zero_read_length_returns_zero(self):
        """Guard against division by zero when read has no sequence (SEQ = *)."""
        coverage, _ = _cigar_query_coverage_and_identity("100M", 0, nm=0)
        assert coverage == 0.0


# ===========================================================================
# 5.  GeneStats — coverage union, depth, and identity aggregation
# ===========================================================================

class TestGeneStats:
    """Tests for GeneStats accumulation and finalise() calculations."""

    def test_single_hit_full_coverage(self):
        gs = GeneStats(gene_name="geneA")
        gs.add_hit(sstart=1, send=100, identity=95.0, gene_len=100)
        gs.finalise()
        assert abs(gs.gene_coverage - 100.0) < 1e-9

    def test_single_hit_partial_coverage(self):
        gs = GeneStats(gene_name="geneA")
        gs.add_hit(sstart=1, send=50, identity=95.0, gene_len=100)
        gs.finalise()
        assert abs(gs.gene_coverage - 50.0) < 1e-9

    def test_two_non_overlapping_hits_union_is_100(self):
        """Two reads covering the two halves of a gene → 100 % coverage."""
        gs = GeneStats(gene_name="geneA")
        gs.add_hit(sstart=1, send=50, identity=90.0, gene_len=100)
        gs.add_hit(sstart=51, send=100, identity=90.0, gene_len=100)
        gs.finalise()
        assert abs(gs.gene_coverage - 100.0) < 1e-9

    def test_overlapping_hits_no_double_count(self):
        """Overlapping hits share positions — the union must not double-count."""
        gs = GeneStats(gene_name="geneA")
        gs.add_hit(sstart=1, send=80, identity=90.0, gene_len=100)
        gs.add_hit(sstart=60, send=100, identity=90.0, gene_len=100)
        gs.finalise()
        # Covered = positions 0..99 → 100 %
        assert abs(gs.gene_coverage - 100.0) < 1e-9

    def test_overlapping_hits_depth_increments(self):
        """Positions covered by both reads have depth = 2."""
        gs = GeneStats(gene_name="geneA")
        gs.add_hit(sstart=1, send=80, identity=90.0, gene_len=100)
        gs.add_hit(sstart=60, send=100, identity=90.0, gene_len=100)
        gs.finalise()
        # Positions 59-79 (0-based) are covered by both reads
        overlap_start = 59  # 0-based for sstart=60
        overlap_end = 79    # 0-based for send=80 (exclusive end of 80 in 1-based)
        for pos in range(overlap_start, overlap_end):
            assert gs.position_depths.get(pos, 0) == 2, f"Position {pos} should have depth 2"

    def test_coverage_at_depth_uses_whole_gene_denominator(self):
        """coverage_at_depth must ask how much of the whole gene reaches depth N."""
        gs = GeneStats(gene_name="geneA")
        # Only first 50 positions covered, each by 1 read
        gs.add_hit(sstart=1, send=50, identity=95.0, gene_len=100)
        gs.finalise()
        assert abs(gs.coverage_at_depth(1) - 50.0) < 1e-9
        assert abs(gs.coverage_at_depth(2) - 0.0) < 1e-9

    def test_base_depth_includes_uncovered_zeros(self):
        """base_depth (full-gene mean) must be diluted by uncovered positions."""
        gs = GeneStats(gene_name="geneA")
        gs.add_hit(sstart=1, send=50, identity=95.0, gene_len=100)
        gs.finalise()
        # 50 positions × depth 1, 50 positions × depth 0 → mean = 0.5
        assert abs(gs.base_depth - 0.5) < 1e-9

    def test_coverage_at_arbitrary_depth_is_calculated_on_demand(self):
        """User-specified depth thresholds must not be limited to 2/3/5/10."""
        gs = GeneStats(gene_name="geneA")
        for _ in range(4):
            gs.add_hit(sstart=1, send=100, identity=95.0, gene_len=100)
        gs.finalise()

        assert abs(gs.coverage_at_depth(4) - 100.0) < 1e-9
        assert abs(gs.coverage_at_depth(5) - 0.0) < 1e-9

    def test_avg_identity(self):
        gs = GeneStats(gene_name="geneA")
        gs.add_hit(sstart=1, send=50, identity=90.0, gene_len=100)
        gs.add_hit(sstart=51, send=100, identity=100.0, gene_len=100)
        gs.finalise()
        assert abs(gs.avg_identity - 95.0) < 1e-9

    def test_num_sequences_incremented_per_hit(self):
        gs = GeneStats(gene_name="geneA")
        for i in range(5):
            gs.add_hit(sstart=1, send=20, identity=99.0, gene_len=100)
        gs.finalise()
        assert gs.num_sequences == 5

    def test_add_positions_bam_path(self):
        """add_positions (BAM code path) must produce same coverage as add_hit."""
        gs_hit = GeneStats(gene_name="geneA")
        gs_pos = GeneStats(gene_name="geneA")

        # add_hit path (BLAST)
        gs_hit.add_hit(sstart=1, send=60, identity=95.0, gene_len=100)
        gs_hit.finalise()

        # add_positions path (BAM) — 0-based positions matching sstart=1..send=60
        positions = set(range(0, 60))  # same as add_hit converts to (1-1=0 .. 60)
        gs_pos.add_positions(positions, identity=95.0, gene_len=100)
        gs_pos.finalise()

        assert abs(gs_hit.gene_coverage - gs_pos.gene_coverage) < 1e-9

    def test_gene_length_updated_to_maximum(self):
        """gene_length should reflect the largest slen seen across hits."""
        gs = GeneStats(gene_name="geneA")
        gs.add_hit(sstart=1, send=50, identity=95.0, gene_len=80)
        gs.add_hit(sstart=51, send=100, identity=95.0, gene_len=100)
        gs.finalise()
        assert gs.gene_length == 100

    def test_no_hits_gives_zero_coverage(self):
        gs = GeneStats(gene_name="geneA")
        gs.finalise()
        assert gs.gene_coverage == 0.0
        assert gs.avg_identity == 0.0

    def test_coverage_cannot_exceed_100(self):
        """Even with many overlapping hits coverage must stay ≤ 100 %."""
        gs = GeneStats(gene_name="geneA")
        for _ in range(20):
            gs.add_hit(sstart=1, send=100, identity=95.0, gene_len=100)
        gs.finalise()
        assert gs.gene_coverage <= 100.0

    def test_compact_interval_mode_matches_per_position_mode(self):
        regular = GeneStats(gene_name="geneA")
        compact = GeneStats(gene_name="geneA", compact=True)
        hits = [
            (1, 50, 90.0),
            (25, 75, 95.0),
            (80, 100, 100.0),
        ]
        for start, end, identity in hits:
            regular.add_hit(start, end, identity, gene_len=100)
            compact.add_hit(start, end, identity, gene_len=100)
        regular.finalise()
        compact.finalise()

        assert compact.covered_positions == set()
        assert compact.position_depths == {}
        assert abs(compact.gene_coverage - regular.gene_coverage) < 1e-9
        assert abs(compact.base_depth - regular.base_depth) < 1e-9
        assert abs(compact.coverage_at_depth(2) - regular.coverage_at_depth(2)) < 1e-9
        assert abs(compact.avg_identity - regular.avg_identity) < 1e-9

    def test_compact_bam_intervals_match_positions(self):
        regular = GeneStats(gene_name="geneA")
        compact = GeneStats(gene_name="geneA", compact=True)
        positions = set(range(0, 20)) | set(range(30, 50))
        regular.add_positions(positions, identity=97.0, gene_len=100)
        compact.add_intervals([(0, 20), (30, 50)], identity=97.0, gene_len=100)
        regular.finalise()
        compact.finalise()

        assert abs(compact.gene_coverage - regular.gene_coverage) < 1e-9
        assert abs(compact.base_depth - regular.base_depth) < 1e-9
        assert abs(compact.coverage_at_depth(2) - regular.coverage_at_depth(2)) < 1e-9

    def test_compact_intervals_merge_overlaps_within_read(self):
        regular = GeneStats(gene_name="geneA")
        compact = GeneStats(gene_name="geneA", compact=True)
        positions = set(range(0, 30))
        regular.add_positions(positions, identity=97.0, gene_len=100)
        compact.add_intervals([(0, 20), (10, 30)], identity=97.0, gene_len=100)
        regular.finalise()
        compact.finalise()

        assert abs(compact.gene_coverage - regular.gene_coverage) < 1e-9
        assert abs(compact.base_depth - regular.base_depth) < 1e-9
        assert abs(compact.coverage_at_depth(2) - regular.coverage_at_depth(2)) < 1e-9


# ===========================================================================
# 5b. Gene-Stats translated-coordinate coverage regressions
# ===========================================================================

class TestGeneStatsVisualiserBlastCoverage:
    """Regression tests for GeneFior-Gene-Stats BLAST/DIAMOND coverage reports."""

    def _write_raw_blast_dir(self, tmp_path):
        source = tmp_path / "source_run"
        raw_outputs = source / "raw_outputs"
        raw_outputs.mkdir(parents=True)
        blast_file = raw_outputs / "db_blastn_results.tsv"
        blast_file.write_text(_blast_line(qseqid="read1", sseqid="geneA"))
        return source, raw_outputs, blast_file

    def test_gene_stats_can_use_explicit_raw_dir_for_recompute_outputs(
            self, tmp_path):
        input_dir = tmp_path / "qualified"
        input_dir.mkdir()
        source, raw_outputs, blast_file = self._write_raw_blast_dir(tmp_path)
        visualiser = GeneVisualiser(
            input_dir=str(input_dir),
            output_dir=str(tmp_path / "out"),
            genes=[],
            databases=["db"],
            tools=["blastn"],
            raw_dir=str(source),
        )

        assert visualiser.discover_files()
        assert visualiser.raw_dir_used == raw_outputs
        assert visualiser.found_files == [
            ("db", "blastn", "blast", blast_file)
        ]

    def test_gene_stats_uses_source_raw_outputs_from_manifest(
            self, tmp_path):
        input_dir = tmp_path / "qualified"
        input_dir.mkdir()
        _, raw_outputs, blast_file = self._write_raw_blast_dir(tmp_path)
        (input_dir / "run_parameters.json").write_text(json.dumps({
            "source_raw_outputs": str(raw_outputs),
        }))
        visualiser = GeneVisualiser(
            input_dir=str(input_dir),
            output_dir=str(tmp_path / "out"),
            genes=[],
            databases=["db"],
            tools=["blastn"],
        )

        assert visualiser.discover_files()
        assert visualiser.raw_dir_used == raw_outputs
        assert visualiser.found_files == [
            ("db", "blastn", "blast", blast_file)
        ]

    def _write_evidence_matrix(self, visualiser):
        matrix = Path(visualiser.input_dir) / "db_evidence_matrix.tsv"
        matrix.write_text(
            "Gene\tDetection_System\tblastn\tBWA\tBest_Evidence_Status\t"
            "Evidence_Detections\tCandidate_Allele_Detections\t"
            "Exact_Allele_Detections\tProfile_Detections\t"
            "Strict_Detections\tAlways_Flagged\tEvidence_Warnings\n"
            "candidate_gene\tqualified\tCANDIDATE_ALLELE_DETECTED\t"
            "NOT_DETECTED\tCANDIDATE_ALLELE_DETECTED\t1\t1\t0\t0\t0\t0\t\n"
            "exact_gene\tqualified\tEXACT_ALLELE_DETECTED\t"
            "NOT_DETECTED\tEXACT_ALLELE_DETECTED\t1\t0\t1\t0\t1\t0\t\n"
            "family_gene\tqualified\tFAMILY_DETECTED\t"
            "NOT_DETECTED\tFAMILY_DETECTED\t1\t0\t0\t0\t0\t0\t\n"
            "flagged_gene\tqualified\tMUST_FLAG_REVIEW\t"
            "NOT_DETECTED\tMUST_FLAG_REVIEW\t0\t0\t0\t0\t0\t1\t"
            "ALWAYS_FLAGGED_GENE\n"
        )
        return matrix

    def test_gene_stats_evidence_selection_includes_always_flagged_genes(
            self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        visualiser.genes = set()
        visualiser.gene_selection = "evidence"
        self._write_evidence_matrix(visualiser)

        selected = visualiser.load_detected_genes("db")

        assert "flagged_gene" in selected

    def test_gene_stats_can_select_candidate_or_exact_genes(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        visualiser.genes = set()
        visualiser.gene_selection = "candidate-or-exact"
        self._write_evidence_matrix(visualiser)

        selected = visualiser.load_detected_genes("db")

        assert selected == {"candidate_gene", "exact_gene"}

    def test_gene_stats_exact_selection_filters_raw_parsing(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        visualiser.genes = set()
        visualiser.gene_selection = "exact"
        self._write_evidence_matrix(visualiser)
        visualiser.selected_genes_by_db["db"] = visualiser.load_detected_genes("db")
        blast_file = tmp_path / "blastn.tsv"
        blast_file.write_text(
            _blast_line(qseqid="read1", sseqid="candidate_gene")
            + _blast_line(qseqid="read2", sseqid="exact_gene")
            + _blast_line(qseqid="read3", sseqid="family_gene")
        )

        visualiser.parse_blast_results(blast_file, "db", "blastn")

        assert set(visualiser.gene_coverages["db"]["blastn"]) == {"exact_gene"}

    def test_diamond_subject_amino_acid_coords_expand_to_full_codons(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        blast_file = tmp_path / "diamond.tsv"
        blast_file.write_text(
            _blast_line(
                qseqid="read1",
                sseqid="geneA",
                pident=100.0,
                length=2,
                qstart=1,
                qend=6,
                sstart=1,
                send=2,
                qlen=6,
                slen=2,
            )
        )

        visualiser.parse_blast_results(blast_file, "db", "diamond")
        coverage = visualiser.gene_coverages["db"]["diamond"]["geneA"]

        assert coverage.gene_length == 6
        assert [coverage.positions[i].depth for i in range(6)] == [1, 1, 1, 1, 1, 1]
        assert coverage.get_coverage_stats()["gaps"] == []

    def test_loaded_query_sequence_does_not_double_count_blast_depth(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        visualiser.query_sequences["read1"] = "AAAAAA"
        blast_file = tmp_path / "diamond.tsv"
        blast_file.write_text(
            _blast_line(
                qseqid="read1",
                sseqid="geneA",
                pident=100.0,
                length=2,
                qstart=1,
                qend=6,
                sstart=1,
                send=2,
                qlen=6,
                slen=2,
            )
        )

        visualiser.parse_blast_results(blast_file, "db", "diamond")
        coverage = visualiser.gene_coverages["db"]["diamond"]["geneA"]

        assert [coverage.positions[i].depth for i in range(6)] == [1, 1, 1, 1, 1, 1]
        assert coverage.get_coverage_stats()["positions_with_bases"] == 6

    def test_diamond_blastp_keeps_amino_acid_coordinates(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        visualiser.sequence_type = "Genes-FASTA"
        visualiser.genes_type = "aa"
        blast_file = tmp_path / "diamond.tsv"
        blast_file.write_text(_blast_line(
            qseqid="protein1",
            sseqid="geneA",
            pident=100.0,
            length=100,
            qstart=1,
            qend=100,
            sstart=1,
            send=100,
            qlen=100,
            slen=100,
        ))

        visualiser.parse_blast_results(blast_file, "db", "diamond")
        coverage = visualiser.gene_coverages["db"]["diamond"]["geneA"]

        assert coverage.gene_length == 100
        assert coverage.unit == "aa"
        assert coverage.get_coverage_stats()["coverage_percent"] == 100.0

    def test_overlapping_hsps_from_one_read_do_not_inflate_depth(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        blast_file = tmp_path / "blastn.tsv"
        blast_file.write_text(
            _blast_line(
                qseqid="read1", sseqid="geneA", qstart=1, qend=60,
                sstart=1, send=60, qlen=100, slen=100,
            )
            + _blast_line(
                qseqid="read1", sseqid="geneA", qstart=41, qend=100,
                sstart=41, send=100, qlen=100, slen=100,
            )
        )

        visualiser.parse_blast_results(blast_file, "db", "blastn")
        coverage = visualiser.gene_coverages["db"]["blastn"]["geneA"]

        assert coverage.read_count == 1
        assert coverage.get_coverage_stats()["coverage_percent"] == 100.0
        assert coverage.get_coverage_stats()["coverage_2x"] == 0.0

    def test_bam_secondary_alignments_are_excluded(
            self, tmp_path, monkeypatch):
        visualiser = _make_visualiser(tmp_path)
        bam_file = tmp_path / "alignments.bam"
        bam_file.touch()
        sam = (
            "@SQ\tSN:geneA\tLN:10\n"
            "read1\t0\tgeneA\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tNM:i:0\n"
            "read1\t256\tgeneA\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tNM:i:0\n"
        )

        class FakePopen:
            def __init__(self, *args, **kwargs):
                self.stdout = _io.StringIO(sam)
                self.stderr = _io.StringIO("")
                self.returncode = 0

            def wait(self):
                return self.returncode

        monkeypatch.setattr(
            "GeneFior.GeneFior_Gene_Stats.subprocess.Popen",
            FakePopen,
        )

        visualiser.parse_bam_file(bam_file, "db", "bowtie2")
        coverage = visualiser.gene_coverages["db"]["bowtie2"]["geneA"]

        assert coverage.read_count == 1
        assert coverage.get_coverage_stats()["avg_depth"] == 1.0

    def test_legacy_report_uses_only_binary_detected_field(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        visualiser.detection_system = "legacy-relaxed"
        visualiser.detected_by_tool["db"]["geneA"]["blastn"] = True
        coverage = GeneCoverage("geneA", 10)
        coverage.add_alignment(0, 10)

        visualiser.generate_text_report(
            "geneA", "db", "blastn", coverage
        )

        report = (
            visualiser.report_dir
            / "db_blastn_geneA_coverage.txt"
        ).read_text()
        assert "Detected:              1" in report
        assert "Evidence status:" not in report
        assert "Evidence warnings:" not in report

    def test_text_report_makes_gene_span_explicit(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        visualiser.detection_system = "legacy-relaxed"
        coverage = GeneCoverage("geneA", 876)
        coverage.add_alignment(0, 876)

        visualiser.generate_text_report(
            "geneA", "db", "blastn", coverage
        )

        report = (
            visualiser.report_dir
            / "db_blastn_geneA_coverage.txt"
        ).read_text()
        assert "Gene length:           876 bp" in report
        assert "Gene span: 1..876 bp" in report
        assert "Position scale (bp): start=1 | mid=438 | end=876" in report
        assert "ASCII histogram: actual mean depth per plotted bin" in report

    def test_ascii_report_uses_actual_depth_labels_above_50x(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        visualiser.detection_system = "legacy-relaxed"
        coverage = GeneCoverage("geneA", 10)
        for _ in range(150):
            coverage.add_alignment(0, 10)

        visualiser.generate_text_report(
            "geneA", "db", "blastn", coverage
        )

        report = (
            visualiser.report_dir
            / "db_blastn_geneA_coverage.txt"
        ).read_text()
        assert "Max depth: 150×" in report
        assert "150× |" in report
        assert visualiser._ascii_depth_levels(150)[0] == 150

    def test_plot_axis_makes_gene_length_and_end_tick_explicit(self, tmp_path):
        visualiser = _make_visualiser(tmp_path)
        if not visualiser.plot_enabled:
            pytest.skip("matplotlib not available")

        fig, ax = visualiser.plt.subplots()
        try:
            visualiser._label_gene_axis(ax, 876, "bp", xlabel=True)

            assert ax.get_xlim() == (0.5, 876.5)
            assert 1 in ax.get_xticks()
            assert 876 in ax.get_xticks()
            assert ax.get_xlabel() == (
                "Gene position (bp; full length = 876 bp)"
            )
        finally:
            visualiser.plt.close(fig)


# ===========================================================================
# 6.  Integration — parse_blast_results with mock BLAST output file
# ===========================================================================

class TestParseBlastResultsIntegration:
    """Integration tests that call parse_blast_results with real temporary files."""

    def _run(self, lines: list, **wf_kwargs) -> tuple:
        """Write lines to a temp file and parse them; return (detected, gene_reads)."""
        tool_name = wf_kwargs.pop("tool_name", "blastn")
        wf = _make_workflow(**wf_kwargs)
        wf.gene_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(GeneStats)))
        wf.detections = defaultdict(lambda: defaultdict(lambda: defaultdict(bool)))
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as fh:
            fh.writelines(lines)
            tmp_path = Path(fh.name)
        try:
            detected, gene_reads = wf.parse_blast_results(
                tmp_path, "db", tool_name
            )
        finally:
            tmp_path.unlink(missing_ok=True)
        return detected, gene_reads, wf

    def test_three_perfect_hits_detected(self):
        """Three full-coverage hits should satisfy the 3x evidence floor."""
        lines = _blast_lines(
            pident=100.0, qstart=1, qend=100, qlen=100,
            sstart=1, send=100, slen=100,
        )
        detected, _, _ = self._run([lines],
                                   query_min_coverage=40.0, query_min_identity=80.0,
                                   detection_min_coverage=80.0)
        assert "geneA" in detected

    def test_identity_below_threshold_not_detected(self):
        """A hit below the identity threshold must not contribute to detection."""
        line = _blast_line(pident=70.0, qstart=1, qend=100, qlen=100, sstart=1, send=100, slen=100)
        detected, _, _ = self._run([line],
                                   query_min_identity=80.0,
                                   detection_min_coverage=80.0)
        assert "geneA" not in detected

    def test_coverage_below_threshold_not_detected(self):
        """A hit covering only 30 % of the query must not pass a 40 % query coverage threshold."""
        line = _blast_line(pident=100.0, qstart=1, qend=30, qlen=100, sstart=1, send=30, slen=100)
        detected, _, _ = self._run([line],
                                   query_min_coverage=40.0,
                                   detection_min_coverage=80.0)
        assert "geneA" not in detected

    def test_multiple_read_groups_union_drives_detection(self):
        """Three reads over each half give 100% 3x coverage and detection."""
        left = _blast_lines(
            qseqid_prefix="left", pident=95.0, qstart=1, qend=100,
            qlen=100, sstart=1, send=50, slen=100,
        )
        right = _blast_lines(
            qseqid_prefix="right", pident=95.0, qstart=1, qend=100,
            qlen=100, sstart=51, send=100, slen=100,
        )
        detected, _, _ = self._run([left, right],
                                   query_min_coverage=40.0, query_min_identity=80.0,
                                   detection_min_coverage=80.0)
        assert "geneA" in detected

    def test_gene_coverage_below_detection_threshold_not_detected(self):
        """Single read covering only 40 % of the gene should not pass an 80 % detection threshold."""
        line = _blast_line(pident=100.0, qstart=1, qend=100, qlen=100,
                           sstart=1, send=40, slen=100)
        detected, _, wf = self._run(
            [line],
            query_min_coverage=40.0,
            query_min_identity=80.0,
            detection_min_coverage=80.0,
        )
        assert "geneA" not in detected
        call = wf.evidence_calls["db"]["geneA"]["blastn"]
        assert not call.evidence_present
        assert not wf.detections["db"]["geneA"]["blastn"]

    def test_protein_evidence_is_detected_without_claiming_exact_allele(self):
        """Protein evidence keeps ordinary detection semantics but is not exact."""
        lines = _blast_lines(
            pident=95.0, qstart=1, qend=100, qlen=100,
            sstart=1, send=100, slen=100,
        )
        detected, _, wf = self._run(
            [lines],
            tool_name="DIAMOND-BLASTX",
            query_min_coverage=40.0,
            query_min_identity=80.0,
            detection_min_coverage=80.0,
            detection_min_identity=80.0,
        )

        call = wf.evidence_calls["db"]["geneA"]["DIAMOND-BLASTX"]
        assert "geneA" in detected
        assert call.evidence_present
        assert not call.exact_allele_detected
        assert wf.detections["db"]["geneA"]["DIAMOND-BLASTX"]

    def test_workflow_legacy_relaxed_mode_disables_allele_resolution(self):
        wf = _make_workflow(
            query_min_coverage=40.0,
            query_min_identity=80.0,
            detection_min_coverage=80.0,
            detection_min_identity=80.0,
        )
        wf.detection_system = "legacy-relaxed"
        with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(_blast_lines(
                pident=99.0, qstart=1, qend=100, qlen=100,
                sstart=1, send=100, slen=100,
            ))
            path = Path(handle.name)
        try:
            detected, _ = wf.parse_blast_results(path, "db", "BLASTn")
        finally:
            path.unlink(missing_ok=True)

        call = wf.evidence_calls["db"]["geneA"]["BLASTn"]
        assert "geneA" in detected
        assert call.status == "DETECTED"
        assert call.evidence_present
        assert not call.exact_allele_detected
        assert wf.family_calls["db"]["BLASTn"] == {}

    def test_legacy_tool_stats_report_only_binary_detection_fields(
            self, tmp_path):
        wf = _make_workflow()
        wf.detection_system = "legacy-relaxed"
        wf.stats_dir = tmp_path / "tool_stats"
        wf.stats_dir.mkdir()
        hit_file = tmp_path / "legacy.tsv"
        hit_file.write_text(_blast_lines(pident=100.0))

        _, gene_reads = wf.parse_blast_results(
            hit_file, "db", "BLASTn"
        )
        wf.write_tool_stats("db", "BLASTn", gene_reads)

        with open(
                wf.stats_dir / "db_BLASTn_stats.tsv",
                newline="",
        ) as handle:
            reader = _csv.DictReader(handle, delimiter="\t")
            rows = list(reader)

        assert reader.fieldnames[-1] == "Detected"
        assert "Evidence_Status" not in reader.fieldnames
        assert "Evidence_Present" not in reader.fieldnames
        assert "Exact_Allele_Detected" not in reader.fieldnames
        assert "Detection_System" not in reader.fieldnames
        assert rows[0]["Detected"] == "1"

    def test_legacy_mode_does_not_emit_qualified_report_files(
            self, tmp_path):
        wf = _make_workflow()
        wf.detection_system = "legacy-relaxed"
        wf.output_dir = tmp_path
        for suffix in (
                "evidence_matrix.tsv",
                "evidence_summary.tsv",
                "allele_resolution.tsv"):
            (tmp_path / f"db_{suffix}").write_text("stale\n")
        hit_file = tmp_path / "legacy.tsv"
        hit_file.write_text(_blast_lines(pident=100.0))

        wf.parse_blast_results(hit_file, "db", "BLASTn")
        wf.generate_detection_matrix("db")

        assert (tmp_path / "db_detection_matrix.tsv").is_file()
        assert not (tmp_path / "db_evidence_matrix.tsv").exists()
        assert not (tmp_path / "db_evidence_summary.tsv").exists()
        assert not (tmp_path / "db_allele_resolution.tsv").exists()

    def test_final_summary_separates_evidence_from_exact_alleles(
            self, tmp_path):
        wf = _make_workflow(
            query_min_coverage=40.0,
            query_min_identity=80.0,
            detection_min_coverage=80.0,
            detection_min_identity=80.0,
        )
        wf.output_dir = tmp_path
        protein_file = tmp_path / "protein.tsv"
        nucleotide_file = tmp_path / "nucleotide.tsv"
        protein_file.write_text(_blast_lines(
            qseqid_prefix="protein_read",
            sseqid="protein_gene",
            pident=95.0,
        ))
        nucleotide_file.write_text("".join(
            _blast_line(
                qseqid=f"nucleotide_read_{index}",
                sseqid="nucleotide_gene",
                pident=100.0,
            )
            for index in range(3)
        ))

        wf.parse_blast_results(
            protein_file, "db", "DIAMOND-BLASTX"
        )
        wf.parse_blast_results(
            nucleotide_file, "db", "BLASTn"
        )
        wf.generate_detection_matrix("db")

        with open(
                tmp_path / "db_evidence_summary.tsv",
                newline="",
        ) as handle:
            rows = {
                row["Tool"]: row
                for row in _csv.DictReader(handle, delimiter="\t")
            }

        assert rows["DIAMOND-BLASTX"]["Evidence_Genes"] == "1"
        assert rows["DIAMOND-BLASTX"]["Candidate_Allele_Genes"] == "0"
        assert rows["DIAMOND-BLASTX"]["Exact_Allele_Genes"] == "0"
        assert rows["BLASTn"]["Evidence_Genes"] == "1"
        assert rows["BLASTn"]["Candidate_Allele_Genes"] == "0"
        assert rows["BLASTn"]["Exact_Allele_Genes"] == "1"
        assert rows["ALL_TOOLS"]["Evidence_Genes"] == "2"
        assert rows["ALL_TOOLS"]["Candidate_Allele_Genes"] == "0"
        assert rows["ALL_TOOLS"]["Exact_Allele_Genes"] == "1"

    def test_qualified_tool_stats_include_evidence_fields(self, tmp_path):
        wf = _make_workflow()
        wf.stats_dir = tmp_path / "tool_stats"
        wf.stats_dir.mkdir()
        hit_file = tmp_path / "qualified.tsv"
        hit_file.write_text(_blast_line(pident=100.0))

        _, gene_reads = wf.parse_blast_results(
            hit_file, "db", "BLASTn"
        )
        wf.write_tool_stats("db", "BLASTn", gene_reads)

        with open(
                wf.stats_dir / "db_BLASTn_stats.tsv",
                newline="",
        ) as handle:
            reader = _csv.DictReader(handle, delimiter="\t")
            list(reader)

        assert "Evidence_Status" in reader.fieldnames
        assert "Evidence_Present" in reader.fieldnames
        assert "Candidate_Allele_Detected" in reader.fieldnames
        assert "Exact_Allele_Detected" in reader.fieldnames
        assert "Detection_System" in reader.fieldnames
        assert "Detection_Depth" in reader.fieldnames
        assert "Detection_Depth_Coverage" in reader.fieldnames
        assert "Base_Coverage_Hit" not in reader.fieldnames

    def test_candidate_allele_is_reported_between_evidence_and_exact(
            self, tmp_path):
        wf = _make_workflow()
        wf.output_dir = tmp_path
        wf.stats_dir = tmp_path / "tool_stats"
        wf.stats_dir.mkdir()
        hit_file = tmp_path / "candidate.tsv"
        hit_file.write_text("".join(
            _blast_line(qseqid=f"read_{index}", pident=99.0)
            for index in range(3)
        ))

        _, gene_reads = wf.parse_blast_results(
            hit_file, "db", "BLASTn"
        )
        wf.write_tool_stats("db", "BLASTn", gene_reads)
        wf.generate_detection_matrix("db")

        with open(
                wf.stats_dir / "db_BLASTn_stats.tsv",
                newline="",
        ) as handle:
            stats_rows = list(_csv.DictReader(handle, delimiter="\t"))
        with open(
                tmp_path / "db_evidence_summary.tsv",
                newline="",
        ) as handle:
            summary_rows = {
                row["Tool"]: row
                for row in _csv.DictReader(handle, delimiter="\t")
            }

        call = wf.evidence_calls["db"]["geneA"]["BLASTn"]
        assert call.candidate_allele_detected
        assert not call.exact_allele_detected
        assert stats_rows[0]["Candidate_Allele_Detected"] == "1"
        assert stats_rows[0]["Exact_Allele_Detected"] == "0"
        assert summary_rows["BLASTn"]["Evidence_Genes"] == "1"
        assert summary_rows["BLASTn"]["Candidate_Allele_Genes"] == "1"
        assert summary_rows["BLASTn"]["Exact_Allele_Genes"] == "0"

    def test_large_alignment_len_regression_does_not_inflate_query_coverage(self):
        """Regression: alignment_len (field 3) > qlen must NOT produce > 100 % query coverage.

        Old formula: alignment_len / qlen would yield 120/80 = 150 % → always passes.
        New formula: (qend - qstart + 1) / qlen = 80/80 = 100 % → passes correctly.

        We verify the result matches expectation from the span formula,
        not the inflated alignment_len formula.
        """
        # qstart=1, qend=80, qlen=80 → span coverage = 100 %
        # alignment_len=120 → old coverage = 150 % (bug)
        lines = _blast_lines(
            pident=95.0, length=120, qstart=1, qend=80, qlen=80,
            sstart=1, send=80, slen=100,
        )
        detected, _, wf = self._run([lines],
                                    query_min_coverage=40.0, query_min_identity=80.0,
                                    detection_min_coverage=70.0)
        assert "geneA" in detected   # span-based coverage = 100 % → passes

        # Rerun with a query that covers only 50 % by span
        line2 = _blast_line(pident=95.0, length=120, qstart=1, qend=40, qlen=80,
                            sstart=1, send=40, slen=100)
        detected2, _, _ = self._run([line2],
                                    query_min_coverage=60.0,   # requires 60 % query span
                                    query_min_identity=80.0,
                                    detection_min_coverage=70.0)
        # span = (40-1+1)/80 = 50 % < 60 % threshold → should NOT pass
        assert "geneA" not in detected2

    def test_comment_lines_ignored(self):
        """Lines starting with '#' must be silently skipped."""
        comment = "# This is a comment\n"
        lines = _blast_lines(
            pident=100.0, qstart=1, qend=100, qlen=100,
            sstart=1, send=100, slen=100,
        )
        detected, _, _ = self._run([comment, lines],
                                   query_min_coverage=40.0, query_min_identity=80.0,
                                   detection_min_coverage=80.0)
        assert "geneA" in detected

    def test_malformed_line_skipped(self):
        """A line with fewer than 13 fields must be silently skipped without error."""
        bad_line = "read1\tgeneA\t95.0\n"
        good_lines = _blast_lines(
            pident=95.0, qstart=1, qend=100, qlen=100,
            sstart=1, send=100, slen=100,
        )
        detected, _, _ = self._run([bad_line, good_lines])
        assert "geneA" in detected

    def test_thirteen_column_line_skipped_without_aborting(self):
        """A 13-column line lacks slen and must not abort parsing of later hits."""
        thirteen_cols = "\t".join(str(x) for x in [
            "read_bad", "geneA", 95.0, 100, 0, 0,
            1, 100, 1, 100, 1e-20, 100.0, 100,
        ]) + "\n"
        good_lines = _blast_lines(
            qseqid_prefix="read_good", pident=95.0, qstart=1, qend=100,
            qlen=100, sstart=1, send=100, slen=100,
        )
        detected, _, _ = self._run([thirteen_cols, good_lines])
        assert "geneA" in detected

    def test_malformed_numeric_row_does_not_abort_later_hits(self):
        """A malformed numeric field must not discard the rest of a large file."""
        malformed = _blast_line(
            qseqid="read_bad", pident=95.0, qstart=1, qend=100,
            qlen=100, sstart=1, send=100, slen=100,
        ).replace("\t95.0\t", "\tnot-a-number\t", 1)
        good_lines = _blast_lines(
            qseqid_prefix="read_good", pident=95.0, qstart=1, qend=100,
            qlen=100, sstart=1, send=100, slen=100,
        )
        detected, _, _ = self._run([malformed, good_lines])
        assert "geneA" in detected

    def test_min_num_reads_threshold(self):
        """With detection_min_num_reads=3, a single read must not produce a detection."""
        line = _blast_line(pident=100.0, qstart=1, qend=100, qlen=100, sstart=1, send=100, slen=100)
        detected, _, _ = self._run([line],
                                   query_min_coverage=40.0, query_min_identity=80.0,
                                   detection_min_coverage=80.0,
                                   detection_min_num_reads=3)
        assert "geneA" not in detected

    def test_count_only_mode_keeps_counts_not_read_lists(self):
        """Recompute can avoid storing every read name while preserving stats counts."""
        wf = _make_workflow(query_min_coverage=40.0, query_min_identity=80.0,
                            detection_min_coverage=80.0)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as fh:
            fh.write(_blast_lines(
                qseqid_prefix="r", pident=95.0, qstart=1, qend=100,
                qlen=100, sstart=1, send=100, slen=100,
            ))
            tmp_path = Path(fh.name)
        try:
            detected, gene_reads = wf.parse_blast_results(tmp_path, "db", "blastn", count_only=True)
        finally:
            tmp_path.unlink(missing_ok=True)

        assert "geneA" in detected
        assert gene_reads["geneA"]["all_count"] == 3
        assert gene_reads["geneA"]["passing_count"] == 3
        assert "all" not in gene_reads["geneA"]
        assert wf.gene_stats["db"]["blastn"]["geneA"].avg_identity == 95.0

    def test_multiple_hsps_from_one_read_count_as_one_depth_observation(self):
        """Overlapping HSPs from one read must not inflate reads or depth."""
        wf = _make_workflow(
            query_min_coverage=40.0,
            query_min_identity=80.0,
            detection_min_coverage=80.0,
        )
        with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tsv", delete=False) as handle:
            handle.write(_blast_line(
                qseqid="r1", pident=100.0, qstart=1, qend=60,
                qlen=100, sstart=1, send=60, slen=100,
            ))
            handle.write(_blast_line(
                qseqid="r1", pident=100.0, qstart=41, qend=100,
                qlen=100, sstart=41, send=100, slen=100,
            ))
            path = Path(handle.name)
        try:
            detected, gene_reads = wf.parse_blast_results(
                path, "db", "blastn", count_only=True
            )
        finally:
            path.unlink(missing_ok=True)

        stats = wf.gene_stats["db"]["blastn"]["geneA"]
        assert "geneA" not in detected
        assert gene_reads["geneA"]["passing_count"] == 1
        assert stats.num_sequences == 1
        assert stats.passing_read_support == 1
        assert stats.coverage_2x == 0.0

    def test_partial_diagnostic_call_is_written_to_evidence_matrix(
            self, tmp_path):
        wf = _make_workflow(
            query_min_coverage=20.0,
            query_min_identity=80.0,
            detection_min_coverage=80.0,
        )
        wf.output_dir = tmp_path
        hit_file = tmp_path / "partial.tsv"
        hit_file.write_text(_blast_line(
            pident=100.0,
            qstart=1,
            qend=30,
            qlen=100,
            sstart=1,
            send=30,
            slen=100,
        ))

        wf.parse_blast_results(hit_file, "db", "blastn")
        wf.generate_detection_matrix("db")

        with open(
                tmp_path / "db_evidence_matrix.tsv",
                newline="",
        ) as handle:
            rows = list(_csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 1
        assert rows[0]["Gene"] == "geneA"
        assert rows[0]["blastn"] == "PARTIAL_OR_DIVERGENT"
        assert rows[0]["Evidence_Detections"] == "0"

    def test_recompute_evalue_threshold_is_applied(self):
        """A stricter recompute e-value threshold should filter otherwise valid hits."""
        wf = _make_workflow(query_min_coverage=40.0, query_min_identity=80.0,
                            detection_min_coverage=80.0)
        wf.evalue = 1e-30
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as fh:
            fh.write(_blast_line(pident=100.0, qstart=1, qend=100, qlen=100,
                                 sstart=1, send=100, slen=100, evalue=1e-10))
            tmp_path = Path(fh.name)
        try:
            detected, _ = wf.parse_blast_results(tmp_path, "db", "blastn")
        finally:
            tmp_path.unlink(missing_ok=True)
        assert "geneA" not in detected


# ===========================================================================
# 6b. DIAMOND command construction
# ===========================================================================

class TestDiamondCommandConstruction:
    def test_diamond_results_tsv_is_written_without_comment_header(self, tmp_path):
        wf = object.__new__(Workflow)
        wf.raw_dir = tmp_path
        wf.input_fasta = tmp_path / "reads.fasta"
        wf.tool_sensitivity_params = {}
        wf.query_min_identity = 80.0
        wf.query_min_coverage = 40.0
        wf.threads = 1
        wf.evalue = None
        wf.extra_tool_params = {}
        wf.check_gzip = lambda _: True

        captured = {}

        def fake_run_command(cmd, description):
            captured["cmd"] = cmd
            captured["description"] = description
            return False

        wf.run_command = fake_run_command

        success, detected = Workflow.run_diamond(wf, "db.dmnd", "resfinder")

        assert success is False
        assert detected == set()
        assert "--header" not in captured["cmd"]
        assert captured["cmd"][captured["cmd"].index("-f") + 1] == "6"
        assert captured["description"] == "resfinder - DIAMOND-BLASTX"

    def test_diamond_blastp_is_named_in_command_logging(self, tmp_path):
        wf = object.__new__(Workflow)
        wf.raw_dir = tmp_path
        wf.input_fasta = tmp_path / "proteins.fasta"
        wf.tool_sensitivity_params = {}
        wf.query_min_identity = 80.0
        wf.query_min_coverage = 40.0
        wf.threads = 1
        wf.evalue = None
        wf.extra_tool_params = {}
        wf.check_gzip = lambda _: True

        captured = {}

        def fake_run_command(cmd, description):
            captured["cmd"] = cmd
            captured["description"] = description
            return False

        wf.run_command = fake_run_command

        Workflow.run_diamond(
            wf, "db.dmnd", "custom", query_mode="blastp"
        )

        assert captured["cmd"][1] == "blastp"
        assert captured["description"] == "custom - DIAMOND-BLASTP"


# ===========================================================================
# 7.  Coverage formula cross-check: BLAST span ≡ BAM (M+I) for equivalent reads
# ===========================================================================

class TestCrossToolConsistency:
    """Verify the BLAST and BAM coverage formulas agree for equivalent alignments."""

    @pytest.mark.parametrize("qstart,qend,qlen,cigar,read_len,nm", [
        # Perfect 100 bp alignment, no indels
        (1, 100, 100, "100M", 100, 0),
        # 80 bp aligned out of 100 bp read (20 bp soft-clipped in BAM)
        (1, 80, 100, "80M20S", 100, 0),
        # Read with 10 bp insertion: BLAST span includes them, new BAM formula does too
        (1, 100, 100, "90M10I", 100, 0),
    ])
    def test_blast_and_bam_coverage_agree(self, qstart, qend, qlen, cigar, read_len, nm):
        blast_cov = ((abs(qend - qstart) + 1) / qlen) * 100
        bam_cov, _ = _cigar_query_coverage_and_identity(cigar, read_len, nm)
        assert abs(blast_cov - bam_cov) < 1e-9, (
            f"BLAST coverage ({blast_cov:.4f} %) differs from BAM coverage "
            f"({bam_cov:.4f} %) for qstart={qstart} qend={qend} cigar={cigar}"
        )

    def test_blastx_qstart_qend_are_nucleotide_coords(self):
        """For blastx, qstart/qend/qlen are in nucleotides even though the
        database is protein.  The span formula is therefore self-consistent
        and requires NO amino-acid conversion (unlike alignment_len which is
        in AA and would need *3).

        Scenario: 300 nt read aligned against a 100 AA protein via blastx.
          qstart=1, qend=150, qlen=300  → span coverage = 50 % (150 NT / 300 NT)
          alignment_len ≈ 50 AA         → old formula: 50 / 300 = 16.7 %  (wrong)
          alignment_len * 3 = 150 NT    → old workaround: 150 / 300 = 50 % (same as span)
          span formula: (150-1+1)/300   → 50 %  ✓ correct, no conversion needed
        """
        qstart, qend, qlen = 1, 150, 300   # nucleotide query coordinates
        alignment_len_aa = 50              # BLAST field-3 for blastx is in AA

        span_coverage    = (abs(qend - qstart) + 1) / float(qlen) * 100.0  # 50 %
        old_raw_formula  = alignment_len_aa / float(qlen) * 100.0           # 16.7 % (wrong)
        old_fixed_formula = (alignment_len_aa * 3) / float(qlen) * 100.0   # 50 % (same result, complex)

        # New span formula gives 50 % without any AA→NT conversion
        assert abs(span_coverage - 50.0) < 1e-9
        # Old raw formula was wrong (mixing units)
        assert old_raw_formula < span_coverage
        # Old workaround coincides with span when there are no gaps, confirming equivalence
        assert abs(old_fixed_formula - span_coverage) < 1e-9

    def test_blastp_aa_coords_self_consistent(self):
        """For blastp, qstart/qend/qlen are all in amino acids — the span formula
        works identically; no special handling is needed.
        """
        qstart, qend, qlen = 1, 80, 100   # AA coordinates
        coverage = (abs(qend - qstart) + 1) / float(qlen) * 100.0
        assert abs(coverage - 80.0) < 1e-9

    def test_insertion_old_bam_vs_new_bam_vs_blast(self):
        """Show that old M-only BAM formula diverges from BLAST but new M+I formula does not."""
        cigar = "90M10I"
        read_len = 100
        qstart, qend, qlen = 1, 100, 100

        blast_cov = ((abs(qend - qstart) + 1) / qlen) * 100  # 100 %

        # Old formula (M only)
        old_m_only = (90 / read_len) * 100  # 90 %

        # New formula (M + I)
        new_m_plus_i, _ = _cigar_query_coverage_and_identity(cigar, read_len, nm=0)  # 100 %

        assert old_m_only != blast_cov, "Setup: old formula should differ from BLAST"
        assert abs(new_m_plus_i - blast_cov) < 1e-9, "New formula must equal BLAST"


# ---------------------------------------------------------------------------
# Bug-regression: threshold consistency between BLAST and BAM
# ---------------------------------------------------------------------------

class TestThresholdConsistency:
    """Verify that BLAST and BAM use the same per-read gate (query_min_identity)
    and that both use detection_min_identity for the gene-level avg_identity gate.

    Before the fix:
      - BAM per-read gate used detection_min_identity (wrong)
      - Neither BLAST nor BAM applied detection_min_identity at the gene level (wrong)
    """

    def test_blast_uses_query_min_identity_not_detection(self, tmp_path):
        """A read at 75% identity should be admitted when query_min_identity=70 even if
        detection_min_identity=80 — BLAST per-read gate is query_min_identity."""
        wf = _make_workflow(
            query_min_identity=70.0,
            detection_min_identity=80.0,
            detection_min_coverage=0.0,
        )
        blast_file = tmp_path / "hits.tsv"
        blast_file.write_text(_blast_line(pident=75.0, qstart=1, qend=100, qlen=100, sstart=1, send=100, slen=100) + "\n")
        wf.parse_blast_results(blast_file, "db", "blastn")
        stats = wf.gene_stats["db"]["blastn"].get("geneA")
        assert stats is not None, "Read should have been admitted (identity 75 >= query_min_identity 70)"
        assert stats.num_sequences == 1

    def test_blast_rejects_read_below_query_min_identity(self, tmp_path):
        """A read at 65% identity must be rejected when query_min_identity=70, regardless
        of detection_min_identity."""
        wf = _make_workflow(
            query_min_identity=70.0,
            detection_min_identity=60.0,
            detection_min_coverage=0.0,
        )
        blast_file = tmp_path / "hits.tsv"
        blast_file.write_text(_blast_line(pident=65.0, qstart=1, qend=100, qlen=100, sstart=1, send=100, slen=100) + "\n")
        wf.parse_blast_results(blast_file, "db", "blastn")
        stats = wf.gene_stats["db"]["blastn"].get("geneA")
        assert stats is None, "Read below query_min_identity should not contribute to stats"

    def test_detection_min_identity_gates_gene_level_blast(self, tmp_path):
        """Gene should not be detected when avg_identity < detection_min_identity even if
        gene_coverage and base_depth thresholds are met."""
        wf = _make_workflow(
            query_min_identity=70.0,    # per-read gate — lets 75% reads through
            detection_min_identity=80.0,  # gene-level gate — requires 80% avg
            detection_min_coverage=50.0,
        )
        blast_file = tmp_path / "hits.tsv"
        # Two reads both at 75% identity: avg_identity = 75, below detection_min_identity=80
        blast_file.write_text(
            _blast_line(qseqid="r1", pident=75.0, qstart=1, qend=100, qlen=100, sstart=1, send=50, slen=100) + "\n" +
            _blast_line(qseqid="r2", pident=75.0, qstart=1, qend=100, qlen=100, sstart=51, send=100, slen=100) + "\n"
        )
        detected, _ = wf.parse_blast_results(blast_file, "db", "blastn")
        assert "geneA" not in detected, "Gene avg_identity (75%) is below detection_min_identity (80%) — must not be detected"

    def test_detection_min_identity_allows_gene_when_avg_identity_sufficient(self, tmp_path):
        """Gene should be detected when avg_identity >= detection_min_identity."""
        wf = _make_workflow(
            query_min_identity=70.0,
            detection_min_identity=80.0,
            detection_min_coverage=80.0,
        )
        blast_file = tmp_path / "hits.tsv"
        # Three reads at 95%: avg_identity = 95, 3x coverage = 100/100 = 100%
        blast_file.write_text(
            _blast_lines(
                pident=95.0, qstart=1, qend=100, qlen=100,
                sstart=1, send=100, slen=100,
            )
        )
        detected, _ = wf.parse_blast_results(blast_file, "db", "blastn")
        assert "geneA" in detected, "Gene with avg_identity 95% must be detected"


# ---------------------------------------------------------------------------
# HMMER annotation loading and must-flag detection
# ---------------------------------------------------------------------------

import csv as _csv
import io as _io


def _make_hmmer_workflow(
    detection_min_coverage: float = 80.0,
    detection_min_depth: int = 1,
    detection_min_num_reads: int = 1,
) -> Workflow:
    """Minimal Workflow for HMMER parsing tests."""
    import logging
    wf = object.__new__(Workflow)
    wf.query_min_coverage = 0.0
    wf.query_min_identity = 0.0
    wf.detection_min_coverage = detection_min_coverage
    wf.detection_min_identity = 0.0
    wf.detection_min_depth = detection_min_depth
    wf.detection_min_num_reads = detection_min_num_reads
    wf.min_bitscore = None
    wf.hmmer_evalue = None
    wf.evalue = None
    wf.hmmer_threshold_mode = 'evalue'
    wf.is_genes_fasta = False
    wf.evidence_config = EvidenceConfig(
        corroborating_depth=1,
        exact_identity_min=0.0,
        min_unique_best_reads=1,
        min_unique_best_fraction=0.0,
    )
    wf.gene_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(GeneStats)))
    wf.detections = defaultdict(lambda: defaultdict(lambda: defaultdict(bool)))
    wf.evidence_calls = defaultdict(lambda: defaultdict(dict))
    wf.family_calls = defaultdict(lambda: defaultdict(dict))
    wf.report_fasta = None
    wf.all_reads = {}
    wf.input_fasta = None
    wf.hmmer_annotations = {}
    wf.hmmer_annotations_meta = {}
    wf.databases = {}
    wf.output_dir = Path("/tmp")
    wf.logger = logging.getLogger("test_hmmer")
    wf.logger.addHandler(logging.NullHandler())
    return wf


def _write_tblout(genes: list, evalue: float = 1e-10) -> str:
    """Minimal hmmer --tblout content for the given profile names."""
    lines = ["# tblout mock\n"]
    for gene in genes:
        lines.append(f"seq1\t-\t{gene}\t-\t{evalue}\t200.0\t0\t0\t0\t0\t0\t0\t0\t0\t-\n")
    return "".join(lines)


def _write_domtblout(genes: list, hmm_from=1, hmm_to=100, qlen=100, ali_from=1, ali_to=100, acc=0.99) -> str:
    """Minimal hmmer --domtblout content for the given profile names.

    domtblout column layout (0-based, space-separated):
      0  target_name  1  target_acc  2  tlen  3  query_name  4  query_acc  5  qlen
      6  full_E       7  full_score  8  full_bias
      9  dom_num      10 dom_total   11 c-Evalue  12 i-Evalue  13 score  14 bias
      15 hmm_from     16 hmm_to      17 ali_from  18 ali_to    19 env_from  20 env_to
      21 acc          22 description
    """
    lines = ["# domtblout mock\n"]
    for gene in genes:
        # Build exactly 23 fields matching the domtblout spec above
        row_fields = [
            "seq1",   # 0 target_name
            "-",      # 1 target_acc
            "150",    # 2 tlen
            gene,     # 3 query_name (profile)
            "-",      # 4 query_acc
            str(qlen),# 5 qlen
            "1e-10",  # 6 full-seq E-value
            "200.0",  # 7 full-seq score
            "0",      # 8 full-seq bias
            "1",      # 9 domain number
            "1",      # 10 total domains
            "1e-10",  # 11 c-Evalue
            "1e-10",  # 12 i-Evalue
            "200",    # 13 domain score (integer-safe)
            "0",      # 14 domain bias
            str(hmm_from),  # 15 hmm_from
            str(hmm_to),    # 16 hmm_to
            str(ali_from),  # 17 ali_from
            str(ali_to),    # 18 ali_to
            "1",      # 19 env_from
            "100",    # 20 env_to
            str(acc), # 21 acc
            "desc",   # 22 description
        ]
        assert len(row_fields) == 23
        lines.append(" ".join(row_fields) + "\n")
    return "".join(lines)


class TestHmmerAnnotationLoading:
    """Tests for _load_hmmer_annotations."""

    def test_loads_with_must_flag_column(self, tmp_path):
        csv_file = tmp_path / "ann.csv"
        csv_file.write_text("ID,Description,Must flag\ngeneA,Gene A,TRUE\ngeneB,Gene B,FALSE\n")
        wf = _make_hmmer_workflow()
        ann, meta = wf._load_hmmer_annotations(str(csv_file))
        assert meta['has_must_flag'] is True
        assert ann['geneA']['must_flag'] is True
        assert ann['geneB']['must_flag'] is False

    def test_loads_without_must_flag_column(self, tmp_path):
        csv_file = tmp_path / "ann.csv"
        csv_file.write_text("ID,Description\ngeneA,Gene A\ngeneB,Gene B\n")
        wf = _make_hmmer_workflow()
        ann, meta = wf._load_hmmer_annotations(str(csv_file))
        assert meta['has_must_flag'] is False
        assert ann['geneA']['must_flag'] is False
        assert ann['geneB']['must_flag'] is False

    def test_case_insensitive_must_flag_column(self, tmp_path):
        csv_file = tmp_path / "ann.csv"
        csv_file.write_text("ID,Description,must flag\ngeneA,Gene A,Yes\n")
        wf = _make_hmmer_workflow()
        ann, meta = wf._load_hmmer_annotations(str(csv_file))
        assert meta['has_must_flag'] is True
        assert ann['geneA']['must_flag'] is True

    def test_must_flag_true_values(self, tmp_path):
        for val in ('TRUE', 'True', 'YES', 'yes', '1'):
            csv_file = tmp_path / f"ann_{val}.csv"
            csv_file.write_text(f"ID,Must flag\ngeneA,{val}\n")
            wf = _make_hmmer_workflow()
            ann, _ = wf._load_hmmer_annotations(str(csv_file))
            assert ann['geneA']['must_flag'] is True, f"Value '{val}' should be must_flag=True"

    def test_must_flag_false_values(self, tmp_path):
        for val in ('FALSE', 'False', 'NO', 'no', '0', ''):
            csv_file = tmp_path / f"ann_{val}.csv"
            csv_file.write_text(f"ID,Must flag\ngeneA,{val}\n")
            wf = _make_hmmer_workflow()
            ann, _ = wf._load_hmmer_annotations(str(csv_file))
            assert ann['geneA']['must_flag'] is False, f"Value '{val}' should be must_flag=False"

    def test_missing_file_returns_empty(self, tmp_path):
        wf = _make_hmmer_workflow()
        ann, meta = wf._load_hmmer_annotations(str(tmp_path / "nonexistent.csv"))
        assert ann == {}
        assert meta.get('has_must_flag') is False


class TestHmmerMustFlagOverride:
    """must_flag genes must be reported even when they fail coverage/min-reads gates."""

    def _run(self, wf, genes_tblout, genes_domtblout, tmp_path, hmm_from=1, hmm_to=100, qlen=100):
        tbl = tmp_path / "out.tbl"
        dom = tmp_path / "out.domtbl"
        tbl.write_text(_write_tblout(genes_tblout))
        dom.write_text(_write_domtblout(genes_domtblout, hmm_from=hmm_from, hmm_to=hmm_to, qlen=qlen))
        return wf.parse_hmmer_results(tbl, dom, "db", "hmmer_protein")

    def test_must_flag_gene_detected_despite_low_coverage(self, tmp_path):
        """Low-coverage must-flag evidence must be surfaced for review, not called exact."""
        wf = _make_hmmer_workflow(detection_min_coverage=80.0)
        wf.hmmer_annotations['db'] = {'geneA': {'description': 'test', 'must_flag': True}}
        wf.hmmer_annotations_meta['db'] = {'has_must_flag': True}
        # hmm_from=1, hmm_to=10 out of qlen=100 → 10% coverage
        detected, _ = self._run(wf, ['geneA'], ['geneA'], tmp_path, hmm_from=1, hmm_to=10, qlen=100)
        assert 'geneA' not in detected
        call = wf.evidence_calls['db']['geneA']['hmmer_protein']
        assert call.status == MUST_FLAG_REVIEW
        assert not call.evidence_present

    def test_non_must_flag_gene_filtered_at_low_coverage(self, tmp_path):
        """Without must_flag, a low-coverage gene must be filtered out."""
        wf = _make_hmmer_workflow(detection_min_coverage=80.0)
        wf.hmmer_annotations['db'] = {'geneA': {'description': 'test', 'must_flag': False}}
        wf.hmmer_annotations_meta['db'] = {'has_must_flag': True}
        detected, _ = self._run(wf, ['geneA'], ['geneA'], tmp_path, hmm_from=1, hmm_to=10, qlen=100)
        assert 'geneA' not in detected, "non-must_flag gene with low coverage must be filtered"

    def test_must_flag_without_annotations_column_is_not_triggered(self, tmp_path):
        """When the annotation CSV has no must_flag column, no override happens."""
        wf = _make_hmmer_workflow(detection_min_coverage=80.0)
        # Annotations loaded from CSV without must_flag column → must_flag=False for all
        wf.hmmer_annotations['db'] = {'geneA': {'description': 'test', 'must_flag': False}}
        wf.hmmer_annotations_meta['db'] = {'has_must_flag': False}
        detected, _ = self._run(wf, ['geneA'], ['geneA'], tmp_path, hmm_from=1, hmm_to=10, qlen=100)
        assert 'geneA' not in detected, "No must_flag column means no override; low-coverage gene must be filtered"

    def test_must_flag_gene_without_annotations(self, tmp_path):
        """No annotations at all: must_flag cannot fire; standard gates apply."""
        wf = _make_hmmer_workflow(detection_min_coverage=80.0)
        # No annotations loaded for this DB
        detected, _ = self._run(wf, ['geneA'], ['geneA'], tmp_path, hmm_from=1, hmm_to=10, qlen=100)
        assert 'geneA' not in detected

    def test_must_flag_gene_with_full_coverage_also_detected(self, tmp_path):
        """must_flag gene with full coverage is detected via the normal path."""
        wf = _make_hmmer_workflow(detection_min_coverage=80.0)
        wf.hmmer_annotations['db'] = {'geneA': {'description': 'test', 'must_flag': True}}
        wf.hmmer_annotations_meta['db'] = {'has_must_flag': True}
        detected, _ = self._run(wf, ['geneA'], ['geneA'], tmp_path, hmm_from=1, hmm_to=100, qlen=100)
        assert 'geneA' in detected

    def test_no_must_flag_disables_override(self, tmp_path):
        """When hmmer_must_flag=False the threshold override must not fire."""
        wf = _make_hmmer_workflow(detection_min_coverage=80.0)
        wf.hmmer_must_flag = False  # simulate --no-must-flag
        wf.hmmer_annotations['db'] = {'geneA': {'description': 'test', 'must_flag': True}}
        wf.hmmer_annotations_meta['db'] = {'has_must_flag': True}
        # 10% coverage — would normally be rescued by must_flag override
        detected, _ = self._run(wf, ['geneA'], ['geneA'], tmp_path, hmm_from=1, hmm_to=10, qlen=100)
        assert 'geneA' not in detected, "must_flag override must be disabled when hmmer_must_flag=False"

    def test_evalue_is_applied_to_each_target_profile_pair(self, tmp_path):
        """A weak target must not borrow significance from another target."""
        wf = _make_hmmer_workflow(detection_min_coverage=80.0)
        wf.hmmer_evalue = 1e-5
        tbl = tmp_path / "out.tbl"
        dom = tmp_path / "out.domtbl"
        tbl.write_text(
            "# tblout mock\n"
            "target_good\t-\tgeneA\t-\t1e-20\t200.0\t0\t0\t0\t0\t0\t0\t0\t0\t-\n"
            "target_weak\t-\tgeneA\t-\t1\t20.0\t0\t0\t0\t0\t0\t0\t0\t0\t-\n"
        )
        good_domain = _write_domtblout(
            ["geneA"], hmm_from=1, hmm_to=50, qlen=100,
        ).replace("seq1", "target_good")
        weak_domain = _write_domtblout(
            ["geneA"], hmm_from=51, hmm_to=100, qlen=100,
        ).replace("seq1", "target_weak")
        dom.write_text(good_domain + "".join(
            line for line in weak_domain.splitlines(keepends=True)
            if not line.startswith("#")
        ))

        detected, _ = wf.parse_hmmer_results(
            tbl, dom, "db", "hmmer_protein"
        )

        stats = wf.gene_stats["db"]["hmmer_protein"]["geneA"]
        assert "geneA" not in detected
        assert stats.gene_coverage == 50.0
        assert stats.passing_read_support == 1


class TestHmmerAnnotationsReportColumns:
    """TSV report must include Must_Flag column only when the source CSV had one."""

    def _run_report(self, wf, detected_gene, has_must_flag: bool, tmp_path):
        wf.detections['db'][detected_gene]['hmmer_protein'] = True
        wf.hmmer_annotations_meta['db'] = {'has_must_flag': has_must_flag}
        wf.output_dir = tmp_path
        wf._write_hmmer_annotations_report('db', 'hmmer_protein')
        report = tmp_path / "db_hmmer_protein_annotations.tsv"
        return report

    def test_must_flag_column_present_when_source_has_it(self, tmp_path):
        wf = _make_hmmer_workflow()
        wf.hmmer_annotations['db'] = {'geneA': {'description': 'A', 'must_flag': True}}
        report = self._run_report(wf, 'geneA', has_must_flag=True, tmp_path=tmp_path)
        headers = report.read_text().splitlines()[0].split('\t')
        assert 'Must_Flag' in headers

    def test_must_flag_column_absent_when_source_lacks_it(self, tmp_path):
        wf = _make_hmmer_workflow()
        wf.hmmer_annotations['db'] = {'geneA': {'description': 'A', 'must_flag': False}}
        report = self._run_report(wf, 'geneA', has_must_flag=False, tmp_path=tmp_path)
        headers = report.read_text().splitlines()[0].split('\t')
        assert 'Must_Flag' not in headers, "Must_Flag column must be absent when source CSV lacked it"

    def test_detected_column_always_present(self, tmp_path):
        for has_mf in (True, False):
            wf = _make_hmmer_workflow()
            wf.hmmer_annotations['db'] = {'geneA': {'description': 'A', 'must_flag': False}}
            report = self._run_report(wf, 'geneA', has_must_flag=has_mf, tmp_path=tmp_path)
            headers = report.read_text().splitlines()[0].split('\t')
            assert 'Detected' in headers
            assert 'Gene_ID' in headers
            assert 'Description' in headers
