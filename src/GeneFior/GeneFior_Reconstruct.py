"""
GeneFior_Reconstruct.py — In-Sample Gene Reconstruction
========================================================
Takes the output of a GeneFior alignment run (SAM/BAM files) **or** a
BLAST / DIAMOND tabular file and reconstructs the 'true' sequence(s) of a
gene as they exist in the sample, rather than in the reference database.

Key features
------------
* **Consensus reconstruction** — pile reads up against the reference gene
  and call the most-supported base at every position.
* **Boundary extension** — attempt to reconstruct ``±boundary_extension``
  (default 20 %) beyond both ends of the reference gene.  Extension windows
  are rounded to codon boundaries (multiples of 3) to maintain reading frame.
  Soft-clipped read bases supply the extension sequence.
* **Length & codon assessment** — after reconstruction the sequence is
  compared to the reference length and assessed for valid ORF structure
  (length divisible by 3, canonical start / stop codons).  Genes that are
  genuinely shorter in the sample but retain a valid ORF are reported as
  ``SHORTER_IN_SAMPLE`` without a grade penalty.  Genes that extend beyond
  the reference are reported as ``LONGER_IN_SAMPLE``.
* **Multi-version detection** — examine the allele-frequency spectrum to
  classify the locus as *single*, *multi* (≥ 2 versions), or *uncertain*.
* **Haplotype separation** — when multiple versions are detected, assign
  reads to haplotypes by scoring their allele profiles against the majority
  and minority allele distributions, then build a separate consensus per
  haplotype.
* **Grading** — assign an overall reconstruction quality grade:

  =========  ================================================================
  Grade      Meaning
  =========  ================================================================
  A          ≥ 95 % coverage, ≥ 95 % read support, low depth-CV, < 5 % N
  B          ≥ 80 % coverage, ≥ 85 % read support
  C          ≥ 60 % coverage, ≥ 70 % read support (or uncertain multi-version)
  C*         Multiple versions detected — consensus is a blend; use haplotypes
  F          < 60 % coverage or < 70 % read support
  =========  ================================================================

  When a gene is ``SHORTER_IN_SAMPLE`` with a valid ORF the coverage
  denominator is adjusted to the reconstructed length so the grade reflects
  quality *over the covered span*, not against the full reference length.

Usage examples
--------------
SAM / BAM mode (output from a GeneFior run)::

    GeneFior-Reconstruct \\
        -output_dir  /path/to/genefior_results/ \\
        -gene        tet(Q)_1_L33696 \\
        -reference_fasta /path/to/resfinder_db.fa

BLAST / DIAMOND mode::

    GeneFior-Reconstruct \\
        -blast_tsv  /path/to/hits.tsv \\
        -gene       tet(Q)_1_L33696 \\
        -reference_fasta /path/to/resfinder_db.fa

Custom boundary extension and codon settings::

    GeneFior-Reconstruct \\
        -output_dir /path/to/genefior_results/ \\
        -gene tet(Q)_1_L33696 \\
        -boundary_extension 0.30 \\
        -start_codons ATG \\
        -stop_codons TAA,TAG,TGA

Outputs (under ``-recon_dir``, default: ``<output_dir>/reconstruction/``)
--------------------------------------------------------------------------
``{gene}_consensus.fasta``
    Reconstructed sample sequence(s).  Includes boundary-extension bases where
    reads extend beyond the reference.  FASTA header tags annotate the
    reference length, length-vs-reference status, codon validity, and
    extension captured at each end.

``{gene}_haplotypes.fasta``
    Per-haplotype FASTA when multi-version signal is detected.

``{gene}_depth.tsv``
    Per-position depth, consensus base, and allele frequencies.

``{gene}_variants.tsv``
    Positions where the sample consensus differs from the reference.

``{gene}_sample_vs_reference.txt``
    Detailed sample-vs-reference divergence report.  Divergence is *expected*
    biology — the reference is used only as an alignment scaffold.

``{gene}_reads.fasta``
    All reads mapped to this gene.

``{gene}_validation.txt``
    Reconstruction validation report (grade, coverage, length assessment).

``{gene}_artifact_validation.txt``
    Deeper biological-plausibility checks (ORF integrity, GC, strand bias …).

``{gene}_assembly_graph.gfa``
    Assembly graph in GFA format for Bandage visualisation.

``{gene}_reconstruction_plot.html``
    Interactive reconstruction overview plot (coverage, allele frequencies,
    divergence track, optional haplotype depth traces).

``reconstruction_report.txt``
    Summary report across all requested genes.
"""

import argparse
import gzip
import logging
import math
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# ── version ────────────────────────────────────────────────────────────────
RECONSTRUCT_VERSION = "v2.2.0"

# ── allow direct execution ──────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

if __package__:
    from .constants import GENEFIOR_VERSION  # type: ignore
    from .GeneFior_Validate import validate_reconstruction  # type: ignore
else:
    try:
        from constants import GENEFIOR_VERSION  # type: ignore
    except ImportError:
        GENEFIOR_VERSION = "unknown"
    try:
        from GeneFior_Validate import validate_reconstruction  # type: ignore
    except ImportError:
        validate_reconstruction = None

# ── logging ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("GeneFior.Reconstruct")

# ── CIGAR regex ─────────────────────────────────────────────────────────────
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

# IUPAC ambiguity codes
_IUPAC = {
    frozenset("AC"): "M", frozenset("AG"): "R", frozenset("AT"): "W",
    frozenset("CG"): "S", frozenset("CT"): "Y", frozenset("GT"): "K",
    frozenset("ACG"): "V", frozenset("ACT"): "H", frozenset("AGT"): "D",
    frozenset("CGT"): "B", frozenset("ACGT"): "N",
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Low-level helpers
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _samtools_available() -> bool:
    return shutil.which("samtools") is not None


def _open_sam_source(bam_path: Path, sam_path: Path):
    """Return (iterator-over-SAM-lines, process-or-None).
    Prefer sorted BAM via samtools; fall back to .sam."""
    if _samtools_available() and bam_path.exists():
        proc = subprocess.Popen(
            ["samtools", "view", "-h", str(bam_path)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, start_new_session=True,
        )
        logger.info(f"  Reading alignments via samtools: {bam_path.name}")
        return proc.stdout, proc
    elif sam_path.exists():
        logger.info(f"  Reading alignments from SAM: {sam_path.name}")
        return open(sam_path, "r"), None
    return None, None


def _load_reference_fasta(fasta_path: str) -> Dict[str, str]:
    """Load all sequences from a FASTA (plain or gzipped)."""
    seqs: Dict[str, str] = {}
    if not fasta_path or not os.path.exists(fasta_path):
        return seqs
    opener = gzip.open if fasta_path.endswith((".gz", ".gzip")) else open
    name = None
    parts: List[str] = []
    with opener(fasta_path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
        if name:
            seqs[name] = "".join(parts)
    return seqs


def _parse_md_tag(md: str) -> list:
    tokens = []
    i = 0
    while i < len(md):
        c = md[i]
        if c == "^":
            j = i + 1
            while j < len(md) and md[j].isalpha():
                j += 1
            tokens.append(("del", md[i + 1:j]))
            i = j
        elif c.isdigit():
            j = i
            while j < len(md) and md[j].isdigit():
                j += 1
            tokens.append(("match", int(md[i:j])))
            i = j
        elif c.isalpha():
            tokens.append(("mismatch", c.upper()))
            i += 1
        else:
            i += 1
    return tokens


def _reconstruct_ref_bases(seq: str, cigar_str: str, md_str: str,
                           ref_start: int) -> Dict[int, str]:
    """Use read sequence + CIGAR + MD tag to reconstruct reference bases."""
    ref_bases: Dict[int, str] = {}
    if not cigar_str or cigar_str == "*" or not md_str or seq == "*":
        return ref_bases
    cigar_ops = CIGAR_RE.findall(cigar_str)
    md_tokens = _parse_md_tag(md_str)
    md_idx = md_offset = 0
    ref_pos = ref_start
    read_pos = 0
    for count_str, op in cigar_ops:
        try:
            length = int(count_str)
        except ValueError:
            logger.warning(f"Invalid CIGAR count: {count_str}")
            continue
            
        if op in ("M", "=", "X"):
            for k in range(length):
                if read_pos >= len(seq):
                    break
                read_base = seq[read_pos]
                ref_base = None
                while md_idx < len(md_tokens):
                    kind, val = md_tokens[md_idx]
                    if kind == "match":
                        if md_offset < val:
                            ref_base = read_base
                            md_offset += 1
                            if md_offset == val:
                                md_idx += 1; md_offset = 0
                            break
                        else:
                            md_idx += 1; md_offset = 0
                    elif kind == "mismatch":
                        ref_base = val
                        md_idx += 1; md_offset = 0
                        break
                    else:
                        md_idx += 1; md_offset = 0
                if ref_base is not None:
                    ref_bases[ref_pos] = ref_base
                ref_pos += 1; read_pos += 1
        elif op == "I":
            read_pos += length
        elif op == "D":
            for _ in range(length):
                while md_idx < len(md_tokens):
                    kind, val = md_tokens[md_idx]
                    if kind == "del":
                        if md_offset < len(val):
                            ref_bases[ref_pos] = val[md_offset]
                            md_offset += 1
                            if md_offset == len(val):
                                md_idx += 1; md_offset = 0
                            break
                        else:
                            md_idx += 1; md_offset = 0
                    else:
                        break
                ref_pos += 1
        elif op == "N":
            ref_pos += length
        elif op == "S":
            read_pos += length

    return ref_bases


def _iupac_ambiguity(votes: Dict[str, int], depth: int,
                     threshold: float = 0.10) -> str:
    present = frozenset(
        b for b, c in votes.items()
        if b in "ACGT" and c / depth >= threshold
    )
    if not present:
        return "N"
    if len(present) == 1:
        return next(iter(present))
    return _IUPAC.get(present, "N")


def _safe_filename(s: str) -> str:
    return re.sub(r"[^\w\-.]", "_", s)


def _pct_identity(seq_a: str, seq_b: str) -> float:
    """Simple pairwise identity over aligned positions (no N in seq_a)."""
    if not seq_a or not seq_b:
        return 0.0
    length = min(len(seq_a), len(seq_b))
    matches = sum(
        1 for i in range(length)
        if seq_a[i] == seq_b[i] and seq_a[i] not in ("N", "-")
    )
    valid = sum(1 for i in range(length) if seq_a[i] not in ("N", "-"))
    return matches / valid * 100.0 if valid > 0 else 0.0


_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Codon-structure helpers
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_START_CODONS = frozenset({"ATG", "GTG", "TTG"})
_STOP_CODONS  = frozenset({"TAA", "TAG", "TGA"})

# Human-readable defaults shown in help text
_DEFAULT_START_CODONS_STR = "ATG,GTG,TTG"
_DEFAULT_STOP_CODONS_STR  = "TAA,TAG,TGA"


def _parse_codon_list(raw: str, kind: str) -> frozenset:
    """Parse a comma-separated codon string, validate each entry, and return a frozenset."""
    codons = frozenset(c.strip().upper() for c in raw.split(",") if c.strip())
    bad = [c for c in codons if len(c) != 3 or not c.isalpha()]
    if bad:
        raise ValueError(
            f"Invalid {kind} codon(s): {bad}. "
            "Each must be exactly 3 letters (e.g. ATG,GTG,TTG)."
        )
    return codons


def _round_up_to_codon(n: int) -> int:
    """Round *n* up to the nearest multiple of 3 (codon boundary)."""
    r = n % 3
    return n if r == 0 else n + (3 - r)


def _assess_codon_structure(
    seq: str,
    start_codons: Optional[frozenset] = None,
    stop_codons:  Optional[frozenset] = None,
) -> dict:
    """
    Check whether *seq* looks like a valid coding sequence.

    Leading / trailing N or '-' characters are stripped before the check so
    that low-coverage termini do not confuse the codon-boundary test.

    *start_codons* and *stop_codons* are frozensets of accepted triplets.
    When ``None`` the module defaults (_START_CODONS / _STOP_CODONS) are used.

    Returns a dict:
      is_codon_multiple  – True when len(stripped) % 3 == 0
      has_start_codon    – True when first codon is in start_codons
      has_stop_codon     – True when last  codon is in stop_codons
      start_codon        – the actual first triplet, "" if absent
      stop_codon         – the actual last  triplet, "" if absent
      is_valid_orf       – all three conditions satisfied
      length_nt          – length of the stripped sequence
    """
    sc = start_codons if start_codons is not None else _START_CODONS
    tc = stop_codons  if stop_codons  is not None else _STOP_CODONS
    s = seq.upper().strip("N").strip("-")
    is_mult   = len(s) % 3 == 0
    start_c   = s[:3]  if len(s) >= 3 else ""
    stop_c    = s[-3:] if len(s) >= 3 else ""
    has_start = start_c in sc
    has_stop  = stop_c  in tc
    return {
        "is_codon_multiple": is_mult,
        "has_start_codon":   has_start,
        "has_stop_codon":    has_stop,
        "start_codon":       start_c if has_start else "",
        "stop_codon":        stop_c  if has_stop  else "",
        "is_valid_orf":      is_mult and has_start and has_stop,
        "length_nt":         len(s),
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Extended-boundary pileup (captures soft-clipped bases beyond the reference)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _build_extended_pileup(
    alignments: List[dict],
    gene_length: int,
    ext_5: int,
    ext_3: int,
) -> Tuple[Dict[int, Dict[str, int]], Dict[int, Dict[str, int]]]:
    """
    Like ``_build_pileup`` but extends the analysis window *ext_5* bp upstream
    and *ext_3* bp downstream.  Soft-clipped read bases are captured as
    potential sequence beyond the reference gene boundaries.

    Extended coordinate convention::

        index 0                              → start of 5ʹ extension
        index ext_5                          → first base of reference gene
        index ext_5 + gene_length            → first base of 3ʹ extension
        index ext_5 + gene_length + ext_3 - 1 → last position

    Returns (base_pileup, ins_pileup) keyed by extended indices.
    """
    ext_total = ext_5 + gene_length + ext_3
    base_pileup: Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    ins_pileup:  Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for aln in alignments:
        seq   = aln["seq"]
        cigar = aln["cigar"]
        pos   = aln["pos"]          # 0-based original reference coordinate
        if not seq or seq == "*" or not cigar or cigar == "*":
            continue
        cigar_ops = CIGAR_RE.findall(cigar)
        if not cigar_ops:
            continue

        # ── leading soft clip: bases BEFORE the mapped region ────────────────
        # CIGAR "5S95M" → positions (pos-5)..(pos-1) in original ref coords.
        leading_soft = int(cigar_ops[0][0]) if cigar_ops[0][1] == "S" else 0
        for k in range(leading_soft):
            ext_p = ext_5 + (pos - leading_soft + k)   # may fall outside window
            if 0 <= ext_p < ext_total:
                rp = k
                if rp < len(seq):
                    b = seq[rp].upper()
                    if b not in "ACGTN":
                        b = "N"
                    base_pileup[ext_p][b] += 1

        # ── main CIGAR processing ─────────────────────────────────────────────
        ref_pos        = ext_5 + pos     # start in extended coordinate space
        read_pos       = 0
        last_ref       = ref_pos
        first_consumed = False           # True once we've passed any leading S

        for count_str, op in cigar_ops:
            length = int(count_str)
            if op in ("M", "=", "X"):
                first_consumed = True
                for k in range(length):
                    ext_p = ref_pos + k
                    if 0 <= ext_p < ext_total:
                        rp = read_pos + k
                        if rp < len(seq):
                            b = seq[rp].upper()
                            if b not in "ACGTN":
                                b = "N"
                            base_pileup[ext_p][b] += 1
                last_ref = ref_pos + length - 1
                ref_pos  += length
                read_pos += length
            elif op == "I":
                first_consumed = True
                ins_seq = seq[read_pos:read_pos + length].upper()
                ins_pileup[last_ref][ins_seq] += 1
                read_pos += length
            elif op == "D":
                first_consumed = True
                for k in range(length):
                    ext_p = ref_pos + k
                    if 0 <= ext_p < ext_total:
                        base_pileup[ext_p]["-"] += 1
                last_ref = ref_pos + length - 1
                ref_pos  += length
            elif op == "N":
                first_consumed = True
                ref_pos += length
            elif op == "S":
                if not first_consumed:
                    # Leading S — already captured above; just advance read_pos
                    first_consumed = True
                    read_pos += length
                else:
                    # Trailing soft clip: bases AFTER the mapped region.
                    # S does NOT consume the reference, so ref_pos stays put.
                    for k in range(length):
                        ext_p = ref_pos + k
                        if 0 <= ext_p < ext_total:
                            rp = read_pos + k
                            if rp < len(seq):
                                b = seq[rp].upper()
                                if b not in "ACGTN":
                                    b = "N"
                                base_pileup[ext_p][b] += 1
                    read_pos += length

    return base_pileup, ins_pileup


def _assess_length_vs_reference(
    final_seq: str,
    ref_length: int,
    captured_5_nt: int,
    captured_3_nt: int,
    start_codons: Optional[frozenset] = None,
    stop_codons:  Optional[frozenset] = None,
) -> dict:
    """
    Compare the reconstructed sequence length to the reference gene length and
    assess codon structure.

    *final_seq* is the full output sequence (may include extension bases and
    internal Ns for low-coverage positions).  Terminal Ns are stripped before
    the codon check.

    *start_codons* / *stop_codons* — frozensets of accepted triplets; ``None``
    falls back to the module defaults.

    Returns a dict with the length status, codon-structure flags, and a
    human-readable note suitable for inclusion in the validation report.
    """
    codon_info = _assess_codon_structure(final_seq, start_codons, stop_codons)
    recon_len  = codon_info["length_nt"]          # stripped length
    diff_nt    = recon_len - ref_length
    diff_pct   = round(diff_nt / ref_length * 100, 1) if ref_length > 0 else 0.0

    if diff_nt == 0:
        status = "MATCHES_REFERENCE"
        note   = f"Reconstructed length equals reference length ({ref_length} nt)."
    elif diff_nt > 0:
        status = "LONGER_IN_SAMPLE"
        note   = (
            f"Sample gene is LONGER than reference by {diff_nt} nt "
            f"({diff_pct:+.1f}%). "
            f"Extension captured: 5ʹ +{captured_5_nt} nt, 3ʹ +{captured_3_nt} nt. "
            "This likely reflects a genuine longer version of the gene in this sample."
        )
    else:
        status = "SHORTER_IN_SAMPLE"
        if codon_info["is_valid_orf"]:
            note = (
                f"Sample gene is SHORTER than reference by {-diff_nt} nt "
                f"({diff_pct:.1f}%). "
                f"Codon structure is VALID "
                f"(start: {codon_info['start_codon']}, stop: {codon_info['stop_codon']}). "
                "This likely reflects a genuine shorter version of the gene in this sample. "
                "Grade is not penalised for length; see coverage metrics for reconstruction quality."
            )
        else:
            issues = []
            if not codon_info["is_codon_multiple"]:
                issues.append(f"length {recon_len} not divisible by 3")
            if not codon_info["has_start_codon"]:
                issues.append(
                    f"no canonical start codon (got {codon_info['start_codon']!r})"
                )
            if not codon_info["has_stop_codon"]:
                issues.append(
                    f"no canonical stop codon (got {codon_info['stop_codon']!r})"
                )
            note = (
                f"Sample gene is SHORTER than reference by {-diff_nt} nt "
                f"({diff_pct:.1f}%). "
                f"Codon structure INVALID: {'; '.join(issues) or 'unknown'}. "
                "May indicate partial coverage or genuine truncation."
            )

    if codon_info["is_valid_orf"] and status == "LONGER_IN_SAMPLE":
        note += (
            f" Codon structure VALID "
            f"(start: {codon_info['start_codon']}, stop: {codon_info['stop_codon']})."
        )

    return {
        "status":            status,
        "recon_length_nt":   recon_len,
        "ref_length_nt":     ref_length,
        "diff_nt":           diff_nt,
        "diff_pct":          diff_pct,
        "captured_5_nt":     captured_5_nt,
        "captured_3_nt":     captured_3_nt,
        "codon_valid":       codon_info["is_valid_orf"],
        "is_codon_multiple": codon_info["is_codon_multiple"],
        "has_start_codon":   codon_info["has_start_codon"],
        "has_stop_codon":    codon_info["has_stop_codon"],
        "start_codon":       codon_info["start_codon"],
        "stop_codon":        codon_info["stop_codon"],
        "note":              note,
    }


def _blast_aligned_to_cigar(qseq_aln: str, sseq_aln: str) -> str:
    """
    Derive a CIGAR string from a pair of gapped BLAST aligned sequences.
      gap in qseq  → deletion in query  (reference advances) → D
      gap in sseq  → insertion in query (query advances)      → I
      both bases   → match/mismatch                           → M
    """
    ops: List[str] = []
    cur_op: Optional[str] = None
    cur_len = 0
    for q, s in zip(qseq_aln, sseq_aln):
        if q == "-":
            op = "D"
        elif s == "-":
            op = "I"
        else:
            op = "M"
        if op == cur_op:
            cur_len += 1
        else:
            if cur_op:
                ops.append(f"{cur_len}{cur_op}")
            cur_op, cur_len = op, 1
    if cur_op:
        ops.append(f"{cur_len}{cur_op}")
    return "".join(ops)


# ── standard BLAST/DIAMOND tabular column names ──────────────────────────────
_BLAST_DEFAULT_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart",  "qend",  "sstart", "send",   "evalue",  "bitscore",
]


def _parse_blast_tabular(
    tsv_path: str,
    gene_name: str,
    query_seqs: dict[str, str],
    field_names: Optional[list[str]] = None,
) -> tuple[int, dict[str, str], list[dict]]:
    """
    Parse a BLAST / DIAMOND tabular output file (format 6 or custom).

    *tsv_path*   — path to the tabular file (plain or gzip).
    *gene_name*  — subject sequence ID to extract hits for.
    *query_seqs* — dict of {qseqid: sequence} loaded from the query FASTA.
                   May be empty if the tabular file contains qseq/sseq columns.
    *field_names* — ordered list of column names.  Defaults to the standard
                   12-column BLAST format 6 list.  Commonly useful extras:
                   ``qseq``, ``sseq`` (aligned sequences), ``slen``
                   (subject length).

    Returns ``(gene_length, read_sequences, alignments)``.
    """
    fields = field_names if field_names else _BLAST_DEFAULT_FIELDS

    def _col(name: str) -> Optional[int]:
        return fields.index(name) if name in fields else None

    i_qseqid = _col("qseqid")
    i_sseqid = _col("sseqid")
    i_qstart = _col("qstart")
    i_qend   = _col("qend")
    i_sstart = _col("sstart")
    i_send   = _col("send")
    i_qseq   = _col("qseq")   # optional — gapped aligned query seq
    i_sseq   = _col("sseq")   # optional — gapped aligned subject seq
    i_slen   = _col("slen")   # optional — subject (gene) length

    gene_length = 0
    read_seqs: Dict[str, str] = {}
    alignments: List[dict] = []

    opener = gzip.open if tsv_path.endswith((".gz", ".gzip")) else open
    with opener(tsv_path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # skip any header row that starts with a recognised field name
            if line.split("\t")[0] in ("qseqid", "query"):
                continue

            cols = line.split("\t")
            if len(cols) < 10:
                continue

            sseqid = cols[i_sseqid] if i_sseqid is not None and i_sseqid < len(cols) else None
            if sseqid != gene_name:
                continue

            qseqid = cols[i_qseqid] if i_qseqid is not None and i_qseqid < len(cols) else cols[0]

            try:
                qstart = int(cols[i_qstart]) if i_qstart is not None else 1
                qend   = int(cols[i_qend])   if i_qend   is not None else 1
                sstart = int(cols[i_sstart]) if i_sstart is not None else 1
                send   = int(cols[i_send])   if i_send   is not None else 1
            except (ValueError, IndexError):
                continue

            # Update gene length estimate
            if i_slen is not None and i_slen < len(cols):
                try:
                    gene_length = max(gene_length, int(cols[i_slen]))
                except ValueError:
                    pass
            gene_length = max(gene_length, max(sstart, send))

            # Reverse-strand: normalise so sstart < send
            reverse_hit = sstart > send
            if reverse_hit:
                sstart, send = send, sstart
                qstart, qend = qend, qstart

            # ── derive sequence + CIGAR ───────────────────────────────────
            if (
                i_qseq is not None and i_qseq < len(cols)
                and i_sseq is not None and i_sseq < len(cols)
                and cols[i_qseq] not in ("", "N/A", "*", "-")
            ):
                # Prefer pre-aligned sequences from the tabular file.
                # BLAST's qseq/sseq values already reflect the aligned hit
                # orientation, so reverse-strand rows must not be reversed again.
                qseq_aln = cols[i_qseq].upper()
                sseq_aln = cols[i_sseq].upper()
                seq = qseq_aln.replace("-", "")
                cigar = _blast_aligned_to_cigar(qseq_aln, sseq_aln)
            else:
                # Fall back to extracting the aligned slice from the query FASTA
                full_seq = query_seqs.get(qseqid, "")
                if not full_seq:
                    logger.warning(
                        f"  No sequence for query '{qseqid}' in query FASTA — "
                        "skipping hit.  Provide -query_fasta or use a tabular "
                        "format that includes qseq/sseq columns."
                    )
                    continue
                qs = min(qstart, qend) - 1  # 0-based
                qe = max(qstart, qend)       # exclusive
                seq = full_seq[qs:qe].upper()
                if reverse_hit:
                    seq = _reverse_complement(seq)
                # Approximate CIGAR — no gap positional info without qseq/sseq
                cigar = f"{len(seq)}M"

            if not seq:
                continue

            read_seqs[qseqid] = seq
            alignments.append({
                "name":  qseqid,
                "flag":  0,
                "pos":   sstart - 1,  # convert to 0-based
                "cigar": cigar,
                "seq":   seq,
                "qual":  "*",
                "md":    None,        # not available from BLAST
                "nm":    None,
            })

    return gene_length, read_seqs, alignments


def _mean_std(values: List[float]) -> Tuple[float, float]:
    if not values:
        return 0.0, 0.0
    n = len(values)
    mean = sum(values) / n
    if n < 2:
        return mean, 0.0
    variance = sum((v - mean) ** 2 for v in values) / (n - 1)
    return mean, math.sqrt(variance)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Pileup helpers (reusable for both full and per-haplotype builds)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _build_pileup(alignments: List[dict], gene_length: int):
    """
    base_pileup[ref_pos] = {base: count}
    ins_pileup[ref_pos_before_ins] = {ins_seq: count}
    """
    base_pileup: Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    ins_pileup:  Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for aln in alignments:
        seq   = aln["seq"]
        cigar = aln["cigar"]
        pos   = aln["pos"]
        if not seq or seq == "*" or not cigar or cigar == "*":
            continue
        ref_pos = pos
        read_pos = 0
        last_ref = pos

        for count_str, op in CIGAR_RE.findall(cigar):
            length = int(count_str)
            if op in ("M", "=", "X"):
                for k in range(length):
                    rp = read_pos + k
                    if rp < len(seq):
                        b = seq[rp].upper()
                        if b not in "ACGTN":
                            b = "N"
                        base_pileup[ref_pos + k][b] += 1
                last_ref = ref_pos + length - 1
                ref_pos  += length
                read_pos += length
            elif op == "I":
                ins_seq = seq[read_pos:read_pos + length].upper()
                ins_pileup[last_ref][ins_seq] += 1
                read_pos += length
            elif op == "D":
                for k in range(length):
                    base_pileup[ref_pos + k]["-"] += 1
                last_ref = ref_pos + length - 1
                ref_pos += length
            elif op == "N":
                ref_pos += length
            elif op == "S":
                read_pos += length

    return base_pileup, ins_pileup


def _call_consensus(
    base_pileup: Dict[int, Dict[str, int]],
    ins_pileup:  Dict[int, Dict[str, int]],
    gene_length: int,
    min_depth: int,
    min_freq: float,
    min_ins_depth: int,
    min_ins_freq: float,
) -> Tuple[str, str, List[dict]]:
    """
    Returns (ref_coord_consensus, sample_seq_with_ins, per_pos_stats).
    """
    consensus_chars = []
    per_pos: List[dict] = []

    for i in range(gene_length):
        votes = base_pileup.get(i, {})
        real = {b: c for b, c in votes.items() if b in "ACGT"}
        depth = sum(real.values())

        # insertion after this position?
        ins_info = ins_pileup.get(i, {})
        best_ins = max(ins_info, key=ins_info.get) if ins_info else None
        best_ins_count = ins_info[best_ins] if best_ins else 0
        ins_freq = best_ins_count / depth if (best_ins and depth > 0) else 0.0
        has_ins = (
            best_ins is not None
            and best_ins_count >= min_ins_depth
            and ins_freq >= min_ins_freq
        )

        if depth < min_depth:
            consensus_chars.append("N")
            per_pos.append({
                "pos": i + 1, "depth": depth,
                "consensus_base": "N", "top_base_freq": 0.0,
                "alt_base": "", "alt_base_freq": 0.0,
                "has_insertion": has_ins,
                "ins_seq": best_ins if has_ins else "",
                "ins_depth": best_ins_count if has_ins else 0,
                "ins_freq": round(ins_freq * 100, 2) if has_ins else 0.0,
            })
            continue

        sorted_votes = sorted(real.items(), key=lambda x: -x[1])
        top_base, top_count = sorted_votes[0]
        top_freq = top_count / depth

        alt_base = alt_freq = ""
        if len(sorted_votes) > 1:
            alt_base, alt_count = sorted_votes[1]
            alt_freq = round(alt_count / depth * 100, 2)

        if top_freq >= min_freq:
            cons_base = top_base
        else:
            cons_base = _iupac_ambiguity(real, depth)

        consensus_chars.append(cons_base)
        per_pos.append({
            "pos": i + 1, "depth": depth,
            "consensus_base": cons_base,
            "top_base_freq": round(top_freq * 100, 2),
            "alt_base": alt_base,
            "alt_base_freq": alt_freq,
            "has_insertion": has_ins,
            "ins_seq": best_ins if has_ins else "",
            "ins_depth": best_ins_count if has_ins else 0,
            "ins_freq": round(ins_freq * 100, 2) if has_ins else 0.0,
        })

    # build sample_seq inserting common insertions
    sample_parts = []
    for base, pinfo in zip(consensus_chars, per_pos):
        sample_parts.append(base)
        if pinfo["has_insertion"]:
            sample_parts.append(pinfo["ins_seq"])

    return "".join(consensus_chars), "".join(sample_parts), per_pos


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Multi-version detection
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _detect_multiple_versions(
    base_pileup: Dict[int, Dict[str, int]],
    gene_length: int,
    min_depth: int = 3,
    min_bimodal_freq: float = 0.20,
    min_informative_sites: int = 3,
) -> Tuple[str, List[dict]]:
    """
    Analyse the allele-frequency spectrum at every covered position.
    Returns (classification, informative_positions).
    """
    informative = []

    for i in range(gene_length):
        votes = base_pileup.get(i, {})
        real  = {b: c for b, c in votes.items() if b in "ACGT"}
        depth = sum(real.values())
        if depth < min_depth or len(real) < 2:
            continue
        sv = sorted(real.items(), key=lambda x: -x[1])
        top_base,  top_count  = sv[0]
        alt_base,  alt_count  = sv[1]
        alt_freq = alt_count / depth
        if alt_freq >= min_bimodal_freq:
            informative.append({
                "pos":      i,
                "top_base": top_base,
                "top_freq": round(top_count / depth, 4),
                "alt_base": alt_base,
                "alt_freq": round(alt_freq, 4),
                "depth":    depth,
            })

    if len(informative) < min_informative_sites:
        return "single", informative

    alt_freqs = [p["alt_freq"] for p in informative]
    mean_alt, _ = _mean_std(alt_freqs)

    if mean_alt >= 0.30:
        classification = "multi"
    elif mean_alt >= 0.15:
        classification = "uncertain"
    else:
        classification = "single"

    return classification, informative


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Haplotype separation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _extract_read_alleles(
    alignments: List[dict],
    informative_positions: List[dict],
) -> Dict[str, Dict[int, str]]:
    """
    For each read, recover which base it carries at each informative position.
    Returns {read_name: {ref_pos(0-based): base}}.
    Only positions in informative_positions are tracked.
    """
    pos_set: Set[int] = {p["pos"] for p in informative_positions}
    read_alleles: Dict[str, Dict[int, str]] = {}

    for aln in alignments:
        seq   = aln["seq"]
        cigar = aln["cigar"]
        pos   = aln["pos"]
        if not seq or seq == "*" or not cigar or cigar == "*":
            continue

        ref_pos  = pos
        read_pos = 0
        local_map: Dict[int, str] = {}

        for count_str, op in CIGAR_RE.findall(cigar):
            length = int(count_str)
            if op in ("M", "=", "X"):
                for k in range(length):
                    rp = ref_pos + k
                    if rp in pos_set:
                        idx = read_pos + k
                        if idx < len(seq):
                            local_map[rp] = seq[idx].upper()
                ref_pos  += length
                read_pos += length
            elif op == "I":
                read_pos += length
            elif op in ("D", "N"):
                ref_pos += length
            elif op == "S":
                read_pos += length

        if local_map:
            read_alleles[aln["name"]] = local_map

    return read_alleles


def _assign_reads_to_haplotypes(
    read_alleles: Dict[str, Dict[int, str]],
    informative_positions: List[dict],
    min_consistency: float = 0.70,
) -> Tuple[List[str], List[str], List[str]]:
    """
    Score each read against two profiles:
      hap1 profile = majority allele at every informative position
      hap2 profile = minority allele at every informative position

    Returns (hap1_reads, hap2_reads, ambiguous_reads).
    """
    pos_alleles: Dict[int, Tuple[str, str]] = {
        p["pos"]: (p["top_base"], p["alt_base"]) for p in informative_positions
    }

    hap1: List[str] = []
    hap2: List[str] = []
    ambiguous: List[str] = []

    for read_name, alleles in read_alleles.items():
        h1 = h2 = 0
        for rp, base in alleles.items():
            if rp not in pos_alleles:
                continue
            top_b, alt_b = pos_alleles[rp]
            if base == top_b:
                h1 += 1
            elif base == alt_b:
                h2 += 1
        total = h1 + h2
        if total == 0:
            ambiguous.append(read_name)
        elif h1 / total >= min_consistency:
            hap1.append(read_name)
        elif h2 / total >= min_consistency:
            hap2.append(read_name)
        else:
            ambiguous.append(read_name)

    return hap1, hap2, ambiguous


def _estimate_copy_ratio(informative_positions: List[dict]) -> str:
    """
    From the mean minor-allele frequencies, estimate the copy number ratio
    (e.g. '1:1', '2:1', '3:1') for the two detected haplotypes.
    """
    if not informative_positions:
        return "unknown"
    alt_freqs = [p["alt_freq"] for p in informative_positions]
    mean_alt, _ = _mean_std(alt_freqs)
    # major:minor = (1-mean_alt) : mean_alt
    # express as integer ratio
    if mean_alt <= 0:
        return "1:0"
    ratio = (1.0 - mean_alt) / mean_alt
    # round to nearest sensible integer ratio
    for major in range(1, 6):
        for minor in range(1, 6):
            if abs(ratio - major / minor) < 0.15:
                return f"{major}:{minor}"
    return f"{1.0 - mean_alt:.2f}:{mean_alt:.2f}"


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Validation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _analyze_divergence(
    consensus: str,
    reference: str,
    per_pos: List[dict],
    variants: List[dict]
) -> Dict:
    """
    Analyze how sample sequence diverges from reference.
    Divergence is EXPECTED - represents real biological variation!
    """
    if not reference or not consensus:
        return {"status": "no_reference"}
    
    min_len = min(len(consensus), len(reference))
    snps = insertions = deletions = identical = ambiguous = 0
    divergence_positions = []
    
    for i in range(min_len):
        cons_base = consensus[i]
        ref_base = reference[i]
        
        if cons_base == "N":
            ambiguous += 1
        elif cons_base == ref_base:
            identical += 1
        else:
            snps += 1
            read_support = depth = None
            if i < len(per_pos):
                read_support = per_pos[i].get("top_base_freq", 0)
                depth = per_pos[i].get("depth", 0)
            
            divergence_positions.append({
                "position": i + 1,
                "reference": ref_base,
                "sample": cons_base,
                "type": "SNP",
                "read_support_pct": read_support,
                "depth": depth
            })
    
    if len(consensus) > len(reference):
        insertions = len(consensus) - len(reference)
    elif len(reference) > len(consensus):
        deletions = len(reference) - len(consensus)
    
    total_comparable = min_len - ambiguous
    identity_pct = (identical / total_comparable * 100) if total_comparable > 0 else 0
    divergence_pct = (snps / total_comparable * 100) if total_comparable > 0 else 0
    
    if divergence_pct < 1:
        divergence_class = "NEARLY_IDENTICAL"
        interpretation = "Sample nearly identical to reference (< 1% divergence)"
    elif divergence_pct < 5:
        divergence_class = "MINOR_VARIANT"
        interpretation = "Sample is a minor variant of reference (1-5% divergence)"
    elif divergence_pct < 15:
        divergence_class = "SIGNIFICANT_VARIANT"
        interpretation = "Sample significantly diverged from reference (5-15% divergence)"
    elif divergence_pct < 30:
        divergence_class = "MAJOR_VARIANT"
        interpretation = "Sample is a major variant or different allele (15-30% divergence)"
    else:
        divergence_class = "HIGHLY_DIVERGENT"
        interpretation = "Sample highly divergent - possible different species/HGT (>30% divergence)"
    
    return {
        "status": "analyzed",
        "identity_pct": round(identity_pct, 2),
        "divergence_pct": round(divergence_pct, 2),
        "snps": snps,
        "insertions": insertions,
        "deletions": deletions,
        "identical_positions": identical,
        "ambiguous_positions": ambiguous,
        "total_comparable": total_comparable,
        "divergence_class": divergence_class,
        "interpretation": interpretation,
        "divergent_positions": divergence_positions[:50],
        "total_divergent_positions": len(divergence_positions)
    }


def _write_divergence_report(
    report_path: Path,
    gene_name: str,
    consensus: str,
    reference: str,
    divergence: Dict,
    per_pos: List[dict]
):
    """Write detailed divergence analysis report."""
    with open(report_path, "w") as fh:
        fh.write("=" * 70 + "\n")
        fh.write("SAMPLE vs REFERENCE DIVERGENCE ANALYSIS\n")
        fh.write("=" * 70 + "\n\n")
        
        fh.write("⚠  IMPORTANT: Divergence from reference is EXPECTED and GOOD!\n\n")
        fh.write("This report shows how the reconstructed SAMPLE sequence differs\n")
        fh.write("from the DATABASE REFERENCE. These differences represent REAL BIOLOGY:\n")
        fh.write("  • Strain/species variation\n")
        fh.write("  • Natural mutations\n")
        fh.write("  • Selective adaptation\n")
        fh.write("  • Allelic differences\n\n")
        fh.write("The reference is an ALIGNMENT SCAFFOLD - consensus from YOUR READS!\n\n")
        fh.write("=" * 70 + "\n\n")
        
        fh.write(f"Gene: {gene_name}\n")
        fh.write(f"Sample length: {len(consensus)} bp\n")
        fh.write(f"Reference length: {len(reference)} bp\n\n")
        
        if divergence["status"] != "analyzed":
            fh.write("No reference available for comparison.\n")
            return
        
        fh.write("-" * 70 + "\n")
        fh.write("DIVERGENCE SUMMARY\n")
        fh.write("-" * 70 + "\n\n")
        
        fh.write(f"Identity:        {divergence['identity_pct']:.2f}%\n")
        fh.write(f"Divergence:      {divergence['divergence_pct']:.2f}%\n\n")
        
        fh.write(f"Classification:  {divergence['divergence_class']}\n")
        fh.write(f"Interpretation:  {divergence['interpretation']}\n\n")
        
        fh.write(f"SNPs:            {divergence['snps']}\n")
        fh.write(f"Insertions:      {divergence['insertions']} bp\n")
        fh.write(f"Deletions:       {divergence['deletions']} bp\n")
        fh.write(f"Identical:       {divergence['identical_positions']}\n")
        fh.write(f"Ambiguous (N):   {divergence['ambiguous_positions']}\n\n")
        
        if divergence['divergent_positions']:
            fh.write("-" * 70 + "\n")
            fh.write("DIVERGENT POSITIONS (sample != reference)\n")
            fh.write("-" * 70 + "\n\n")
            fh.write("Pos\tRef\tSample\tType\tDepth\tSupport%\n")
            fh.write("-" * 70 + "\n")
            
            for div in divergence['divergent_positions']:
                pos = div['position']
                ref = div['reference']
                samp = div['sample']
                dtype = div['type']
                depth = div.get('depth', 'N/A')
                support = div.get('read_support_pct', 'N/A')
                fh.write(f"{pos}\t{ref}\t{samp}\t{dtype}\t{depth}\t{support}\n")
            
            if divergence['total_divergent_positions'] > 50:
                fh.write(f"\n... and {divergence['total_divergent_positions'] - 50} more\n")
            fh.write("\n")
        
        fh.write("-" * 70 + "\n")
        fh.write("METHODOLOGY\n")
        fh.write("-" * 70 + "\n\n")
        fh.write("Sample sequence reconstruction:\n")
        fh.write("  1. Reference used as ALIGNMENT TARGET\n")
        fh.write("  2. Reads aligned to reference coordinates\n")
        fh.write("  3. Consensus called from READ PILEUP (NOT reference!)\n")
        fh.write("  4. Most common read base = sample consensus\n\n")
        fh.write("All differences are supported by actual sample reads.\n")
        fh.write("This is the TRUE sample sequence, not reference-forced!\n\n")


def _validate_reconstruction(
    consensus: str,
    per_pos: List[dict],
    reference: Optional[str],
    base_pileup: Dict[int, Dict[str, int]],
    gene_length: int,
    min_depth: int,
    mv_classification: str,
    informative_positions: List[dict],
    haplotype_label: str = "",
    length_assessment: Optional[dict] = None,
) -> dict:
    """
    Compute validation metrics and assign an overall grade.

    When *length_assessment* indicates the sample gene is genuinely shorter
    than the reference but still codon-valid, the coverage denominator is
    adjusted to the reconstructed length so the grade reflects quality over
    the actual covered span rather than the uncovered reference tail.

    Returns a dict with all metrics and a 'grade' key.
    """
    depths = [p["depth"] for p in per_pos]
    covered = sum(1 for d in depths if d >= min_depth)
    cov_pct = covered / gene_length * 100 if gene_length > 0 else 0.0

    # ── Coverage adjustment for shorter-but-codon-valid genes ────────────────
    # The user explicitly selected these genes for reconstruction and accepts
    # that the sample version may be genuinely shorter than the reference.
    # When codon structure is intact, don't penalise coverage against the full
    # reference length — use the reconstructed span instead.
    if (length_assessment is not None
            and length_assessment.get("status") == "SHORTER_IN_SAMPLE"
            and length_assessment.get("codon_valid", False)):
        recon_len = length_assessment.get("recon_length_nt", gene_length)
        if 0 < recon_len < gene_length:
            recon_covered = sum(
                1 for p in per_pos
                if p["depth"] >= min_depth and p["pos"] <= recon_len
            )
            if recon_len > 0:
                cov_pct = recon_covered / recon_len * 100
                covered = recon_covered

    mean_depth, std_depth = _mean_std([float(d) for d in depths])
    cv_depth = (std_depth / mean_depth * 100) if mean_depth > 0 else 0.0

    n_pct = consensus.count("N") / max(len(consensus), 1) * 100

    # Read support score: at called (non-N) positions, what fraction of reads
    # agree with the consensus base?
    support_num = support_den = 0
    for pinfo in per_pos:
        cb = pinfo["consensus_base"]
        if cb == "N":
            continue
        votes = base_pileup.get(pinfo["pos"] - 1, {})
        real  = {b: c for b, c in votes.items() if b in "ACGT"}
        depth = sum(real.values())
        if depth == 0:
            continue
        support_num += real.get(cb, 0)
        support_den += depth
    read_support_pct = support_num / support_den * 100 if support_den > 0 else 0.0

    # Pairwise identity vs reference
    ref_identity = None
    if reference:
        ref_identity = _pct_identity(consensus, reference)

    # Allele frequency score (mean minor-allele freq at informative sites)
    alt_freqs = [p["alt_freq"] for p in informative_positions]
    mean_alt_freq, _ = _mean_std(alt_freqs) if alt_freqs else (0.0, 0.0)

    # ── Grade ────────────────────────────────────────────────────────────────
    # A: >95% coverage, >95% read support, CV < 80%, single version
    # B: >80% coverage, >85% read support
    # C: >60% coverage or mixed signal
    # F: <60% coverage or very low support
    if cov_pct >= 95 and read_support_pct >= 95 and cv_depth < 80 and n_pct < 5:
        grade = "A"
    elif cov_pct >= 80 and read_support_pct >= 85:
        grade = "B"
    elif cov_pct >= 60 and read_support_pct >= 70:
        grade = "C"
    else:
        grade = "F"

    # Downgrade if mixed signal and this is the blended consensus
    if mv_classification == "multi" and not haplotype_label:
        if grade in ("A", "B"):
            grade = "C*"   # * = multiple versions present, consensus is blended

    interpretation = {
        "A":  "High-confidence reconstruction.  Single version detected.",
        "B":  "Good reconstruction, minor gaps or low depth in some regions.",
        "C":  "Moderate confidence.  Check depth and coverage uniformity.",
        "C*": "Multiple versions detected — this is the blended consensus.  "
              "See haplotype FASTA for individual versions.",
        "F":  "Low confidence.  Insufficient coverage for reliable reconstruction.",
    }.get(grade, "")

    return {
        "haplotype_label":  haplotype_label or "consensus",
        "gene_length":      gene_length,
        "covered_positions": covered,
        "coverage_pct":     round(cov_pct, 2),
        "mean_depth":       round(mean_depth, 1),
        "cv_depth_pct":     round(cv_depth, 1),
        "n_pct":            round(n_pct, 2),
        "read_support_pct": round(read_support_pct, 2),
        "ref_identity_pct": round(ref_identity, 2) if ref_identity is not None else None,
        "mv_classification": mv_classification,
        "n_informative_sites": len(informative_positions),
        "mean_alt_freq_pct": round(mean_alt_freq * 100, 2),
        "grade":            grade,
        "interpretation":   interpretation,
        "length_assessment": length_assessment,
    }


def _write_validation_report(path: Path, validations: List[dict], gene_name: str):
    with open(path, "w") as fh:
        fh.write(f"RECONSTRUCTION VALIDATION REPORT\n")
        fh.write(f"Gene: {gene_name}\n")
        fh.write(f"Generated: {datetime.now().isoformat()}\n")
        fh.write("=" * 70 + "\n\n")

        for v in validations:
            fh.write(f"  Sequence: {v['haplotype_label']}\n")
            fh.write(f"  Grade:    {v['grade']}  — {v['interpretation']}\n\n")
            fh.write(f"  Coverage\n")
            fh.write(f"    Gene length:            {v['gene_length']} bp\n")
            fh.write(f"    Covered positions:      {v['covered_positions']} ({v['coverage_pct']}%)\n")
            fh.write(f"    N positions:            {v['n_pct']}%\n")
            fh.write(f"  Depth\n")
            fh.write(f"    Mean depth:             {v['mean_depth']}x\n")
            fh.write(f"    Depth CV:               {v['cv_depth_pct']}%  "
                     f"({'uniform' if v['cv_depth_pct'] < 50 else 'uneven'})\n")
            fh.write(f"  Read support\n")
            fh.write(f"    Support score:          {v['read_support_pct']}%\n")
            if v["ref_identity_pct"] is not None:
                fh.write(f"  Reference comparison\n")
                fh.write(f"    Identity vs reference:  {v['ref_identity_pct']}%\n")
            fh.write(f"  Multi-version analysis\n")
            fh.write(f"    Classification:         {v['mv_classification']}\n")
            fh.write(f"    Informative sites:      {v['n_informative_sites']}\n")
            fh.write(f"    Mean minor-allele freq: {v['mean_alt_freq_pct']}%\n")

            # ── Length-vs-reference assessment ────────────────────────────────
            la = v.get("length_assessment")
            if la:
                fh.write(f"  Length vs reference\n")
                fh.write(f"    Status:                 {la['status']}\n")
                fh.write(f"    Reconstructed length:   {la['recon_length_nt']} nt\n")
                fh.write(f"    Reference length:       {la['ref_length_nt']} nt\n")
                fh.write(f"    Difference:             {la['diff_nt']:+d} nt  "
                         f"({la['diff_pct']:+.1f}%)\n")
                fh.write(f"    Extension captured 5ʹ:  {la['captured_5_nt']} nt\n")
                fh.write(f"    Extension captured 3ʹ:  {la['captured_3_nt']} nt\n")
                fh.write(f"    Codon-valid ORF:        {'YES' if la['codon_valid'] else 'NO'}\n")
                if la["has_start_codon"] or la["has_stop_codon"]:
                    fh.write(
                        f"    Start / stop codons:    "
                        f"{la['start_codon'] or 'absent'} / "
                        f"{la['stop_codon']  or 'absent'}\n"
                    )
                fh.write(f"    Note: {la['note']}\n")

            fh.write("\n" + "-" * 70 + "\n\n")

        fh.write(
            "GRADE KEY\n"
            "  A   High-confidence, single clean version\n"
            "  B   Good, minor issues\n"
            "  C   Moderate confidence\n"
            "  C*  Multiple versions detected — blended consensus; use haplotypes\n"
            "  F   Insufficient coverage\n"
            "\n"
            "LENGTH STATUS\n"
            "  MATCHES_REFERENCE  — Reconstructed length equals reference length\n"
            "  SHORTER_IN_SAMPLE  — Sample gene shorter than reference\n"
            "                       (grade adjusted when codon structure is valid)\n"
            "  LONGER_IN_SAMPLE   — Sample gene extends beyond reference boundaries\n"
            "                       (captured via boundary extension of reads)\n"
        )


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Main reconstruction class
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class GeneReconstructor:

    def __init__(
        self,
        genefior_output_dir: str,
        gene_name: str,
        recon_dir: str,
        reference_fasta: Optional[str] = None,
        tool: Optional[str] = None,
        min_depth: int = 3,
        min_freq: float = 0.5,
        min_insertion_depth: int = 3,
        min_insertion_freq: float = 0.5,
        min_bimodal_freq: float = 0.20,
        min_informative_sites: int = 3,
        hap_min_reads: int = 5,
        max_hap_n_pct: float = 30.0,
        emit_reads_fasta: bool = True,
        emit_plots: bool = True,
        hap_on_uncertain: bool = True,
        boundary_extension: float = 0.20,
        codon_aware: bool = True,
        start_codons: Optional[str] = None,
        stop_codons: Optional[str] = None,
        # BLAST / DIAMOND mode
        blast_tsv: Optional[str] = None,
        query_fasta: Optional[str] = None,
        blast_format: Optional[list[str]] = None,
    ):
        # Validate parameters
        if min_depth < 1:
            raise ValueError(f"min_depth must be >= 1, got {min_depth}")
        if not (0.0 < min_freq <= 1.0):
           raise ValueError(f"min_freq must be in (0.0, 1.0], got {min_freq}")
        if not (0.0 < min_insertion_freq <= 1.0):
            raise ValueError(f"min_insertion_freq must be in (0.0, 1.0], got {min_insertion_freq}")
        if not (0.0 < min_bimodal_freq < 0.5):
            raise ValueError(f"min_bimodal_freq must be in (0.0, 0.5), got {min_bimodal_freq}")
        if min_informative_sites < 1:
            raise ValueError(f"min_informative_sites must be >= 1, got {min_informative_sites}")
        if hap_min_reads < 1:
            raise ValueError(f"hap_min_reads must be >= 1, got {hap_min_reads}")
        if not (0.0 <= max_hap_n_pct <= 100.0):
            raise ValueError(f"max_hap_n_pct must be in [0.0, 100.0], got {max_hap_n_pct}")
        
        self.genefior_dir        = Path(genefior_output_dir) if genefior_output_dir else Path(".")
        self.gene_name           = gene_name
        self.recon_dir           = Path(recon_dir)
        self.reference_fasta     = reference_fasta
        self.preferred_tool      = tool
        self.min_depth           = min_depth
        self.min_freq            = min_freq
        self.min_ins_depth       = min_insertion_depth
        self.min_ins_freq        = min_insertion_freq
        self.min_bimodal_freq    = min_bimodal_freq
        self.min_informative_sites = min_informative_sites
        self.hap_min_reads       = hap_min_reads
        self.max_hap_n_pct       = max_hap_n_pct
        self.emit_reads_fasta    = emit_reads_fasta
        self.emit_plots          = emit_plots
        self.hap_on_uncertain    = hap_on_uncertain
        self.boundary_extension  = max(0.0, boundary_extension)
        self.codon_aware         = codon_aware
        # Parse codon sets — accept either a comma-separated string or None
        self.start_codons: frozenset = (
            _parse_codon_list(start_codons, "start")
            if start_codons else _START_CODONS
        )
        self.stop_codons: frozenset = (
            _parse_codon_list(stop_codons, "stop")
            if stop_codons else _STOP_CODONS
        )
        # BLAST / DIAMOND
        self.blast_tsv           = blast_tsv
        self.query_fasta         = query_fasta
        self.blast_format        = blast_format   # list of column names or None

        self.raw_dir   = self.genefior_dir / "raw_outputs"
        self.stats_dir = self.genefior_dir / "tool_stats"

    # ── file discovery ───────────────────────────────────────────────────────

    def _find_alignment_files(self) -> List[Tuple[str, Path, Path]]:
        patterns = [
            ("Bowtie2",  "bowtie2_results_sorted.bam",  "bowtie2_results.sam"),
            ("BWA",      "bwa_results_sorted.bam",       "bwa_results.sam"),
            ("Minimap2", "minimap2_results_sorted.bam",  "minimap2_results.sam"),
        ]
        seen: Set[str] = set()
        candidates: List[Tuple[str, Path, Path]] = []

        for label, bam_sfx, sam_sfx in patterns:
            for bam in sorted(self.raw_dir.glob(f"*{bam_sfx}")):
                sam = bam.parent / bam.name.replace("_sorted.bam", ".sam").replace(
                    "_results_sorted", "_results"
                )
                key = str(sam)
                if key not in seen:
                    seen.add(key)
                    candidates.append((label, bam, sam))
            for sam in sorted(self.raw_dir.glob(f"*{sam_sfx}")):
                key = str(sam)
                if key not in seen:
                    bam = sam.parent / sam.name.replace("_results.sam", "_results_sorted.bam")
                    seen.add(key)
                    candidates.append((label, bam, sam))

        if self.preferred_tool:
            pref = self.preferred_tool.lower()
            filtered = [(t, b, s) for t, b, s in candidates if t.lower().startswith(pref)]
            if filtered:
                return filtered
            logger.warning(f"Tool '{self.preferred_tool}' not found; using all available.")
        return candidates

    def _find_stats_files(self, tool_label: str) -> List[Path]:
        """Return ALL stats TSVs for this tool (may come from multiple databases)."""
        return sorted(self.stats_dir.glob(f"*_{tool_label}_stats.tsv"))

    def _load_gene_stats(self, stats_paths) -> Optional[dict]:
        """Search one or more stats TSVs for the gene row, trying each in turn."""
        if not stats_paths:
            return None
        if isinstance(stats_paths, Path):
            stats_paths = [stats_paths]
        for stats_path in stats_paths:
            if not stats_path or not stats_path.exists():
                continue
            with open(stats_path) as fh:
                header = fh.readline().strip().split("\t")
                for line in fh:
                    row = dict(zip(header, line.strip().split("\t")))
                    if row.get("Gene") == self.gene_name:
                        return row
        return None

    # ── SAM parsing ──────────────────────────────────────────────────────────

    def _parse_sam_for_gene(self, sam_iter) -> Tuple[int, Dict[str, str], List[dict]]:
        gene_length = 0
        read_sequences: Dict[str, str] = {}
        alignments: List[dict] = []

        for line in sam_iter:
            if isinstance(line, bytes):
                line = line.decode()
            line = line.rstrip("\n")
            if line.startswith("@SQ"):
                for part in line.split("\t"):
                    if part.startswith("SN:") and part[3:] == self.gene_name:
                        pass  # found our gene, keep reading for LN
                parts = line.split("\t")
                sn = ln = None
                for p in parts:
                    if p.startswith("SN:"):
                        sn = p[3:]
                    elif p.startswith("LN:"):
                        try:
                            ln = int(p[3:])
                        except ValueError:
                            pass
                if sn == self.gene_name and ln:
                    gene_length = ln
                continue
            if line.startswith("@"):
                continue

            fields = line.split("\t")
            if len(fields) < 11:
                continue
            if fields[2] != self.gene_name:
                continue

            flag = int(fields[1])
            if flag & 0x4:            # unmapped
                continue
            if flag & 0x100 or flag & 0x800:  # secondary / supplementary
                continue

            rname = fields[0]
            if flag & 0x40:
                rname += "/1"
            elif flag & 0x80:
                rname += "/2"

            try:
                pos = int(fields[3]) - 1
            except ValueError:
                pos = 0

            seq  = fields[9]
            md   = nm = None
            for tag in fields[11:]:
                if tag.startswith("MD:Z:"):
                    md = tag[5:]
                elif tag.startswith("NM:i:"):
                    try:
                        nm = int(tag[5:])
                    except ValueError:
                        nm = 0

            if seq and seq != "*":
                read_sequences[rname] = seq

            alignments.append({
                "name": rname, "flag": flag,
                "pos": pos, "cigar": fields[5],
                "seq": seq, "qual": fields[10],
                "md": md, "nm": nm,
            })

        return gene_length, read_sequences, alignments

    # ── reference reconstruction from MD tags ────────────────────────────────

    def _reconstruct_reference(self, alignments: List[dict], gene_length: int) -> str:
        ref_vote: Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        for aln in alignments:
            seq, cigar, md, pos = aln["seq"], aln["cigar"], aln["md"], aln["pos"]
            if not seq or seq == "*" or not md or not cigar or cigar == "*":
                continue
            for rp, base in _reconstruct_ref_bases(seq, cigar, md, pos).items():
                ref_vote[rp][base] += 1
        return "".join(
            max(ref_vote[i], key=ref_vote[i].get) if ref_vote.get(i) else "N"
            for i in range(gene_length)
        )

    # ── variant calling ───────────────────────────────────────────────────────

    def _compute_variants(self, consensus: str, reference: str,
                          per_pos: List[dict]) -> List[dict]:
        variants = []
        for i in range(min(len(consensus), len(reference))):
            c, r = consensus[i], reference[i]
            if c == r:
                continue
            pi = per_pos[i] if i < len(per_pos) else {}
            depth = pi.get("depth", 0)
            vtype = (
                "low_coverage" if depth < self.min_depth
                else "ambiguous"  if c not in "ACGT"
                else "SNP"
            )
            variants.append({
                "pos": i + 1, "ref_base": r, "sample_base": c,
                "type": vtype, "depth": depth,
                "freq_pct": pi.get("top_base_freq", 0.0),
            })
        for pinfo in per_pos:
            if pinfo.get("has_insertion"):
                variants.append({
                    "pos": pinfo["pos"], "ref_base": "-",
                    "sample_base": pinfo["ins_seq"], "type": "insertion",
                    "depth": pinfo["ins_depth"], "freq_pct": pinfo["ins_freq"],
                })
        return sorted(variants, key=lambda v: v["pos"])

    # ── output writers ────────────────────────────────────────────────────────

    def _write_consensus_fasta(self, path: Path, sequences: List[tuple]):
        """sequences = [(header_suffix, sequence, stats_row[, length_assessment])]

        The optional 4th element is the ``_assess_length_vs_reference`` dict.
        When present, extra annotations are added to the FASTA header:
          [ref_len=N] [vs_ref=STATUS] [codon_valid=YES/NO] [ext_5=Xnt] [ext_3=Ynt]
        """
        with open(path, "w") as fh:
            for entry in sequences:
                label, seq, stats_row = entry[0], entry[1], entry[2]
                la = entry[3] if len(entry) > 3 else None
                n_pct = seq.count("N") / max(len(seq), 1) * 100
                reads  = stats_row.get("Num_Sequences_Mapped", "?") if stats_row else "?"
                cov    = stats_row.get("Gene_Coverage",        "?") if stats_row else "?"
                avg_id = stats_row.get("Avg_Identity",         "?") if stats_row else "?"

                # Build base header
                header = (
                    f">{self.gene_name}_{label} "
                    f"[len={len(seq)}]"
                )
                if la:
                    header += (
                        f" [ref_len={la['ref_length_nt']}]"
                        f" [vs_ref={la['status']}]"
                        f" [codon_valid={'YES' if la['codon_valid'] else 'NO'}]"
                        f" [ext_5={la['captured_5_nt']}nt]"
                        f" [ext_3={la['captured_3_nt']}nt]"
                    )
                header += (
                    f" [N_pct={n_pct:.1f}]"
                    f" [reads={reads}] [gene_cov={cov}%] [avg_id={avg_id}%]"
                )
                fh.write(header + "\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")

    def _write_depth_table(self, path: Path, per_pos: List[dict]):
        with open(path, "w") as fh:
            fh.write(
                "Position\tDepth\tConsensus_Base\tTop_Base_Freq_pct"
                "\tAlt_Base\tAlt_Base_Freq_pct"
                "\tHas_Insertion\tIns_Seq\tIns_Depth\tIns_Freq_pct\n"
            )
            for p in per_pos:
                fh.write(
                    f"{p['pos']}\t{p['depth']}\t{p['consensus_base']}"
                    f"\t{p['top_base_freq']}"
                    f"\t{p['alt_base']}\t{p['alt_base_freq']}"
                    f"\t{'Y' if p['has_insertion'] else 'N'}"
                    f"\t{p['ins_seq']}\t{p['ins_depth']}\t{p['ins_freq']}\n"
                )

    def _write_variants_tsv(self, path: Path, variants: List[dict]):
        with open(path, "w") as fh:
            fh.write("Position\tRef_Base\tSample_Base\tType\tDepth\tFreq_pct\n")
            for v in variants:
                fh.write(
                    f"{v['pos']}\t{v['ref_base']}\t{v['sample_base']}"
                    f"\t{v['type']}\t{v['depth']}\t{v['freq_pct']}\n"
                )

    def _write_reads_fasta(self, path: Path, read_seqs: Dict[str, str]):
        with open(path, "w") as fh:
            for name, seq in read_seqs.items():
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")

    # ── visualisation ─────────────────────────────────────────────────────────

    def _plot_reconstruction(
        self,
        tool_label: str,
        tool_dir: Path,
        gene_safe: str,
        per_pos: List[dict],
        variants: List[dict],
        informative_pos: List[dict],
        haplotypes: list,
        mv_class: str,
        val: dict,
        reference: Optional[str] = None,
        divergence: Optional[Dict] = None,
    ):
        """
        Write a self-contained SVG/HTML report visualising the reconstruction.

        Panel 1 — Coverage depth (bar chart, coloured by depth level)
        Panel 2 — Allele-frequency plot at variable positions
        Panel 3 — Sample vs Reference divergence track  (NEW!)
        Panel 4 — Gene-map annotation strip
        Panel 5 — Per-haplotype depth traces (only when haplotypes exist)

        Pure Python, no external dependencies.
        
        The divergence panel clearly distinguishes:
          GREEN  — sample matches reference
          RED    — sample SNP vs reference  
          ORANGE — informative (bimodal) position
          GREY   — N / insufficient depth
          PURPLE — insertion in sample
        """
        gene_length = val["gene_length"]
        if gene_length == 0:
            return

        has_divergence = (reference is not None and divergence is not None
                          and divergence.get("status") == "analyzed")

        # ── layout constants ──────────────────────────────────────────
        ML, MR, MT, MB = 65, 20, 10, 30  # margins (left/right/top/bottom per panel)
        W          = 1300
        PLOT_W     = W - ML - MR
        P1_H       = 200   # depth panel
        P2_H       = 130   # allele-freq panel
        P3_H       = 60    # divergence track (NEW)
        P4_H       = 44    # gene-map strip
        P5_H       = 130   # haplotype depth panel (conditional)
        TITLE_H    = 70
        GAP        = 10    # between panels
        LABEL_H    = 22

        scale = PLOT_W / gene_length  # px per bp

        def gx(pos_1based):
            """Map 1-based reference position to SVG x."""
            return ML + (pos_1based - 1) * scale

        # accumulated vertical positions
        panels = []  # (y_start, height, label)
        y = TITLE_H
        panels.append((y, P1_H, "depth")); y += P1_H + GAP + LABEL_H
        panels.append((y, P2_H, "allele")); y += P2_H + GAP + LABEL_H
        if has_divergence:
            panels.append((y, P3_H, "divergence")); y += P3_H + GAP + LABEL_H
        panels.append((y, P4_H, "genemap")); y += P4_H + GAP + LABEL_H
        if haplotypes:
            panels.append((y, P5_H, "haplotypes")); y += P5_H + GAP + LABEL_H
        TOTAL_H = y + MB

        # panel index helpers
        def _panel(label):
            for entry in panels:
                if entry[2] == label:
                    return entry[0], entry[1]
            return None, None

        lines: List[str] = []

        def e(tag):
            return f"</{tag}>"

        def rect(x, y, w, h, fill, opacity=1.0, rx=0, extra=""):
            return (f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(w,0.3):.2f}" '
                    f'height="{h:.2f}" fill="{fill}" opacity="{opacity}" rx="{rx}" {extra}/>')

        def text(x, y, s, size=11, anchor="start", weight="normal", fill="#222"):
            return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
                    f'font-family="monospace" text-anchor="{anchor}" '
                    f'font-weight="{weight}" fill="{fill}">{s}</text>')

        def line(x1, y1, x2, y2, stroke="#888", w=1, dash=""):
            d = f'stroke-dasharray="{dash}"' if dash else ""
            return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" stroke-width="{w}" {d}/>'

        # ── helpers for run-length encoding ────────────────────────────
        def _runs(items, key_fn):
            """Group consecutive items with same key into (key, start, count) tuples."""
            runs = []
            for item in items:
                k = key_fn(item)
                if runs and runs[-1][0] == k:
                    runs[-1] = (k, runs[-1][1], runs[-1][2] + 1)
                else:
                    runs.append([k, item["pos"], 1])
            return runs

        # ── depth colours ──────────────────────────────────────────────
        def _depth_col(d):
            if d < self.min_depth:  return "#d73027"
            if d < 10:              return "#fee08b"
            if d < 30:              return "#74c476"
            return                         "#238b45"

        # ── SVG header ────────────────────────────────────────────────
        lines.append(
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{W}" height="{TOTAL_H}" '
            f'viewBox="0 0 {W} {TOTAL_H}" '
            f'style="background:#fff;font-family:sans-serif">'
        )
        lines.append('<style>text{font-family:Arial,sans-serif}</style>')

        # ── title bar ─────────────────────────────────────────────────
        grade      = val["grade"]
        grade_cols = {"A": "#1a9850", "B": "#4dac26", "C": "#fdae61",
                      "C*": "#f46d43", "F": "#d73027"}
        gcol       = grade_cols.get(grade, "#555")
        lines.append(rect(0, 0, W, TITLE_H, "#1e3a5f"))
        lines.append(text(12, 22, f"{self.gene_name}   ·   {tool_label}", 15,
                          fill="#fff", weight="bold"))
        badge_x = 12
        div_badge_items = []
        if has_divergence:
            d_pct = divergence["divergence_pct"]
            d_cls = divergence["divergence_class"].replace("_", " ")
            div_badge_items = [
                (f"vs ref: {d_pct:.1f}% divergence ({d_cls})", None,
                 "#e05c2e" if d_pct > 5 else "#a8522a"),
            ]
        for label, _, col in [
            (f"Grade: {grade}", None, gcol),
            (f"Coverage: {val['coverage_pct']}%", None, "#4393c3"),
            (f"Mean depth: {val['mean_depth']}×", None, "#74c476"),
            (f"Read support: {val['read_support_pct']}%", None, "#9ecae1"),
            (f"Multi-version: {mv_class}", None, "#fdae61" if mv_class != "single" else "#aaa"),
        ] + div_badge_items:
            lines.append(rect(badge_x, 30, len(label) * 7.5 + 10, 16, col, rx=4))
            lines.append(text(badge_x + 5, 41.5, label, 9, fill="#fff", weight="bold"))
            badge_x += len(label) * 7.5 + 16

        # ── Panel 1: Coverage depth ────────────────────────────────────
        y1, h1 = _panel("depth")
        max_dep = max((p["depth"] for p in per_pos), default=1) or 1
        pad_dep = max_dep * 1.05
        inner_h = h1

        # background
        lines.append(rect(ML, y1, PLOT_W, h1, "#f8f8f8"))

        # depth bars (run-length encode by colour for speed)
        for p in per_pos:
            bh = p["depth"] / pad_dep * inner_h
            lines.append(rect(gx(p["pos"]), y1 + inner_h - bh,
                               scale, bh, _depth_col(p["depth"])))

        # N-runs shading (grey overlay)
        in_n = n_start = None
        for p in per_pos:
            if p["consensus_base"] == "N":
                if in_n is None:
                    in_n = True; n_start = p["pos"]
            else:
                if in_n is not None:
                    lines.append(rect(gx(n_start), y1, (p["pos"] - n_start) * scale,
                                      h1, "#969696", opacity=0.22))
                    in_n = None
        if in_n is not None:
            lines.append(rect(gx(n_start), y1, (gene_length + 1 - n_start) * scale,
                              h1, "#969696", opacity=0.22))

        # min-depth line
        min_y = y1 + inner_h - (self.min_depth / pad_dep * inner_h)
        lines.append(line(ML, min_y, ML + PLOT_W, min_y, "#d73027", 1, "4,3"))

        # y-axis ticks
        for pct in [0.25, 0.50, 0.75, 1.0]:
            ty = y1 + inner_h * (1 - pct)
            tv = int(pad_dep * pct)
            lines.append(line(ML - 4, ty, ML, ty, "#888", 1))
            lines.append(text(ML - 7, ty + 4, str(tv), 8, anchor="end", fill="#555"))
        lines.append(line(ML, y1, ML, y1 + h1, "#888", 1))

        # panel label
        panel_label_y = y1 + h1 + GAP + LABEL_H - 6
        lines.append(text(ML, panel_label_y, "Coverage depth", 10, fill="#333", weight="bold"))
        depth_legend = [
            ("#238b45", "≥30×"), ("#74c476", "10–29×"),
            ("#fee08b", f"{self.min_depth}–9×"), ("#d73027", f"<{self.min_depth}× (N)"),
        ]
        lx = ML + 130
        for col, lbl in depth_legend:
            lines.append(rect(lx, panel_label_y - 10, 12, 10, col))
            lines.append(text(lx + 14, panel_label_y, lbl, 8, fill="#444"))
            lx += len(lbl) * 6.5 + 24

        # ── Panel 2: Allele frequencies ────────────────────────────────
        y2, h2 = _panel("allele")
        informative_set = {p["pos"] for p in informative_pos}

        lines.append(rect(ML, y2, PLOT_W, h2, "#f8f8f8"))
        # 50% guide line
        lines.append(line(ML, y2 + h2 * 0.5, ML + PLOT_W, y2 + h2 * 0.5, "#ccc", 1))

        bar_w  = max(scale * 1.8, 2.0)
        for p in per_pos:
            af = p["alt_base_freq"]
            if not isinstance(af, (int, float)) or af < 1.0:
                continue
            tf   = p["top_base_freq"]
            bx   = gx(p["pos"]) - bar_w / 2
            # top-allele bar (blue)
            th   = tf / 100 * h2
            lines.append(rect(bx, y2 + h2 - th, bar_w, th, "#4393c3"))
            # alt-allele bar (red) on top
            ah   = af / 100 * h2
            lines.append(rect(bx, y2 + h2 - th - ah, bar_w, ah, "#d6604d"))
            # orange diamond for informative sites
            if p["pos"] in informative_set:
                cx = bx + bar_w / 2
                cy = y2 + 6
                s2 = 5
                lines.append(
                    f'<polygon points="{cx},{cy - s2} {cx + s2},{cy} '
                    f'{cx},{cy + s2} {cx - s2},{cy}" fill="#f97a00"/>'
                )

        # y-axis
        for pct in [0.25, 0.5, 0.75, 1.0]:
            ty = y2 + h2 * (1 - pct)
            lines.append(line(ML - 4, ty, ML, ty, "#888", 1))
            lines.append(text(ML - 7, ty + 4, f"{int(pct*100)}%", 8, anchor="end", fill="#555"))
        lines.append(line(ML, y2, ML, y2 + h2, "#888", 1))

        panel_label_y = y2 + h2 + GAP + LABEL_H - 6
        lines.append(text(ML, panel_label_y,
                          "Allele frequencies at variable positions  (◆ = informative site)",
                          10, fill="#333", weight="bold"))
        for col, lbl in [("#4393c3", "Major allele"), ("#d6604d", "Minor allele"),
                          ("#f97a00", "Informative site")]:
            lines.append(rect(lx := (ML + 350 + (ML + 350 < lx) * 0), panel_label_y - 10, 12, 10, col))

        lx = ML + 350
        for col, lbl in [("#4393c3", "Major allele"), ("#d6604d", "Minor allele"),
                          ("#f97a00", "Informative site")]:
            lines.append(rect(lx, panel_label_y - 10, 12, 10, col))
            lines.append(text(lx + 14, panel_label_y, lbl, 8, fill="#444"))
            lx += len(lbl) * 6.5 + 24

        # ── Panel 3: Sample vs Reference Divergence ────────────────────
        if has_divergence:
            yd, hd = _panel("divergence")
            ref_seq = reference  # full reference string

            # Build per-position state for divergence track
            # State: "match", "snp", "informative", "insertion", "N", "no_ref"
            div_pos_set = {d["position"] for d in divergence.get("divergent_positions", [])}

            # Two row tracks inside the panel:
            # Row A = Reference bases (top third)
            # Row B = Sample bases (bottom third)
            # Middle = colour-coded divergence bar
            row_pad   = 4
            bar_top   = yd + row_pad
            bar_h     = hd - row_pad * 2
            mid_y     = bar_top + bar_h / 2
            row_a_h   = bar_h * 0.35     # ref row
            row_b_h   = bar_h * 0.35     # sample row
            bar_strip = bar_h * 0.30     # colour band between them

            # Background
            lines.append(rect(ML, yd, PLOT_W, hd, "#f8f8f8"))

            # Horizontal separators
            lines.append(line(ML, bar_top + row_a_h, ML + PLOT_W,
                              bar_top + row_a_h, "#ddd", 0.5))
            lines.append(line(ML, bar_top + row_a_h + bar_strip,
                              ML + PLOT_W, bar_top + row_a_h + bar_strip, "#ddd", 0.5))

            # Colour-coded strip RLE by divergence state
            def _div_state(p):
                pos = p["pos"]
                cb  = p["consensus_base"]
                if cb == "N":
                    return "N"
                if pos in informative_set:
                    return "informative"
                if pos in div_pos_set:
                    return "snp"
                return "match"

            _DIV_COLORS = {
                "match":       "#4dac26",
                "snp":         "#d73027",
                "informative": "#f97a00",
                "insertion":   "#9c6fc7",
                "N":           "#bababa",
            }

            div_runs = _runs(per_pos, _div_state)
            strip_y  = bar_top + row_a_h
            for state, start_pos, count in div_runs:
                col = _DIV_COLORS.get(state, "#aaa")
                lines.append(rect(gx(start_pos), strip_y,
                                  count * scale, bar_strip, col))

            # Only render individual base labels if scale >= 6px/bp
            if scale >= 6:
                for p in per_pos:
                    px    = gx(p["pos"]) + scale * 0.5
                    rbase = ref_seq[p["pos"] - 1] if p["pos"] - 1 < len(ref_seq) else "?"
                    sbase = p["consensus_base"]
                    # Reference row
                    rfill = "#555"
                    lines.append(text(px, bar_top + row_a_h - 2, rbase, 7,
                                      anchor="middle", fill=rfill))
                    # Sample row
                    sfill = "#d73027" if sbase != rbase and sbase != "N" else (
                            "#aaa" if sbase == "N" else "#1a7a1a")
                    lines.append(text(px, bar_top + row_a_h + bar_strip + row_b_h - 1,
                                      sbase, 7, anchor="middle", fill=sfill))
            else:
                # At small scale: tooltip rectangles on SNP positions
                for p in per_pos:
                    if p["pos"] not in div_pos_set:
                        continue
                    rbase = ref_seq[p["pos"] - 1] if p["pos"] - 1 < len(ref_seq) else "?"
                    sbase = p["consensus_base"]
                    tip = f"{p['pos']}: ref={rbase} → sample={sbase}"
                    px = gx(p["pos"])
                    lines.append(
                        f'<rect x="{px:.2f}" y="{bar_top:.1f}" '
                        f'width="{max(scale,1):.2f}" height="{bar_h:.1f}" '
                        f'fill="transparent">'
                        f'<title>{tip}</title></rect>'
                    )

            # Row labels on left axis
            lines.append(text(ML - 5, bar_top + row_a_h - 2,
                               "Ref", 7, anchor="end", fill="#555"))
            lines.append(text(ML - 5, bar_top + row_a_h + bar_strip + row_b_h - 1,
                               "Sample", 7, anchor="end", fill="#555"))

            # y-axis bar
            lines.append(line(ML, bar_top, ML, bar_top + bar_h, "#888", 1))

            # Panel label + statistics
            div_pct  = divergence["divergence_pct"]
            div_cls  = divergence["divergence_class"].replace("_", " ")
            n_snps   = divergence["snps"]
            panel_label_y = yd + hd + GAP + LABEL_H - 6
            lines.append(text(
                ML, panel_label_y,
                f"Sample vs Reference divergence  —  "
                f"{div_pct:.1f}% divergence ({n_snps} SNPs)  ·  {div_cls}"
                f"  ⚠ Note: divergence is EXPECTED BIOLOGY, not an error.",
                9, fill="#333", weight="bold"
            ))
            lx = ML + 350
            for col, lbl in [
                ("#4dac26", "Matches ref"),
                ("#d73027", "SNP vs ref"),
                ("#f97a00", "Informative SNP"),
                ("#9c6fc7", "Insertion"),
                ("#bababa", "N (insufficient depth)")
            ]:
                lines.append(rect(lx, panel_label_y - 10, 12, 10, col))
                lines.append(text(lx + 14, panel_label_y, lbl, 8, fill="#444"))
                lx += len(lbl) * 6.5 + 24

        # ── Panel 4: Gene map ──────────────────────────────────────────
        y3, h3 = _panel("genemap")
        snp_set = {v["pos"] for v in variants if v["type"] == "SNP"}
        ins_set = {v["pos"] for v in variants if v["type"] == "insertion"}


        def _map_col(p):
            if p["consensus_base"] == "N": return "#bababa"
            if p["pos"] in informative_set:  return "#f97a00"
            if p["pos"] in snp_set:          return "#d73027"
            if p["pos"] in ins_set:          return "#fe8e59"
            return "#4dac26"

        # run-length encode gene map
        runs = _runs(per_pos, _map_col)
        for col, start_pos, count in runs:
            lines.append(rect(gx(start_pos), y3 + 2, count * scale, h3 - 4, col))

        # x-axis ticks for gene map
        tick_step = max(1, (gene_length // 10) // 100 * 100 or 100)
        for tick_pos in range(1, gene_length + 1, tick_step):
            tx = gx(tick_pos)
            lines.append(line(tx, y3 + h3 + 2, tx, y3 + h3 + 8, "#888", 1))
            lines.append(text(tx, y3 + h3 + 18, str(tick_pos), 8, anchor="middle", fill="#555"))

        panel_label_y = y3 + h3 + GAP + LABEL_H - 6
        lines.append(text(ML, panel_label_y, "Gene map", 10, fill="#333", weight="bold"))
        lx = ML + 90
        for col, lbl in [("#4dac26", "Matches ref"), ("#d73027", "SNP"),
                          ("#f97a00", "Informative"), ("#fe8e59", "Insertion"),
                          ("#bababa", "Uncovered")]:
            lines.append(rect(lx, panel_label_y - 10, 12, 10, col))
            lines.append(text(lx + 14, panel_label_y, lbl, 8, fill="#444"))
            lx += len(lbl) * 6.5 + 24

        # ── Panel 5: Haplotype depths (optional) ──────────────────────
        yh, hh = _panel("haplotypes")
        if yh is not None and haplotypes:
            hap_colors = ["#2171b5", "#cb181d", "#238b45", "#984ea3"]

            lines.append(rect(ML, yh, PLOT_W, hh, "#f8f8f8"))

            max_hap_dep = max(
                (p["depth"] for hap in haplotypes for _, _, _, pp, _ in [hap] for p in pp),
                default=1
            ) or 1
            pad_hap = max_hap_dep * 1.05

            for hi, (hap_label_h, h_cons, h_sample, h_pp, h_val) in enumerate(haplotypes):
                col = hap_colors[hi % len(hap_colors)]
                pts_str = " ".join(
                    f"{gx(p['pos']):.1f},{yh + hh - p['depth'] / pad_hap * hh:.1f}"
                    for p in h_pp
                )
                if pts_str:
                    lines.append(
                        f'<polyline points="{pts_str}" fill="none" '
                        f'stroke="{col}" stroke-width="1.2" opacity="0.85"/>'
                    )
                # legend badge
                bx2 = ML + hi * 280
                lines.append(rect(bx2, yh + hh + GAP - 2, 12, 10, col))
                lines.append(text(bx2 + 15, yh + hh + GAP + 8,
                                  f"{hap_label_h}  grade={h_val['grade']}  "
                                  f"cov={h_val['coverage_pct']}%",
                                  9, fill="#333"))

            # y-axis
            for pct in [0.5, 1.0]:
                ty = yh + hh * (1 - pct)
                lines.append(line(ML - 4, ty, ML, ty, "#888", 1))
                lines.append(text(ML - 7, ty + 4, str(int(pad_hap * pct)), 8,
                                  anchor="end", fill="#555"))
            lines.append(line(ML, yh, ML, yh + hh, "#888", 1))
            lines.append(text(10, yh + hh // 2, "Depth", 9, fill="#555"))

            panel_label_y = yh + hh + GAP + LABEL_H - 6
            lines.append(text(ML, panel_label_y,
                               "Per-haplotype coverage depth", 10, fill="#333", weight="bold"))

        lines.append("</svg>")

        html = (
            "<!DOCTYPE html><html><head>"
            f'<meta charset="utf-8"><title>{self.gene_name} — {tool_label} reconstruction</title>'
            "</head><body style='margin:0;padding:10px;background:#eee'>"
            + "\n".join(lines)
            + "</body></html>"
        )

        plot_path = tool_dir / f"{gene_safe}_reconstruction_plot.html"
        try:
            plot_path.write_text(html, encoding="utf-8")
            logger.info(f"  Plot saved: {plot_path.name}  (open in any browser)")
        except Exception as exc:
            logger.warning(f"  Could not save plot: {exc}")

    # ── assembly graph generation ─────────────────────────────────────────────

    def _write_assembly_graph(
        self,
        tool_label: str,
        tool_dir: Path,
        gene_safe: str,
        consensus: str,
        per_pos: List[dict],
        variants: List[dict],
        informative_pos: List[dict],
        haplotypes: list,
        mv_class: str,
        base_pileup: Dict[int, Dict[str, int]],
        alignments: List[dict],
    ):
        """
        Generate GFA (Graphical Fragment Assembly) format file for Bandage visualization.

        Creates an assembly graph showing:
        - Segments for high-quality consensus regions
        - Variant bubbles at SNP/informative positions
        - Coverage information
        - Haplotype paths (when available)
        - Read alignments

        GFA format specification: https://github.com/GFA-spec/GFA-spec
        """
        gfa_path = tool_dir / f"{gene_safe}_assembly_graph.gfa"

        informative_set = {p["pos"] for p in informative_pos}
        variant_set = {v["pos"] for v in variants}

        lines = []
        lines.append("H\tVN:Z:1.0")
        lines.append(f"# GeneFior-Reconstruct {RECONSTRUCT_VERSION} Assembly Graph")
        lines.append(f"# Gene: {self.gene_name}")
        lines.append(f"# Tool: {tool_label}")
        lines.append(f"# Multi-version: {mv_class}")
        lines.append(f"# Consensus length: {len(consensus)} bp")

        # Strategy: Create segments for contiguous non-variant regions
        # Create variant bubbles at informative/variant positions
        # Link segments together

        segments = []  # (seg_id, start_pos, end_pos, sequence, avg_depth)
        edges = []     # (from_seg, to_seg, orientation)
        paths = {}     # {path_name: [(seg_id, orientation), ...]}

        # ── Build segments ────────────────────────────────────────────────
        seg_id = 1
        in_segment = False
        seg_start = None
        seg_seq = []
        seg_depths = []

        for i, p in enumerate(per_pos):
            pos = p["pos"]  # 1-based
            is_variant = pos in variant_set or pos in informative_set

            if not in_segment and not is_variant:
                # Start new segment
                in_segment = True
                seg_start = pos
                seg_seq = [p["consensus_base"]]
                seg_depths = [p["depth"]]

            elif in_segment and not is_variant:
                # Continue segment
                seg_seq.append(p["consensus_base"])
                seg_depths.append(p["depth"])

            elif in_segment and is_variant:
                # Close segment before variant
                if seg_seq:
                    avg_d = sum(seg_depths) / len(seg_depths) if seg_depths else 0
                    segments.append((
                        f"s{seg_id}",
                        seg_start,
                        pos - 1,
                        "".join(seg_seq),
                        avg_d
                    ))
                    seg_id += 1
                in_segment = False

                # Create variant bubble (two alternative segments)
                votes = base_pileup.get(pos - 1, {})  # 0-based for pileup
                real = {b: c for b, c in votes.items() if b in "ACGT"}
                if len(real) >= 2:
                    sorted_votes = sorted(real.items(), key=lambda x: -x[1])
                    # Major allele path
                    segments.append((
                        f"s{seg_id}",
                        pos,
                        pos,
                        sorted_votes[0][0],
                        sorted_votes[0][1]
                    ))
                    major_id = f"s{seg_id}"
                    seg_id += 1
                    # Minor allele path (if significant)
                    segments.append((
                        f"s{seg_id}",
                        pos,
                        pos,
                        sorted_votes[1][0],
                        sorted_votes[1][1]
                    ))
                    minor_id = f"s{seg_id}"
                    seg_id += 1
                else:
                    # Single allele at this position
                    segments.append((
                        f"s{seg_id}",
                        pos,
                        pos,
                        p["consensus_base"],
                        p["depth"]
                    ))
                    seg_id += 1

        # Close final segment
        if in_segment and seg_seq:
            avg_d = sum(seg_depths) / len(seg_depths) if seg_depths else 0
            segments.append((
                f"s{seg_id}",
                seg_start,
                per_pos[-1]["pos"],
                "".join(seg_seq),
                avg_d
            ))

        # ── Write segments (S lines) ──────────────────────────────────────
        for seg_id_str, start, end, seq, depth in segments:
            # GFA S line: S <seg_id> <sequence> [tags]
            # Tags: LN:i:length, DP:f:depth, SC:i:start_coord, EC:i:end_coord
            lines.append(
                f"S\t{seg_id_str}\t{seq}\t"
                f"LN:i:{len(seq)}\t"
                f"DP:f:{depth:.1f}\t"
                f"SC:i:{start}\t"
                f"EC:i:{end}"
            )

        # ── Create edges (L lines) ────────────────────────────────────────
        # Link consecutive segments
        for i in range(len(segments) - 1):
            seg1_id = segments[i][0]
            seg2_id = segments[i + 1][0]
            # GFA L line: L <from> <from_orient> <to> <to_orient> <overlap>
            # We use + orientation and 0 overlap (segments are adjacent, not overlapping)
            lines.append(f"L\t{seg1_id}\t+\t{seg2_id}\t+\t0M")

        # ── Create paths (P lines) ────────────────────────────────────────
        # Main consensus path (only main segments, not alternates)
        main_segments = [s for s in segments if "_alt" not in s[0]]
        consensus_path_segs = ",".join(f"{s[0]}+" for s in main_segments)
        lines.append(f"P\tconsensus\t{consensus_path_segs}\t*\t"
                   f"LN:i:{len(consensus)}\t"
                   f"FC:i:{len([p for p in per_pos if p['depth'] >= 1])}")

        # Haplotype paths (if available)
        if haplotypes:
            for hap_idx, (hap_label, h_cons, h_sample, h_pp, h_val) in enumerate(haplotypes, 1):
                # Create path with metadata
                mean_d = h_val.get('mean_depth', 0.0)
                lines.append(f"P\t{hap_label}\t{consensus_path_segs}\t*\t"
                           f"CV:f:{h_val['coverage_pct']}\t"
                           f"GR:Z:{h_val['grade']}\t"
                           f"DP:f:{mean_d}")

        # ── Add read alignments as containment records (A lines) ──────────
        # Sample up to 100 reads to avoid huge files
        sampled_alns = alignments[:100] if len(alignments) > 100 else alignments
        for aln in sampled_alns:
            read_name = aln["name"]
            # Find which segments this read overlaps
            # For simplicity, just annotate with position range
            pos = aln["pos"] + 1  # 1-based
            cigar = aln.get("cigar", "")
            if cigar and cigar != "*":
                # Estimate alignment length from CIGAR
                aln_len = sum(
                    int(m.group(1)) for m in CIGAR_RE.finditer(cigar)
                    if m.group(2) in "MDN="
                )
                lines.append(f"# Read: {read_name} pos:{pos}-{pos + aln_len}")

        # ── Write GFA file ────────────────────────────────────────────────
        try:
            with open(gfa_path, "w") as fh:
                fh.write("\n".join(lines) + "\n")
            logger.info(f"  Assembly graph saved: {gfa_path.name}  (view in Bandage)")
        except Exception as exc:
            logger.warning(f"  Could not save assembly graph: {exc}")

    # ── shared analysis pipeline ──────────────────────────────────────────────

    def _run_analysis_pipeline(
        self,
        tool_label: str,
        alignments: List[dict],
        read_seqs: Dict[str, str],
        gene_length: int,
        gene_stats_row: Optional[dict],
        used_ref: Optional[str],
        ref_source: str,
    ) -> Optional[dict]:
        """
        Given a set of alignments for one tool/source, run the full
        pileup → consensus → multi-version → haplotype → validation →
        output-write pipeline.  Returns the per-tool result dict or None
        if there was nothing to process.
        """
        gene_safe = _safe_filename(self.gene_name)

        # ── pileup ────────────────────────────────────────────────────
        base_pileup, ins_pileup = _build_pileup(alignments, gene_length)

        # ── single consensus ──────────────────────────────────────────
        consensus_ref, sample_seq, per_pos = _call_consensus(
            base_pileup, ins_pileup, gene_length,
            self.min_depth, self.min_freq,
            self.min_ins_depth, self.min_ins_freq,
        )
        covered = sum(1 for p in per_pos if p["depth"] >= self.min_depth)
        logger.info(
            f"  Consensus: {covered}/{gene_length} positions covered "
            f"({covered / gene_length * 100:.1f}%)"
        )

        # ── reference ─────────────────────────────────────────────────
        if not used_ref:
            used_ref   = self._reconstruct_reference(alignments, gene_length)
            ref_source = "MD-tag reconstruction"
            logger.info("  Reference reconstructed from MD tags.")

        # ── multi-version detection ───────────────────────────────────
        mv_class, informative_pos = _detect_multiple_versions(
            base_pileup, gene_length,
            min_depth=self.min_depth,
            min_bimodal_freq=self.min_bimodal_freq,
            min_informative_sites=self.min_informative_sites,
        )
        logger.info(
            f"  Multi-version classification: {mv_class.upper()} "
            f"({len(informative_pos)} informative sites)"
        )
        if mv_class != "single" and informative_pos:
            copy_ratio = _estimate_copy_ratio(informative_pos)
            logger.info(f"  Estimated copy ratio: {copy_ratio}")

        # ── haplotype separation ──────────────────────────────────────
        haplotypes = []

        if mv_class == "multi" and len(informative_pos) >= self.min_informative_sites:
            attempt_hap = True
            hap_note = "multi"
        elif mv_class == "uncertain" and self.hap_on_uncertain and len(informative_pos) >= self.min_informative_sites:
            attempt_hap = True
            hap_note = "uncertain"
        else:
            attempt_hap = False
            hap_note = ""

        if attempt_hap:
            if hap_note == "uncertain":
                logger.info(
                    "  Attempting haplotype separation on 'uncertain' classification "
                    "(use -no_hap_uncertain to disable)."
                )
            read_alleles = _extract_read_alleles(alignments, informative_pos)
            hap1_names, hap2_names, amb_names = _assign_reads_to_haplotypes(
                read_alleles, informative_pos
            )
            logger.info(
                f"  Haplotype separation → "
                f"hap1: {len(hap1_names)} reads | "
                f"hap2: {len(hap2_names)} reads | "
                f"ambiguous: {len(amb_names)} reads"
            )

            for hap_label, hap_names in [("haplotype_1", hap1_names), ("haplotype_2", hap2_names)]:
                if len(hap_names) < self.hap_min_reads:
                    logger.warning(
                        f"  {hap_label} has < {self.hap_min_reads} reads "
                        f"({len(hap_names)}) — skipping per-haplotype consensus."
                    )
                    continue
                hap_set  = set(hap_names)
                hap_alns = [a for a in alignments if a["name"] in hap_set]
                h_bp, h_ip = _build_pileup(hap_alns, gene_length)
                h_cons, h_sample, h_pp = _call_consensus(
                    h_bp, h_ip, gene_length,
                    self.min_depth, self.min_freq,
                    self.min_ins_depth, self.min_ins_freq,
                )
                h_val = _validate_reconstruction(
                    h_cons, h_pp, used_ref, h_bp,
                    gene_length, self.min_depth,
                    mv_class, informative_pos,
                    haplotype_label=hap_label,
                )
                haplotypes.append((hap_label, h_cons, h_sample, h_pp, h_val))
                logger.info(f"  {hap_label} consensus grade: {h_val['grade']}")

        # ── variants ──────────────────────────────────────────────────
        variants = self._compute_variants(consensus_ref, used_ref, per_pos)
        snps  = sum(1 for v in variants if v["type"] == "SNP")
        insns = sum(1 for v in variants if v["type"] == "insertion")
        low_c = sum(1 for v in variants if v["type"] == "low_coverage")
        logger.info(
            f"  Variants vs {ref_source}: {snps} SNPs, "
            f"{insns} insertions, {low_c} low-cov positions"
        )

        # ── divergence analysis ───────────────────────────────────────────
        divergence = _analyze_divergence(consensus_ref, used_ref, per_pos, variants)
        if divergence.get("status") == "analyzed":
            logger.info(
                f"  Sample divergence: {divergence['divergence_pct']:.1f}% from reference "
                f"({divergence['divergence_class']}) — {divergence['snps']} SNPs"
            )
            if divergence["divergence_pct"] > 15:
                logger.info("  Note: divergence >15% is expected biology, not an error.")

        # ── boundary extension ────────────────────────────────────────────────
        # Reconstruct ±boundary_extension% beyond each reference-gene end by
        # capturing soft-clipped read bases.  The window is codon-rounded by
        # default to stay in-frame.
        logger.info(
            f"  Codon settings — aware: {self.codon_aware}  "
            f"start: {','.join(sorted(self.start_codons))}  "
            f"stop: {','.join(sorted(self.stop_codons))}"
        )
        if self.boundary_extension > 0:
            raw_ext = max(3, round(gene_length * self.boundary_extension))
            ext_5   = _round_up_to_codon(raw_ext) if self.codon_aware else raw_ext
            ext_3   = ext_5
            codon_note = (
                f" (codon-rounded from {raw_ext} nt)" if self.codon_aware and ext_5 != raw_ext
                else " (codon-rounding disabled)" if not self.codon_aware
                else ""
            )
            logger.info(
                f"  Boundary extension: ±{ext_5} nt each end "
                f"({self.boundary_extension * 100:.0f}% of {gene_length} nt){codon_note}"
            )
            ext_bp, ext_ip = _build_extended_pileup(alignments, gene_length, ext_5, ext_3)
            ext_cons_ref, _ext_sample, _ext_pp = _call_consensus(
                ext_bp, ext_ip, ext_5 + gene_length + ext_3,
                self.min_depth, self.min_freq,
                self.min_ins_depth, self.min_ins_freq,
            )
            ext_5_seq     = ext_cons_ref[:ext_5].lstrip("N")
            ext_3_seq     = ext_cons_ref[ext_5 + gene_length:].rstrip("N")
            captured_5_nt = len(ext_5_seq)
            captured_3_nt = len(ext_3_seq)
            if captured_5_nt or captured_3_nt:
                logger.info(
                    f"  Extension captured: 5ʹ +{captured_5_nt} nt | 3ʹ +{captured_3_nt} nt"
                )
            else:
                logger.info("  No sequence captured beyond reference boundaries.")
        else:
            ext_5_seq, ext_3_seq  = "", ""
            captured_5_nt, captured_3_nt = 0, 0

        # Final output sequence: extension wings + reference-region consensus.
        # The reference-region consensus (from the canonical pileup) is preserved
        # for analytical consistency; only the flanking extensions are new.
        final_seq = ext_5_seq + consensus_ref + ext_3_seq

        # ── length vs reference assessment ────────────────────────────────────
        length_assessment = _assess_length_vs_reference(
            final_seq, gene_length, captured_5_nt, captured_3_nt,
            start_codons=self.start_codons,
            stop_codons=self.stop_codons,
        )
        la = length_assessment
        logger.info(
            f"  Length: {la['status']} — "
            f"{la['recon_length_nt']} nt reconstructed vs {la['ref_length_nt']} nt reference "
            f"({la['diff_pct']:+.1f}%)"
        )
        if la["status"] == "SHORTER_IN_SAMPLE":
            if la["codon_valid"]:
                logger.info(
                    f"  ✓ Valid shorter ORF "
                    f"(start: {la['start_codon']}, stop: {la['stop_codon']}); "
                    "grade not penalised for length"
                )
            else:
                issues = [
                    msg for msg, fail in [
                        ("not divisible by 3",        not la["is_codon_multiple"]),
                        ("no canonical start codon",  not la["has_start_codon"]),
                        ("no canonical stop codon",   not la["has_stop_codon"]),
                    ] if fail
                ]
                logger.info(f"  ⚠ Shorter gene — codon issues: {', '.join(issues)}")
        elif la["status"] == "LONGER_IN_SAMPLE":
            cv_tag = (
                f"valid ORF (start: {la['start_codon']}, stop: {la['stop_codon']})"
                if la["codon_valid"] else "ORF structure not fully valid"
            )
            logger.info(f"  Gene extends beyond reference — {cv_tag}")

        # ── validation & grade ────────────────────────────────────────────────
        # A single validation call that incorporates the length assessment so
        # shorter-but-codon-valid genes are not penalised for uncovered reference tails.
        val = _validate_reconstruction(
            consensus_ref, per_pos, used_ref,
            base_pileup, gene_length, self.min_depth,
            mv_class, informative_pos,
            length_assessment=length_assessment,
        )
        logger.info(
            f"  Grade: {val['grade']} — "
            f"read-support {val['read_support_pct']}%, "
            f"coverage {val['coverage_pct']}%, "
            f"depth CV {val['cv_depth_pct']}%"
        )
        tool_dir = self.recon_dir / gene_safe / tool_label.lower()
        tool_dir.mkdir(parents=True, exist_ok=True)

        # Consensus FASTA — extension bases prepended/appended where captured.
        # A second entry is added when insertions shift the sample sequence.
        final_sample_with_ins = ext_5_seq + sample_seq + ext_3_seq
        fasta_entries = [("sample_consensus", final_seq, gene_stats_row, length_assessment)]
        if sample_seq != consensus_ref:
            fasta_entries.append(
                ("sample_with_insertions", final_sample_with_ins,
                 gene_stats_row, length_assessment)
            )
        self._write_consensus_fasta(tool_dir / f"{gene_safe}_consensus.fasta", fasta_entries)

        if haplotypes:
            hap_entries = [
                (lbl, h_cons, gene_stats_row)
                for lbl, h_cons, h_sample, h_pp, h_val in haplotypes
            ]
            self._write_consensus_fasta(tool_dir / f"{gene_safe}_haplotypes.fasta", hap_entries)

        self._write_depth_table(tool_dir / f"{gene_safe}_depth.tsv", per_pos)

        if variants:
            self._write_variants_tsv(tool_dir / f"{gene_safe}_variants.tsv", variants)

        if divergence and divergence.get("status") == "analyzed":
            _write_divergence_report(
                tool_dir / f"{gene_safe}_sample_vs_reference.txt",
                self.gene_name, consensus_ref, used_ref, divergence, per_pos,
            )

        if self.emit_reads_fasta and read_seqs:
            self._write_reads_fasta(tool_dir / f"{gene_safe}_reads.fasta", read_seqs)

        all_vals = [val] + [h_val for *_, h_val in haplotypes]
        _write_validation_report(tool_dir / f"{gene_safe}_validation.txt", all_vals, self.gene_name)

        # Deeper biological-plausibility check (ORF integrity, GC content, etc.)
        if validate_reconstruction is not None:
            try:
                artifact_report = validate_reconstruction(
                    consensus=consensus_ref,
                    reference=used_ref,
                    per_pos=per_pos,
                    alignments=alignments,
                    variants=variants,
                    gene_name=self.gene_name,
                    gene_type="coding",
                )
                if artifact_report:
                    (tool_dir / f"{gene_safe}_artifact_validation.txt").write_text(
                        artifact_report.summary_report(), encoding="utf-8"
                    )
                    status_val = artifact_report.overall_status.value
                    logger.info(
                        f"  Artifact validation: {status_val} "
                        f"(score: {artifact_report.overall_score:.0f}/100)"
                    )
                    if status_val == "FAIL":
                        logger.warning("  ⚠  Artifact validation FAILED — sequence may not be biologically real.")
                    elif status_val == "WARNING":
                        logger.warning("  ⚠  Artifact validation warnings — review artifact report.")
            except Exception as exc:
                logger.warning(f"  Could not run artifact validation: {exc}")

        if self.emit_plots:
            self._plot_reconstruction(
                tool_label=tool_label, tool_dir=tool_dir, gene_safe=gene_safe,
                per_pos=per_pos, variants=variants, informative_pos=informative_pos,
                haplotypes=haplotypes, mv_class=mv_class, val=val,
                reference=used_ref, divergence=divergence,
            )

        self._write_assembly_graph(
            tool_label=tool_label, tool_dir=tool_dir, gene_safe=gene_safe,
            consensus=consensus_ref, per_pos=per_pos, variants=variants,
            informative_pos=informative_pos, haplotypes=haplotypes,
            mv_class=mv_class, base_pileup=base_pileup, alignments=alignments,
        )

        logger.info(f"  Outputs written to: {tool_dir}/")

        return {
            "alignments":         len(alignments),
            "reads":              len(read_seqs),
            "gene_length":        gene_length,
            "covered_positions":  covered,
            "coverage_pct":       round(covered / gene_length * 100, 2),
            "snps":               snps,
            "insertions":         insns,
            "low_cov_positions":  low_c,
            "ref_source":         ref_source,
            "mv_classification":  mv_class,
            "n_informative_sites": len(informative_pos),
            "n_haplotypes":       len(haplotypes),
            "grade":              val["grade"],
            "read_support_pct":   val["read_support_pct"],
            "length_status":      length_assessment["status"],
            "codon_valid":        length_assessment["codon_valid"],
            "recon_length_nt":    length_assessment["recon_length_nt"],
            "output_dir":         str(tool_dir),
        }

    # ── BLAST / DIAMOND mode ──────────────────────────────────────────────────

    def _run_blast_mode(self) -> dict:
        """Reconstruct using BLAST / DIAMOND tabular output instead of SAM/BAM."""
        logger.info(f"\n{'=' * 70}")
        logger.info(f"Reconstructing (BLAST/DIAMOND mode): {self.gene_name}")
        logger.info(f"{'=' * 70}")
        logger.info(f"  Tabular input: {self.blast_tsv}")

        # Load query sequences (needed when qseq/sseq are not in the tabular file)
        query_seqs: Dict[str, str] = {}
        if self.query_fasta:
            query_seqs = _load_reference_fasta(self.query_fasta)
            logger.info(f"  Query FASTA loaded: {len(query_seqs)} sequences")
        else:
            logger.info(
                "  No -query_fasta provided — will use qseq/sseq columns from "
                "the tabular file.  If those are absent, hits will be skipped."
            )

        # Load reference sequence for validation / variant calling
        ref_seqs = _load_reference_fasta(self.reference_fasta) if self.reference_fasta else {}
        ref_seq  = ref_seqs.get(self.gene_name)
        if self.reference_fasta:
            if ref_seq:
                logger.info(f"  Reference loaded: {len(ref_seq)} bp")
            else:
                for rn, rs in ref_seqs.items():
                    if self.gene_name in rn or rn in self.gene_name:
                        ref_seq = rs
                        logger.info(f"  Reference matched via '{rn}': {len(ref_seq)} bp")
                        break
                if not ref_seq:
                    logger.warning("  No reference match in FASTA — reference identity will be skipped.")

        # Parse the tabular file
        gene_length, read_seqs, alignments = _parse_blast_tabular(
            self.blast_tsv,
            self.gene_name,
            query_seqs,
            field_names=self.blast_format,
        )

        # If reference_fasta was provided, prefer its gene length over max(send)
        if ref_seq and len(ref_seq) > 0:
            gene_length = len(ref_seq)

        if not alignments:
            logger.warning(
                f"  No hits found for '{self.gene_name}' in {self.blast_tsv}."
            )
            return {"gene": self.gene_name, "status": "error",
                    "error": "no hits in BLAST/DIAMOND tabular file"}

        logger.info(
            f"  {len(alignments)} hits loaded, {len(read_seqs)} unique queries, "
            f"gene length: {gene_length} bp"
        )

        ref_source = "provided FASTA" if ref_seq else "estimated from alignments"

        result = self._run_analysis_pipeline(
            tool_label="blast",
            alignments=alignments,
            read_seqs=read_seqs,
            gene_length=gene_length,
            gene_stats_row=None,
            used_ref=ref_seq,
            ref_source=ref_source,
        )

        if result is None:
            return {"gene": self.gene_name, "status": "error", "error": "analysis failed"}
        return {"gene": self.gene_name, "status": "ok", "tools": {"blast": result}}

    # ── main run ──────────────────────────────────────────────────────────────

    def run(self) -> dict:
        # Dispatch to BLAST mode when a tabular file has been provided
        if self.blast_tsv:
            return self._run_blast_mode()

        # ── SAM / BAM mode (default) ──────────────────────────────────
        logger.info(f"\n{'=' * 70}")
        logger.info(f"Reconstructing: {self.gene_name}")
        logger.info(f"{'=' * 70}")

        # Load reference sequence
        ref_seqs = _load_reference_fasta(self.reference_fasta) if self.reference_fasta else {}
        ref_seq  = ref_seqs.get(self.gene_name)
        if self.reference_fasta:
            if ref_seq:
                logger.info(f"  Reference loaded: {len(ref_seq)} bp")
            else:
                for rn, rs in ref_seqs.items():
                    if self.gene_name in rn or rn in self.gene_name:
                        ref_seq = rs
                        logger.info(f"  Reference matched via '{rn}': {len(ref_seq)} bp")
                        break
                if not ref_seq:
                    logger.warning("  No reference match — using MD-tag reconstruction.")

        align_files = self._find_alignment_files()
        if not align_files:
            logger.error("  No alignment files found in raw_outputs/")
            return {"gene": self.gene_name, "status": "error", "error": "no alignment files"}

        tool_results = {}

        for tool_label, bam_path, sam_path in align_files:
            logger.info(f"  ── Tool: {tool_label} ──")

            stats_path     = self._find_stats_files(tool_label)
            gene_stats_row = self._load_gene_stats(stats_path)
            if gene_stats_row:
                logger.info(
                    f"  GeneFior stats — reads: {gene_stats_row.get('Num_Sequences_Mapped')} | "
                    f"gene_cov: {gene_stats_row.get('Gene_Coverage')}% | "
                    f"avg_id: {gene_stats_row.get('Avg_Identity')}% | "
                    f"detected: {gene_stats_row.get('Detected')}"
                )

            sam_iter, proc = _open_sam_source(bam_path, sam_path)
            if sam_iter is None:
                logger.warning(f"  No SAM/BAM found — skipping {tool_label}.")
                continue
            try:
                gene_length, read_seqs, alignments = self._parse_sam_for_gene(sam_iter)
            finally:
                try:
                    sam_iter.close()
                except Exception:
                    pass
                if proc:
                    try:
                        proc.wait()
                    except Exception:
                        pass

            if gene_length == 0 or not alignments:
                logger.warning(
                    f"  Gene not found in SAM/no alignments for {tool_label} — skipping."
                )
                continue

            logger.info(
                f"  {len(alignments)} alignments, {len(read_seqs)} reads, "
                f"gene length: {gene_length} bp"
            )

            ref_source = "provided FASTA" if ref_seq else "MD-tag reconstruction"
            result = self._run_analysis_pipeline(
                tool_label=tool_label,
                alignments=alignments,
                read_seqs=read_seqs,
                gene_length=gene_length,
                gene_stats_row=gene_stats_row,
                used_ref=ref_seq,
                ref_source=ref_source,
            )
            if result:
                tool_results[tool_label] = result

        return {"gene": self.gene_name, "status": "ok", "tools": tool_results}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Report
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _write_summary_report(report_path: Path, results: list, args: argparse.Namespace):
    with open(report_path, "w") as fh:
        fh.write("=" * 70 + "\n")
        fh.write("GENEFIOR GENE RECONSTRUCTION REPORT\n")
        fh.write(f"Generated:       {datetime.now().isoformat()}\n")
        if getattr(args, "blast_tsv", None):
            fh.write(f"Input mode:      BLAST/DIAMOND tabular\n")
            fh.write(f"Tabular file:    {args.blast_tsv}\n")
            if getattr(args, "query_fasta", None):
                fh.write(f"Query FASTA:     {args.query_fasta}\n")
        else:
            fh.write(f"GeneFior output: {args.output_dir}\n")
        fh.write(f"Min depth:       {args.min_depth}\n")
        fh.write(f"Min freq:        {args.min_freq}\n")
        fh.write(
            f"Reference FASTA: "
            f"{args.reference_fasta or 'Not provided (MD-tag reconstruction)'}\n"
        )
        fh.write(
            f"Boundary ext:    {args.boundary_extension * 100:.0f}%  "
            f"codon-aware: {not args.no_codon_aware}  "
            f"start: {args.start_codons}  stop: {args.stop_codons}\n"
        )
        fh.write("=" * 70 + "\n\n")

        for res in results:
            fh.write(f"Gene: {res['gene']}\n")
            fh.write(f"  Status: {res.get('status', '?')}\n")
            if res.get("status") == "error":
                fh.write(f"  Error: {res.get('error')}\n\n")
                continue
            for tool, info in res.get("tools", {}).items():
                fh.write(f"  Tool: {tool}\n")
                fh.write(f"    Reads mapped:          {info['reads']}\n")
                fh.write(f"    Gene length:           {info['gene_length']} bp\n")
                fh.write(
                    f"    Covered positions:     {info['covered_positions']} / "
                    f"{info['gene_length']} ({info['coverage_pct']}%)\n"
                )
                fh.write(f"    SNPs (vs reference):   {info['snps']}\n")
                fh.write(f"    Insertions:            {info['insertions']}\n")
                fh.write(f"    Low-cov positions:     {info['low_cov_positions']}\n")
                fh.write(f"    Reference source:      {info['ref_source']}\n")
                fh.write(
                    f"    Multi-version:         {info['mv_classification']} "
                    f"({info['n_informative_sites']} informative sites)\n"
                )
                fh.write(f"    Haplotypes produced:   {info['n_haplotypes']}\n")
                fh.write(
                    f"    Grade:                 {info['grade']}  "
                    f"(read-support {info['read_support_pct']}%)\n"
                )
                fh.write(
                    f"    Length vs reference:   {info.get('length_status', 'N/A')}  "
                    f"({info.get('recon_length_nt', '?')} nt reconstructed)  "
                    f"codon-valid: {info.get('codon_valid', 'N/A')}\n"
                )
                fh.write(f"    Output:                {info['output_dir']}\n")
            fh.write("\n")

        fh.write(
            "GRADE KEY\n"
            "  A   High-confidence, single clean version\n"
            "  B   Good, minor issues\n"
            "  C   Moderate confidence\n"
            "  C*  Multiple versions detected — blended consensus; see haplotypes file\n"
            "  F   Insufficient coverage\n"
            "\n"
            "LENGTH STATUS\n"
            "  MATCHES_REFERENCE  — Reconstructed length equals reference\n"
            "  SHORTER_IN_SAMPLE  — Sample gene shorter; grade adjusted when codon-valid\n"
            "  LONGER_IN_SAMPLE   — Sample gene extends beyond reference (captured)\n"
        )


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CLI
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    parser = argparse.ArgumentParser(
        prog="GeneFior-Reconstruct",
        description=(
            f"GeneFior-Reconstruct {RECONSTRUCT_VERSION} — "
            "Reconstruct in-sample gene sequences from a GeneFior run or from "
            "BLAST / DIAMOND tabular output, including detection of multiple gene versions."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    src = parser.add_argument_group(
        "Input source (use ONE of the two blocks below)"
    )
    src.add_argument(
        "-output_dir",
        help="Path to a completed GeneFior output directory (SAM/BAM mode).",
    )
    src.add_argument(
        "-blast_tsv",
        help=(
            "BLAST or DIAMOND tabular output file (format 6) to use instead of "
            "SAM/BAM files.  When this flag is set -output_dir is not required."
        ),
    )

    req = parser.add_argument_group("Required")
    req.add_argument(
        "-gene", required=True, action="append", dest="genes",
        help=(
            "Gene name to reconstruct (exactly as in the stats/detection files "
            "or as the subject ID in the BLAST/DIAMOND output). "
            "Repeat the flag or use a comma-separated list for multiple genes."
        ),
    )

    opt = parser.add_argument_group("Optional")
    opt.add_argument(
        "-reference_fasta",
        help=(
            "Gene reference FASTA used in the original GeneFior / BLAST run.  "
            "Enables accurate reference length, reference identity scoring, and "
            "variant calling vs the true reference sequence."
        ),
    )
    opt.add_argument(
        "-query_fasta",
        help=(
            "BLAST / DIAMOND mode only.  Original query FASTA file (reads or "
            "contigs) used in the BLAST/DIAMOND search.  Required when the "
            "tabular file does not include qseq/sseq columns (standard format 6 "
            "does not include them by default)."
        ),
    )
    opt.add_argument(
        "-blast_format",
        help=(
            "BLAST / DIAMOND mode only.  Comma-separated list of column names "
            "matching the fields in -blast_tsv when a non-default format was "
            "used.  Default (standard format 6): "
            "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore.  "
            "Add 'qseq' and/or 'sseq' to enable gap-aware CIGAR derivation; "
            "add 'slen' to obtain exact gene length."
        ),
    )
    opt.add_argument("-recon_dir",
                     help="Output directory (default: <output_dir>/reconstruction/).")
    opt.add_argument("-tool",
                     help="SAM/BAM mode only.  Alignment tool: bowtie2 | bwa | minimap2 (default: all).")
    opt.add_argument("-min_depth", type=int, default=3,
                     help="Min read depth to call a consensus base (default: 3).")
    opt.add_argument("-min_freq", type=float, default=0.5,
                     help="Min base frequency for unambiguous call (default: 0.5).")
    opt.add_argument("-min_ins_depth", type=int, default=3,
                     help="Min insertion support depth (default: 3).")
    opt.add_argument("-min_ins_freq", type=float, default=0.5,
                     help="Min insertion frequency (default: 0.5).")
    opt.add_argument("-min_bimodal_freq", type=float, default=0.20,
                     help="Min minor-allele frequency to flag a site as informative "
                          "for multi-version detection (default: 0.20).")
    opt.add_argument("-min_informative_sites", type=int, default=3,
                     help="Min number of informative sites required to classify "
                          "a gene as 'multi-version' (default: 3).")
    opt.add_argument("-hap_min_reads", type=int, default=5,
                     help="Min reads required per haplotype group to produce a "
                          "per-haplotype consensus (default: 5).")
    opt.add_argument("-no_reads_fasta", action="store_true",
                     help="Do not write the reads FASTA file.")
    opt.add_argument("-no_plots", action="store_true",
                     help="Do not generate reconstruction PNG plots "
                          "(requires matplotlib + numpy).")
    opt.add_argument("-no_hap_uncertain", action="store_true",
                     help="Do not attempt haplotype separation when the "
                          "multi-version classification is 'uncertain' "
                          "(only separate when classified as 'multi').")

    # ── Boundary extension & codon-structure options ─────────────────────────
    bnd = parser.add_argument_group(
        "Boundary extension & codon structure",
        description=(
            "Controls how far past the reference gene boundaries the "
            "reconstruction attempts to reach, and which codons are treated "
            "as valid gene starts / stops when assessing the reconstructed "
            "sequence."
        ),
    )
    bnd.add_argument(
        "-boundary_extension", type=float, default=0.20,
        metavar="FRAC",
        help=(
            "Fraction of the reference gene length to extend the reconstruction "
            "window at both the 5ʹ and 3ʹ ends (default: 0.20 = 20%%).  "
            "Soft-clipped read bases are used to populate the extension region. "
            "Set to 0 to disable.  Example: -boundary_extension 0.30"
        ),
    )
    bnd.add_argument(
        "-no_codon_aware", action="store_true",
        help=(
            "Disable codon-aware extension rounding (default: enabled). "
            "By default the extension window is rounded UP to the nearest "
            "multiple of 3 nt so it aligns to a codon boundary.  Use this "
            "flag to use the raw (un-rounded) window size instead."
        ),
    )
    bnd.add_argument(
        "-start_codons",
        default=_DEFAULT_START_CODONS_STR,
        metavar="CODONS",
        help=(
            f"Comma-separated list of triplets accepted as gene start codons "
            f"(default: {_DEFAULT_START_CODONS_STR}).  "
            "Used when assessing whether a shorter or longer reconstructed "
            "sequence is a valid ORF.  "
            "Example: -start_codons ATG,GTG"
        ),
    )
    bnd.add_argument(
        "-stop_codons",
        default=_DEFAULT_STOP_CODONS_STR,
        metavar="CODONS",
        help=(
            f"Comma-separated list of triplets accepted as gene stop codons "
            f"(default: {_DEFAULT_STOP_CODONS_STR}).  "
            "Used when assessing whether a shorter or longer reconstructed "
            "sequence is a valid ORF.  "
            "Example: -stop_codons TAA,TAG,TGA"
        ),
    )

    opt.add_argument("-v", "--version", action="version",
                     version=f"GeneFior-Reconstruct {RECONSTRUCT_VERSION}")

    args = parser.parse_args()

    # ── validation ────────────────────────────────────────────────────────────
    if not args.blast_tsv and not args.output_dir:
        parser.error("Provide either -output_dir (SAM/BAM mode) or -blast_tsv (BLAST/DIAMOND mode).")
    if args.blast_tsv and not os.path.isfile(args.blast_tsv):
        parser.error(f"-blast_tsv file not found: {args.blast_tsv}")
    if args.blast_tsv and args.query_fasta and not os.path.isfile(args.query_fasta):
        parser.error(f"-query_fasta file not found: {args.query_fasta}")
    # Validate codon lists early so the user gets a clear error message
    try:
        _parse_codon_list(args.start_codons, "start")
    except ValueError as exc:
        parser.error(f"-start_codons: {exc}")
    try:
        _parse_codon_list(args.stop_codons, "stop")
    except ValueError as exc:
        parser.error(f"-stop_codons: {exc}")

    # ── flatten gene list ─────────────────────────────────────────────────────
    gene_names = list(dict.fromkeys(
        g.strip() for raw in args.genes for g in raw.split(",") if g.strip()
    ))

    base_dir   = args.output_dir or os.path.dirname(os.path.abspath(args.blast_tsv))
    recon_dir  = Path(args.recon_dir) if args.recon_dir else Path(base_dir) / "reconstruction"
    recon_dir.mkdir(parents=True, exist_ok=True)

    # ── file logger ───────────────────────────────────────────────────────────
    fh = logging.FileHandler(recon_dir / "GeneFior_Reconstruct.log")
    fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(fh)

    logger.info(f"GeneFior-Reconstruct {RECONSTRUCT_VERSION}")
    logger.info(f"Genes requested: {gene_names}")
    logger.info(f"Output dir:      {recon_dir}")
    if args.blast_tsv:
        logger.info("Mode: BLAST/DIAMOND tabular")
    else:
        if not _samtools_available():
            logger.info("samtools not on PATH — parsing .sam files directly.")
    logger.info(
        f"Boundary extension: {args.boundary_extension * 100:.0f}%  "
        f"codon-aware: {not args.no_codon_aware}  "
        f"start codons: {args.start_codons}  "
        f"stop codons: {args.stop_codons}"
    )

    # ── parse blast_format ────────────────────────────────────────────────────
    blast_format = None
    if args.blast_format:
        blast_format = [f.strip() for f in args.blast_format.split(",") if f.strip()]

    # ── run ───────────────────────────────────────────────────────────────────
    all_results = []
    for gene_name in gene_names:
        rec = GeneReconstructor(
            genefior_output_dir=args.output_dir or ".",
            gene_name=gene_name,
            recon_dir=str(recon_dir),
            reference_fasta=args.reference_fasta,
            tool=args.tool,
            min_depth=args.min_depth,
            min_freq=args.min_freq,
            min_insertion_depth=args.min_ins_depth,
            min_insertion_freq=args.min_ins_freq,
            min_bimodal_freq=args.min_bimodal_freq,
            min_informative_sites=args.min_informative_sites,
            hap_min_reads=args.hap_min_reads,
            emit_reads_fasta=not args.no_reads_fasta,
            emit_plots=not args.no_plots,
            hap_on_uncertain=not args.no_hap_uncertain,
            boundary_extension=args.boundary_extension,
            codon_aware=not args.no_codon_aware,
            start_codons=args.start_codons,
            stop_codons=args.stop_codons,
            blast_tsv=args.blast_tsv,
            query_fasta=getattr(args, "query_fasta", None),
            blast_format=blast_format,
        )
        all_results.append(rec.run())

    _write_summary_report(recon_dir / "reconstruction_report.txt", all_results, args)
    logger.info(f"\nSummary report: {recon_dir / 'reconstruction_report.txt'}")
    logger.info("Thank you for using GeneFior-Reconstruct!")
    logger.info("Documentation: https://github.com/NickJD/GeneFior")


if __name__ == "__main__":
    main()

