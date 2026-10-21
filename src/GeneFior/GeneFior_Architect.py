"""Compare GeneFior reconstructions with an assembly or database reference.

Architect is an offline adjudication tool: it aligns each reconstructed
sequence either to matched sample assemblies or directly to its reference
database sequence, then reports identity, aligned fraction, likely failure,
and practical next steps.
It discovers samples and reconstruction files at run time, so an AMRFior run
that is still adding sample directories can be evaluated again later.

Requires the optional ``architect`` dependency (Biopython).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional

try:
    from Bio.Align import PairwiseAligner
except ImportError:  # fail with a helpful CLI message rather than at import time
    PairwiseAligner = None


ARCHITECT_VERSION = "0.3.0"
FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas"}
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass
class AlignmentResult:
    score: float
    identity_pct: float
    query_coverage_pct: float
    target_coverage_pct: float
    query_length: int
    aligned_query_bases: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    strand: str
    gaps: int


def open_text(path: Path):
    return gzip.open(path, "rt") if path.name.endswith(".gz") else path.open()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: Optional[str] = None
    pieces: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(pieces).upper()))
                name, pieces = line[1:].split()[0], []
            elif name is not None:
                pieces.append(line)
    if name is not None:
        records.append((name, "".join(pieces).upper()))
    return records


def make_aligner():
    aligner = PairwiseAligner(mode="local")
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    return aligner


def align_to_contig(query: str, target: str, aligner, sensitive_fallback: bool = True) -> Optional[AlignmentResult]:
    """Find the best local alignment using both query orientations."""
    if not query or not target:
        return None
    best = None
    for strand, oriented in (("+", query), ("-", query.translate(RC_TABLE)[::-1])):
        # Seed exact k-mers first so dynamic programming runs on a small locus
        # window rather than a multi-megabase chromosome.  At least two
        # independently placed seeds must agree on an approximate start.
        k = min(19, max(11, len(oriented) // 12))
        step = max(k, 30)
        seed_positions = list(range(0, max(1, len(oriented) - k + 1), step))
        last = max(0, len(oriented) - k)
        if not seed_positions or seed_positions[-1] != last:
            seed_positions.append(last)
        starts = Counter()
        for qpos in seed_positions:
            seed = oriented[qpos:qpos + k]
            if len(seed) != k or any(base not in "ACGT" for base in seed):
                continue
            tpos = target.find(seed)
            found = 0
            while tpos >= 0 and found < 6:
                starts[tpos - qpos] += 1
                tpos = target.find(seed, tpos + 1)
                found += 1
        candidate_ranges = []
        for offset, votes in starts.most_common(8):
            if votes < 2:
                continue
            start = max(0, offset - 120)
            end = min(len(target), offset + len(oriented) + 120)
            if end > start:
                candidate_ranges.append((start, end))
        # If sparse long seeds miss a divergent allele, retry with dense short
        # seeds. This remains a bounded-window search over chromosome-sized
        # assemblies, while tolerating substantial nucleotide divergence.
        # Always add sensitive seeds to cross-check sparse exact candidates.
        fallback_k = min(9, max(4, len(oriented) // 45))
        fallback_starts = Counter()
        fallback_step = 1 if len(oriented) <= 700 else max(8, fallback_k)
        positions = ([] if candidate_ranges and not sensitive_fallback else
                    list(range(0, max(1, len(oriented) - fallback_k + 1), fallback_step)))
        last = max(0, len(oriented) - fallback_k)
        if positions and positions[-1] != last:
            positions.append(last)
        for qpos in positions:
            seed = oriented[qpos:qpos + fallback_k]
            if len(seed) != fallback_k or any(base not in "ACGT" for base in seed):
                continue
            if oriented.count(seed) > 1:
                continue
            # Short seeds can occur at several homologous loci. Vote for each
            # bounded occurrence, but skip highly repetitive seeds.
            matches = []
            tpos = target.find(seed)
            while tpos >= 0 and len(matches) <= 12:
                matches.append(tpos)
                tpos = target.find(seed, tpos + 1)
            if len(matches) > 12:
                continue
            for tpos in matches:
                fallback_starts[tpos - qpos] += 1
        clustered = []
        for offset in fallback_starts:
            support = sum(votes for nearby, votes in fallback_starts.items()
                          if abs(nearby - offset) <= 30)
            clustered.append((support, offset))
        seen_offsets = set()
        for support, offset in sorted(clustered, reverse=True):
            if support < 2 or any(abs(offset - prior) <= 30 for prior in seen_offsets):
                continue
            seen_offsets.add(offset)
            start = max(0, offset - 120)
            end = min(len(target), offset + len(oriented) + 120)
            if end > start:
                candidate_ranges.append((start, end))
            if len(seen_offsets) >= (64 if sensitive_fallback and len(oriented) <= 700 else 12):
                break
        if not candidate_ranges and len(target) <= 1_000:
            candidate_ranges = [(0, len(target))]

        for window_start, window_end in candidate_ranges:
            window = target[window_start:window_end]
            try:
                alignment = aligner.align(window, oriented)[0]
            except (IndexError, ValueError):
                continue
            coords = alignment.coordinates
            q_aligned = t_aligned = matches = gaps = 0
            for i in range(coords.shape[1] - 1):
                dt = int(coords[0, i + 1] - coords[0, i])
                dq = int(coords[1, i + 1] - coords[1, i])
                if dt and dq:
                    span = min(dt, dq)
                    target_part = window[int(coords[0, i]):int(coords[0, i]) + span]
                    query_part = oriented[int(coords[1, i]):int(coords[1, i]) + span]
                    q_aligned += span
                    t_aligned += span
                    matches += sum(a == b and a in "ACGT" for a, b in zip(target_part, query_part))
                    gaps += abs(dt - dq)
                else:
                    gaps += dt + dq
            coverage = (q_aligned / len(query) * 100) if query else 0.0
            target_coverage = (t_aligned / len(target) * 100) if target else 0.0
            identity = (matches / q_aligned * 100) if q_aligned else 0.0
            oriented_q_start = int(coords[1, 0])
            oriented_q_end = int(coords[1, -1])
            q_start, q_end = (
                (oriented_q_start, oriented_q_end) if strand == "+"
                else (len(query) - oriented_q_end, len(query) - oriented_q_start)
            )
            result = AlignmentResult(
                score=float(alignment.score), identity_pct=identity,
                query_coverage_pct=coverage, target_coverage_pct=target_coverage,
                query_length=len(query),
                aligned_query_bases=q_aligned,
                query_start=q_start, query_end=q_end,
                target_start=window_start + int(coords[0, 0]),
                target_end=window_start + int(coords[0, -1]),
                strand=strand, gaps=gaps,
            )
            if best is None or (result.score, result.query_coverage_pct, result.identity_pct) > (
                best.score, best.query_coverage_pct, best.identity_pct
            ):
                best = result
    return best


def locate_reference(reference: str, records: list[tuple[str, str]], aligner,
                     sensitive_fallback: bool = True):
    best = None
    for name, sequence in records:
        hit = align_to_contig(reference, sequence, aligner, sensitive_fallback=False)
        # Pay for dense short-seed candidates only when the sparse exact-seed
        # search misses or returns a short local fragment on a short query.
        if sensitive_fallback and len(reference) <= 700 and (hit is None or hit.query_coverage_pct < 65):
            retry = align_to_contig(reference, sequence, aligner, sensitive_fallback=True)
            if retry and (hit is None or retry.score > hit.score):
                hit = retry
        if hit and (best is None or hit.score > best[0].score):
            best = (hit, name, sequence)
    return best


def get_reference(gene: str, fasta: Path) -> Optional[str]:
    records = read_fasta(fasta)
    by_name = dict(records)
    if gene in by_name:
        return by_name[gene]
    safe_gene = re.sub(r"[^\w\-.]", "_", gene)
    for name, seq in records:
        if re.sub(r"[^\w\-.]", "_", name) == safe_gene:
            return seq
    return None


def find_assemblies(truth_root: Path, sample_id: str) -> list[tuple[str, Path]]:
    isolate = sample_id.split("_")[0]
    found: list[tuple[str, Path]] = []
    long_dir = truth_root / "long_reads" / isolate
    for filename in ("assembly.fasta", "assembly.fna"):
        candidate = long_dir / filename
        if candidate.is_file():
            found.append(("long_read", candidate))
            break

    enhanced_dir = truth_root / "enhanced_contigs"
    # Prefer the direct sample assembly; annotation subdirectories are a
    # fallback for datasets that package the assembly there instead.
    candidates = []
    if enhanced_dir.is_dir():
        candidates.extend(p for p in enhanced_dir.iterdir()
                          if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES
                          and p.stem == sample_id)
        if not candidates:
            candidates.extend(p for p in enhanced_dir.glob(f"{sample_id}*")
                              if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES)
        if not candidates:
            for subdir in enhanced_dir.glob(f"{sample_id}*"):
                if subdir.is_dir():
                    candidates.extend(p for p in subdir.iterdir()
                                      if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES)
    if candidates:
        candidates.sort(key=lambda p: (p.stem != sample_id, p.suffix.lower() != ".fasta", p.name))
        found.append(("enhanced_contig", candidates[0]))
    return found


def classify(result: Optional[AlignmentResult], ref_hit, n_pct: float,
             reference_available: bool = True) -> tuple[str, str]:
    if not reference_available:
        return "UNASSESSABLE", "The reference sequence for this gene was not found in the supplied FASTA."
    if ref_hit is None:
        return "SEARCH_UNRESOLVED", "No reference-like locus was found after sensitive seeded search; this is not evidence of gene absence."
    if ref_hit.query_coverage_pct < 65 or ref_hit.identity_pct < 75:
        return "UNASSESSABLE", "The reference gene was not located confidently in this assembly."
    if result is None or result.query_coverage_pct < 65 or result.identity_pct < 80:
        return "FAILED", "The assembly supports a reference-like locus, but the reconstruction does not match it sufficiently."
    if n_pct > 10 or result.query_coverage_pct < 90 or result.identity_pct < 95:
        return "REVIEW", "The reconstruction is partial, ambiguous, or substantially different; inspect coverage and allele evidence."
    if result.query_coverage_pct >= 95 and result.identity_pct >= 99:
        return "PASS", "High identity and coverage against the matched assembly."
    return "LIKELY_PASS", "The reconstruction matches the assembly well; review remaining differences and boundaries."


def recommendations(status: str, n_pct: float, result: Optional[AlignmentResult], ref_hit) -> str:
    actions = []
    if status in {"NOT_CONFIRMED", "SEARCH_UNRESOLVED"}:
        actions.append("review whole-assembly reconstruction hits and raw read evidence; a search miss does not establish gene absence")
    elif status == "UNASSESSABLE":
        actions.append("check reference sequence and assembly fragmentation at the expected locus")
    elif status == "FAILED":
        actions.append("inspect per-position depth, read support, and alignment evidence at the truth locus")
        actions.append("retry with an alternate mapper or adjusted minimum depth if support is weak")
    elif status == "REVIEW":
        if n_pct > 10:
            actions.append("investigate low-depth N positions and consider additional reads")
        if result and result.query_coverage_pct < 90:
            actions.append("review the uncovered gene ends and boundary-extension settings")
        if result and result.identity_pct < 95:
            actions.append("inspect variants and allele ambiguity; do not force the reference allele")
    else:
        actions.append("retain this reconstruction and record the assembly source used for validation")
    if ref_hit and ref_hit.query_coverage_pct < 90:
        actions.append("assembly locus is itself partial; treat this comparison as lower confidence")
    return "; ".join(dict.fromkeys(actions))


def evaluate_file(fasta_path: Path, gene: str, sample_id: str,
                  reference_fasta: Path, assemblies: list[tuple[str, Path]], aligner,
                  assembly_cache: dict, locus_cache: dict):
    reference = get_reference(gene, reference_fasta)
    records = read_fasta(fasta_path)
    rows = []
    for record_name, sequence in records:
        # A second with-insertions record is a diagnostic alternative to the
        # primary consensus; keep both but label their role explicitly.
        if not sequence:
            continue
        n_pct = sequence.count("N") / len(sequence) * 100
        for assembly_kind, assembly_path in assemblies:
            assembly_key = str(assembly_path)
            if assembly_key not in assembly_cache:
                assembly_cache[assembly_key] = read_fasta(assembly_path)
            contigs = assembly_cache[assembly_key]
            locus_key = (assembly_key, gene)
            if locus_key not in locus_cache:
                locus_cache[locus_key] = locate_reference(reference, contigs, aligner) if reference else None
            ref_hit = locus_cache[locus_key]
            # Search the reconstruction independently across every contig. A
            # gene can be present away from the expected reference locus.
            global_key = (assembly_key, "global", sequence)
            if global_key not in locus_cache:
                locus_cache[global_key] = locate_reference(
                    sequence, contigs, aligner, sensitive_fallback=False
                )
            global_hit = locus_cache[global_key]
            if global_hit:
                global_alignment, global_contig_name, _ = global_hit
            else:
                global_alignment, global_contig_name = None, ""
            if ref_hit:
                ref_alignment, contig_name, contig_seq = ref_hit
                margin = max(100, len(reference) // 4)
                left = max(0, ref_alignment.target_start - margin)
                right = min(len(contig_seq), ref_alignment.target_end + margin)
                truth_region = contig_seq[left:right]
                result = align_to_contig(sequence, truth_region, aligner)
                ref_locus_start = ref_alignment.target_start + 1
                ref_locus_end = ref_alignment.target_end
                if result:
                    recon_locus_start = left + result.target_start + 1
                    recon_locus_end = left + result.target_end
                    boundary_match = result.strand == ref_alignment.strand
                    if boundary_match:
                        orient = 1 if ref_alignment.strand == "+" else -1
                        ref_5 = ref_locus_start if orient == 1 else ref_locus_end
                        ref_3 = ref_locus_end if orient == 1 else ref_locus_start
                        recon_5 = recon_locus_start if orient == 1 else recon_locus_end
                        recon_3 = recon_locus_end if orient == 1 else recon_locus_start
                        delta_5 = (ref_5 - recon_5) * orient
                        delta_3 = (recon_3 - ref_3) * orient
                        unaligned_5 = result.query_start if result.strand == "+" else len(sequence) - result.query_end
                        unaligned_3 = len(sequence) - result.query_end if result.strand == "+" else result.query_start
                    else:
                        delta_5 = delta_3 = unaligned_5 = unaligned_3 = None
                else:
                    recon_locus_start = recon_locus_end = None
                    delta_5 = delta_3 = unaligned_5 = unaligned_3 = None
                    boundary_match = False
            else:
                contig_name, result, ref_alignment = "", None, None
                ref_locus_start = ref_locus_end = None
                recon_locus_start = recon_locus_end = None
                delta_5 = delta_3 = unaligned_5 = unaligned_3 = None
                boundary_match = False
            status, explanation = classify(result, ref_alignment, n_pct,
                                           reference_available=reference is not None)
            rows.append({
                "comparison_mode": "assembly",
                "sample_id": sample_id, "gene": gene,
                "reconstruction_file": str(fasta_path), "record": record_name,
                "assembly_kind": assembly_kind, "assembly_file": str(assembly_path),
                "comparison_target": f"{assembly_kind}:{contig_name}",
                "truth_contig": contig_name, "status": status,
                "reconstruction_global_contig": global_contig_name,
                "reconstruction_global_start_1based": global_alignment.target_start + 1 if global_alignment else None,
                "reconstruction_global_end_1based": global_alignment.target_end if global_alignment else None,
                "reconstruction_global_identity_pct": round(global_alignment.identity_pct, 3) if global_alignment else None,
                "reconstruction_global_coverage_pct": round(global_alignment.query_coverage_pct, 3) if global_alignment else None,
                "reconstruction_length": len(sequence), "n_pct": round(n_pct, 3),
                "identity_pct": round(result.identity_pct, 3) if result else None,
                "reconstruction_coverage_pct": round(result.query_coverage_pct, 3) if result else None,
                "truth_reference_identity_pct": round(ref_alignment.identity_pct, 3) if ref_alignment else None,
                "truth_reference_coverage_pct": round(ref_alignment.query_coverage_pct, 3) if ref_alignment else None,
                "truth_reference_strand": ref_alignment.strand if ref_alignment else "",
                "truth_reference_locus_start_1based": ref_locus_start,
                "truth_reference_locus_end_1based": ref_locus_end,
                "reconstruction_locus_start_1based": recon_locus_start,
                "reconstruction_locus_end_1based": recon_locus_end,
                "reconstruction_strand": result.strand if result else "",
                "gene_5prime_extension_bp": delta_5,
                "gene_3prime_extension_bp": delta_3,
                "unaligned_5prime_sequence_bp": unaligned_5,
                "unaligned_3prime_sequence_bp": unaligned_3,
                "boundary_strands_match": boundary_match,
                "aligned_bases": result.aligned_query_bases if result else 0,
                "indel_bases": result.gaps if result else None,
                "explanation": explanation,
                "recommended_action": recommendations(status, n_pct, result, ref_alignment),
            })
    return rows


def classify_database(result: Optional[AlignmentResult], n_pct: float,
                      reference_available: bool = True) -> tuple[str, str]:
    """Classify a reconstruction against its named database sequence."""
    if not reference_available:
        return "UNASSESSABLE", "The reference sequence for this gene was not found in the supplied FASTA."
    if result is None or result.target_coverage_pct < 65 or result.identity_pct < 80:
        return "FAILED", "The reconstruction does not sufficiently match the selected reference database sequence."
    if n_pct > 10 or result.target_coverage_pct < 90 or result.identity_pct < 95:
        return "REVIEW", "The reconstruction is partial, ambiguous, or substantially different from the selected reference."
    if result.target_coverage_pct >= 95 and result.identity_pct >= 99:
        max_expected_flank = (result.query_length - result.aligned_query_bases)
        if max_expected_flank > result.query_length * 0.45:
            return "REVIEW", "The reference gene matches, but substantial extra reconstruction sequence cannot be validated by database-only comparison."
        return "PASS", "High reference-gene identity and coverage; database-only comparison cannot validate flanking sequence."
    return "LIKELY_PASS", "The reconstruction matches the selected reference well; review remaining differences and boundaries."


def evaluate_database_file(fasta_path: Path, gene: str, sample_id: str,
                           reference_fasta: Path, aligner):
    """Compare consensus records directly to their named reference sequence."""
    reference = get_reference(gene, reference_fasta)
    rows = []
    for record_name, sequence in read_fasta(fasta_path):
        if not sequence:
            continue
        result = align_to_contig(sequence, reference, aligner) if reference else None
        n_pct = sequence.count("N") / len(sequence) * 100
        status, explanation = classify_database(result, n_pct, reference is not None)
        if result:
            # Coordinates are on the target reference, which is the comparison
            # target in this mode; convert the aligned span to 1-based inclusive.
            locus_start, locus_end = result.target_start + 1, result.target_end
            if result.strand == "+":
                delta_5, delta_3 = -result.target_start, len(reference) - result.target_end
            else:
                delta_5, delta_3 = result.target_start - len(reference), result.target_end
        else:
            locus_start = locus_end = delta_5 = delta_3 = None
        rows.append({
            "comparison_mode": "database",
            "sample_id": sample_id, "gene": gene,
            "reconstruction_file": str(fasta_path), "record": record_name,
            "assembly_kind": "", "assembly_file": str(reference_fasta),
            "comparison_target": gene, "truth_contig": "", "status": status,
            "reconstruction_length": len(sequence), "n_pct": round(n_pct, 3),
            "identity_pct": round(result.identity_pct, 3) if result else None,
            "reconstruction_coverage_pct": round(result.query_coverage_pct, 3) if result else None,
            "database_reference_coverage_pct": round(result.target_coverage_pct, 3) if result else None,
            "truth_reference_identity_pct": 100.0 if reference else None,
            "truth_reference_coverage_pct": 100.0 if reference else None,
            "truth_reference_strand": "+" if reference else "",
            "truth_reference_locus_start_1based": 1 if reference else None,
            "truth_reference_locus_end_1based": len(reference) if reference else None,
            "reconstruction_locus_start_1based": locus_start,
            "reconstruction_locus_end_1based": locus_end,
            "reconstruction_strand": result.strand if result else "",
            "gene_5prime_extension_bp": delta_5,
            "gene_3prime_extension_bp": delta_3,
            "unaligned_5prime_sequence_bp": result.query_start if result else None,
            "unaligned_3prime_sequence_bp": len(sequence) - result.query_end if result else None,
            "boundary_strands_match": bool(result),
            "aligned_bases": result.aligned_query_bases if result else 0,
            "indel_bases": result.gaps if result else None,
            "explanation": explanation,
            "recommended_action": (
                "review reference sequence identity and reconstruction boundaries"
                if status in {"REVIEW", "FAILED", "UNASSESSABLE"}
                else "retain this reconstruction and record the database reference used"
            ),
        })
    return rows


def discover_reconstructions(reconstruction_dir: Path, sample_filter: Optional[str]):
    found = []
    for path in sorted(reconstruction_dir.rglob("*_consensus.fasta")):
        if path.parent.name.lower() not in {"bwa", "bowtie2", "minimap2", "blast"}:
            continue
        gene = path.parent.parent.name
        # Sample-level layout is root/sample/gene/tool; flat layout requires
        # the caller to pass --sample-id.
        sample_id = path.parents[2].name if len(path.parents) > 2 else ""
        if sample_filter and sample_id != sample_filter:
            continue
        found.append((sample_id, gene, path))
    return found


def discover_candidates(amrfior_root: Path):
    """Rank genes with read-level evidence so new completed runs join automatically."""
    rank = {
        "EXACT_ALLELE_DETECTED": 1,
        "EXACT_PROTEIN_DETECTED": 1,
        "CANDIDATE_ALLELE_DETECTED": 2,
        "PROFILE_DETECTED": 3,
        "MIXED_OR_MOSAIC": 4,
        "ALLELE_LIKE": 5,
        "FAMILY_DETECTED": 6,
    }
    candidates = []
    if not amrfior_root.is_dir():
        return candidates
    for sample_dir in sorted(path for path in amrfior_root.iterdir() if path.is_dir()):
        evidence_path = sample_dir / "resfinder_evidence_matrix.tsv"
        if not evidence_path.is_file():
            continue  # AMRFior may still be writing this sample.
        detection_path = sample_dir / "resfinder_detection_matrix.tsv"
        detections = {}
        if detection_path.is_file():
            with detection_path.open(newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    try:
                        detections[row.get("Gene", "")] = int(row.get("Total_Detections", 0) or 0)
                    except ValueError:
                        pass
        with evidence_path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                status = row.get("Best_Evidence_Status", "")
                if status not in rank:
                    continue
                try:
                    evidence = int(row.get("Evidence_Detections", 0) or 0)
                    allele = int(row.get("Candidate_Allele_Detections", 0) or 0)
                except ValueError:
                    evidence = allele = 0
                # Family-only labels become noisy when many homologues hit;
                # keep them only when several independent tools support them.
                if status == "FAMILY_DETECTED" and evidence < 3:
                    continue
                if evidence == 0 and allele == 0:
                    continue
                candidates.append({
                    "sample_id": sample_dir.name,
                    "gene": row.get("Gene", ""),
                    "priority": rank[status],
                    "evidence_status": status,
                    "evidence_count": evidence,
                    "candidate_allele_count": allele,
                    "total_detections": detections.get(row.get("Gene", ""), 0),
                    "warnings": row.get("Evidence_Warnings", ""),
                    "recommendation": (
                        "reconstruct and validate" if rank[status] <= 3
                        else "reconstruct as an ambiguity or failure challenge"
                    ),
                })
    return sorted(candidates, key=lambda row: (
        row["priority"], -row["evidence_count"], -row["candidate_allele_count"],
        row["sample_id"], row["gene"],
    ))


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare GeneFior reconstructions with an assembly or database reference."
    )
    parser.add_argument("--reconstruction-dir", required=True, type=Path,
                        help="Root containing sample/gene/tool reconstruction outputs.")
    parser.add_argument("--comparison-mode", choices=("assembly", "database"), default="assembly",
                        help="Compare to matched sample assemblies (default) or directly to reference sequences.")
    parser.add_argument("--truth-root", type=Path,
                        help="Dataset root containing long_reads/ and enhanced_contigs/ (assembly mode).")
    parser.add_argument("--reference-fasta", type=Path,
                        help="Gene reference FASTA used to locate the locus in assembly mode.")
    parser.add_argument("--database-fasta", type=Path,
                        help="Original nucleotide reference FASTA used by GeneFior; direct comparison target in database mode.")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample-id", help="Evaluate one sample; default discovers all sample directories.")
    parser.add_argument("--amrfior-root", type=Path,
                        help="Optional root with per-sample AMRFior output directories; emits a live candidate-priority table.")
    args = parser.parse_args(argv)
    if PairwiseAligner is None:
        parser.error("Architect alignment requires Biopython. Install GeneFior with the 'architect' extra.")
    if not args.reconstruction_dir.is_dir():
        parser.error(f"reconstruction directory does not exist: {args.reconstruction_dir}")
    if args.comparison_mode == "assembly":
        if args.truth_root is None or not args.truth_root.is_dir():
            parser.error("assembly mode requires --truth-root pointing to a dataset with sample assemblies")
        if args.reference_fasta is None or not args.reference_fasta.is_file():
            parser.error("assembly mode requires --reference-fasta to locate each gene in the assemblies")
        if args.database_fasta is not None:
            parser.error("--database-fasta is only used in database mode")
        comparison_fasta = args.reference_fasta
    else:
        if args.truth_root is not None or args.reference_fasta is not None:
            parser.error("database mode uses --database-fasta directly; omit --truth-root and --reference-fasta")
        if args.database_fasta is None or not args.database_fasta.is_file():
            parser.error("database mode requires --database-fasta pointing to the original nucleotide FASTA used by GeneFior")
        comparison_fasta = args.database_fasta
    args.output_dir.mkdir(parents=True, exist_ok=True)
    aligner = make_aligner()
    discovered = discover_reconstructions(args.reconstruction_dir, args.sample_id)
    if not discovered:
        parser.error("No *_consensus.fasta files found under sample/gene/tool reconstruction directories.")

    all_rows = []
    missing_assemblies = []
    assembly_cache = {}
    locus_cache = {}
    for sample_id, gene, path in discovered:
        if not sample_id:
            sample_id = args.sample_id or ""
        if not sample_id:
            missing_assemblies.append({"reconstruction": str(path), "reason": "sample id not inferred; use --sample-id"})
            continue
        if args.comparison_mode == "database":
            all_rows.extend(evaluate_database_file(path, gene, sample_id,
                                                   comparison_fasta, aligner))
        else:
            assemblies = find_assemblies(args.truth_root, sample_id)
            if not assemblies:
                missing_assemblies.append({"sample_id": sample_id, "reconstruction": str(path), "reason": "no matching truth assembly found"})
                continue
            all_rows.extend(evaluate_file(path, gene, sample_id, comparison_fasta, assemblies,
                                          aligner, assembly_cache, locus_cache))

    tsv_path = args.output_dir / "architect_evaluation.tsv"
    fields = list(all_rows[0]) if all_rows else ["sample_id", "gene", "status", "explanation", "recommended_action"]
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)
    summary = {
        "architect_version": ARCHITECT_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reconstruction_root": str(args.reconstruction_dir),
        "comparison_mode": args.comparison_mode,
        "truth_root": str(args.truth_root) if args.truth_root else None,
        "reference_fasta": str(args.reference_fasta) if args.reference_fasta else None,
        "database_fasta": str(args.database_fasta) if args.database_fasta else None,
        "comparison_fasta": str(comparison_fasta),
        "evaluated_rows": len(all_rows),
        "status_counts": dict(Counter(row["status"] for row in all_rows)),
        "missing_assemblies": missing_assemblies,
        "evaluations": all_rows,
    }
    json_path = args.output_dir / "architect_evaluation.json"
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Evaluated {len(all_rows)} reconstruction/{args.comparison_mode} comparisons.")
    print("Statuses:", ", ".join(f"{k}={v}" for k, v in sorted(summary["status_counts"].items())) or "none")
    print(f"TSV: {tsv_path}")
    print(f"JSON: {json_path}")
    if missing_assemblies:
        print(f"Skipped {len(missing_assemblies)} item(s) without a matched assembly or inferred sample id.")
    if args.amrfior_root:
        candidates = discover_candidates(args.amrfior_root)
        candidate_path = args.output_dir / "architect_candidates.tsv"
        candidate_fields = ["sample_id", "gene", "priority", "evidence_status", "evidence_count",
                            "candidate_allele_count", "total_detections", "warnings", "recommendation"]
        with candidate_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=candidate_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(candidates)
        print(f"Candidate genes: {len(candidates)} in {candidate_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
