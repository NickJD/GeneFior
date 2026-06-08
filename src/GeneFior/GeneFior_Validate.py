"""
GeneFior_Validate.py — Reconstruction Validation Framework
===========================================================

Performs a suite of biological-plausibility and artifact-detection checks on
a reconstructed gene sequence produced by GeneFior-Reconstruct.  The goal is
to flag sequences that may represent assembly artefacts or chimeras rather
than real biological variation.

Checks performed
----------------
1. **ORF integrity** (coding sequences only)
   Searches all six reading frames for the longest uninterrupted ORF,
   counts in-frame stop codons, and scores completeness against the
   expected gene length.

2. **GC content**
   Flags GC% that falls outside the normal biological range (30–70 %) or
   deviates strongly from the reference sequence.

3. **Repeat structure**
   Detects tandem-repeat units (≥ 10 bp, ≥ 5 copies) that may indicate
   polymerase slippage or assembly collapse.

4. **Coverage uniformity**
   Computes the coefficient of variation (CV) of per-position depth and
   identifies sudden coverage discontinuities that could signal chimeric
   assemblies.

5. **Strand bias**
   Reports the forward / reverse read ratio.  Severe imbalance may
   indicate strand-specific artefacts or alignment to the wrong strand.

6. **Reference identity**
   Compares the consensus to the reference (if provided) and warns when
   identity is low enough to suggest the wrong gene or contamination.

7. **Variant quality**
   Checks that called variants are well-supported (depth ≥ 10, frequency
   ≥ 70 %) and flags positions with low evidence.

Each check returns a :class:`ValidationResult` with a PASS / WARNING / FAIL
status, a confidence score (0–100), a human-readable message, and optional
recommendations.

Usage
-----
::

    from GeneFior_Validate import validate_reconstruction

    report = validate_reconstruction(
        consensus  = reconstructed_seq,
        reference  = reference_seq,      # optional
        per_pos    = per_position_data,
        alignments = alignment_list,
        variants   = variant_list,
        gene_name  = "tet(Q)_1_L33696",
        gene_type  = "coding",
    )
    print(report.summary_report())
"""

import re
import math
import statistics
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from enum import Enum


class ValidationStatus(Enum):
    """Validation result status levels"""
    PASS = "PASS"
    WARNING = "WARNING"
    FAIL = "FAIL"
    UNKNOWN = "UNKNOWN"


@dataclass
class ValidationResult:
    """Result of a single validation check"""
    name: str
    status: ValidationStatus
    score: float  # 0-100, higher is better
    message: str
    details: Dict = field(default_factory=dict)
    recommendations: List[str] = field(default_factory=list)


@dataclass
class ValidationReport:
    """Complete validation report for a reconstruction"""
    gene_name: str
    consensus_length: int
    results: List[ValidationResult] = field(default_factory=list)
    overall_status: ValidationStatus = ValidationStatus.UNKNOWN
    overall_score: float = 0.0

    def add_result(self, result: ValidationResult):
        """Add a validation result to the report"""
        self.results.append(result)
        self._update_overall()

    def _update_overall(self):
        """Update overall status and score based on all results"""
        if not self.results:
            self.overall_status = ValidationStatus.UNKNOWN
            self.overall_score = 0.0
            return

        # Overall score is weighted average
        self.overall_score = sum(r.score for r in self.results) / len(self.results)

        # Overall status is worst case
        statuses = [r.status for r in self.results]
        if ValidationStatus.FAIL in statuses:
            self.overall_status = ValidationStatus.FAIL
        elif ValidationStatus.WARNING in statuses:
            self.overall_status = ValidationStatus.WARNING
        else:
            self.overall_status = ValidationStatus.PASS

    def summary_report(self) -> str:
        """Generate human-readable summary report"""
        lines = []
        lines.append("=" * 70)
        lines.append("RECONSTRUCTION VALIDATION REPORT")
        lines.append("=" * 70)
        lines.append(f"Gene: {self.gene_name}")
        lines.append(f"Length: {self.consensus_length} bp")
        lines.append(f"Overall Status: {self.overall_status.value}")
        lines.append(f"Overall Score: {self.overall_score:.1f}/100")
        lines.append("=" * 70)
        lines.append("")

        for result in self.results:
            symbol = {
                ValidationStatus.PASS: "✓",
                ValidationStatus.WARNING: "⚠",
                ValidationStatus.FAIL: "✗",
                ValidationStatus.UNKNOWN: "?"
            }[result.status]

            lines.append(f"{symbol} {result.name}: {result.status.value} ({result.score:.0f}/100)")
            lines.append(f"  {result.message}")

            if result.details:
                for key, val in result.details.items():
                    lines.append(f"    • {key}: {val}")

            if result.recommendations:
                lines.append("  Recommendations:")
                for rec in result.recommendations:
                    lines.append(f"    → {rec}")
            lines.append("")

        return "\n".join(lines)


class ReconstructionValidator:
    """
    Validation engine for a reconstructed gene sequence.

    Parameters
    ----------
    consensus : str
        The reconstructed consensus sequence (will be upper-cased internally).
    reference : str, optional
        Reference sequence from the database or MD-tag reconstruction.
    per_pos : list of dict, optional
        Per-position depth and base-frequency data as produced by
        ``_call_consensus`` in GeneFior_Reconstruct.
    alignments : list of dict, optional
        SAM-derived alignment records for the gene.
    variants : list of dict, optional
        Called variants from ``_compute_variants``.
    gene_name : str
        Display name used in the validation report.
    gene_type : str
        ``"coding"`` (default) or another type such as ``"rRNA"``.
        ORF-integrity checks are only run for coding genes.
    """

    def __init__(
        self,
        consensus: str,
        reference: Optional[str] = None,
        per_pos: Optional[List[dict]] = None,
        alignments: Optional[List[dict]] = None,
        variants: Optional[List[dict]] = None,
        gene_name: str = "unknown",
        gene_type: str = "coding",
    ):
        self.consensus  = consensus.upper()
        self.reference  = reference.upper() if reference else None
        self.per_pos    = per_pos or []
        self.alignments = alignments or []
        self.variants   = variants or []
        self.gene_name  = gene_name
        self.gene_type  = gene_type

        self.report = ValidationReport(
            gene_name=gene_name,
            consensus_length=len(consensus),
        )

    # ─────────────────────────────────────────────────────────────────────
    # 1. Biological plausibility
    # ─────────────────────────────────────────────────────────────────────

    def validate_orf_integrity(self) -> ValidationResult:
        """
        Validate open reading frame integrity for coding sequences.
        Checks for internal stop codons, frameshift mutations, etc.
        """
        if self.gene_type != "coding":
            return ValidationResult(
                name="ORF Integrity",
                status=ValidationStatus.PASS,
                score=100.0,
                message="Not applicable for non-coding sequences"
            )

        # Remove N's for ORF analysis
        clean_seq = self.consensus.replace("N", "")

        if len(clean_seq) < 30:
            return ValidationResult(
                name="ORF Integrity",
                status=ValidationStatus.WARNING,
                score=50.0,
                message="Sequence too short for reliable ORF analysis",
                details={"clean_length": len(clean_seq)}
            )

        # Find longest ORF in all 6 frames
        best_orf_len = 0
        best_frame = 0
        internal_stops = 0

        for frame in range(3):
            length, stops = self._find_orf_in_frame(clean_seq[frame:])
            if length > best_orf_len:
                best_orf_len = length
                best_frame = frame
                internal_stops = stops

        # Check reverse complement
        rev_comp = self._reverse_complement(clean_seq)
        for frame in range(3):
            length, stops = self._find_orf_in_frame(rev_comp[frame:])
            if length > best_orf_len:
                best_orf_len = length
                best_frame = -frame - 1
                internal_stops = stops

        # Scoring
        expected_length = len(clean_seq) - (len(clean_seq) % 3)
        coverage = best_orf_len / expected_length if expected_length > 0 else 0

        if coverage >= 0.95 and internal_stops == 0:
            status = ValidationStatus.PASS
            score = 100.0
            message = f"Complete ORF found (frame {best_frame}, {best_orf_len} bp)"
        elif coverage >= 0.80 and internal_stops <= 1:
            status = ValidationStatus.WARNING
            score = 75.0
            message = f"ORF mostly intact with minor issues (coverage {coverage*100:.1f}%)"
        elif coverage >= 0.50:
            status = ValidationStatus.WARNING
            score = 50.0
            message = f"Partial ORF detected ({internal_stops} internal stops)"
        else:
            status = ValidationStatus.FAIL
            score = 25.0
            message = f"ORF severely disrupted (coverage {coverage*100:.1f}%, {internal_stops} stops)"

        recommendations = []
        if internal_stops > 0:
            recommendations.append("Check for sequencing errors at stop codon positions")
        if coverage < 0.95:
            recommendations.append("Verify gene boundaries and check for frameshifts")
        if best_orf_len < 300:
            recommendations.append("Gene may be truncated or annotation may be incorrect")

        return ValidationResult(
            name="ORF Integrity",
            status=status,
            score=score,
            message=message,
            details={
                "orf_length": best_orf_len,
                "frame": best_frame,
                "coverage": f"{coverage*100:.1f}%",
                "internal_stops": internal_stops
            },
            recommendations=recommendations
        )

    def validate_gc_content(self) -> ValidationResult:
        """
        Validate GC content is within reasonable biological range.
        Extreme GC content may indicate contamination or artifact.
        """
        clean_seq = self.consensus.replace("N", "")

        if len(clean_seq) < 100:
            return ValidationResult(
                name="GC Content",
                status=ValidationStatus.WARNING,
                score=50.0,
                message="Sequence too short for GC analysis"
            )

        gc_count = clean_seq.count("G") + clean_seq.count("C")
        gc_pct = (gc_count / len(clean_seq)) * 100

        # Compare to reference if available
        ref_gc_pct = None
        gc_diff = None
        if self.reference:
            ref_clean = self.reference.replace("N", "")
            ref_gc = ref_clean.count("G") + ref_clean.count("C")
            ref_gc_pct = (ref_gc / len(ref_clean)) * 100 if ref_clean else None
            gc_diff = abs(gc_pct - ref_gc_pct) if ref_gc_pct else None

        # Biological range: most bacterial genes 30-70% GC
        # Extremes: 25-75% possible but unusual
        if 30 <= gc_pct <= 70:
            status = ValidationStatus.PASS
            score = 100.0
            message = f"GC content within normal range ({gc_pct:.1f}%)"
        elif 25 <= gc_pct <= 75:
            status = ValidationStatus.WARNING
            score = 75.0
            message = f"GC content unusual but possible ({gc_pct:.1f}%)"
        else:
            status = ValidationStatus.FAIL
            score = 25.0
            message = f"GC content highly unusual ({gc_pct:.1f}%) - possible artifact"

        # Large deviation from reference is suspicious
        if gc_diff and gc_diff > 10:
            status = ValidationStatus.WARNING
            score = min(score, 60.0)
            message += f" - differs significantly from reference ({ref_gc_pct:.1f}%)"

        details = {"gc_percent": f"{gc_pct:.1f}%"}
        if ref_gc_pct:
            details["reference_gc"] = f"{ref_gc_pct:.1f}%"
            details["difference"] = f"{gc_diff:.1f}%"

        recommendations = []
        if gc_pct < 25 or gc_pct > 75:
            recommendations.append("Check for contamination or adapter sequences")
            recommendations.append("Verify gene identity with BLAST")
        if gc_diff and gc_diff > 15:
            recommendations.append("Large GC deviation from reference - verify sequence identity")

        return ValidationResult(
            name="GC Content",
            status=status,
            score=score,
            message=message,
            details=details,
            recommendations=recommendations
        )

    def validate_repeat_structure(self) -> ValidationResult:
        """
        Detect unusual repeat structures that may indicate assembly artifacts.
        """
        # Look for tandem repeats (same sequence repeated multiple times)
        max_repeat_len = min(50, len(self.consensus) // 4)
        suspicious_repeats = []

        for repeat_len in range(10, max_repeat_len + 1):
            for i in range(len(self.consensus) - repeat_len * 2):
                unit = self.consensus[i:i + repeat_len]
                if "N" in unit:
                    continue

                # Count consecutive repeats
                repeat_count = 1
                pos = i + repeat_len
                while pos + repeat_len <= len(self.consensus):
                    if self.consensus[pos:pos + repeat_len] == unit:
                        repeat_count += 1
                        pos += repeat_len
                    else:
                        break

                # 5+ repeats of 10+bp is suspicious unless it's a known repeat region
                if repeat_count >= 5 and repeat_len >= 10:
                    suspicious_repeats.append({
                        "unit": unit[:20] + "..." if len(unit) > 20 else unit,
                        "length": repeat_len,
                        "count": repeat_count,
                        "position": i
                    })

        if not suspicious_repeats:
            return ValidationResult(
                name="Repeat Structure",
                status=ValidationStatus.PASS,
                score=100.0,
                message="No unusual repeat structures detected"
            )
        elif len(suspicious_repeats) == 1 and suspicious_repeats[0]["count"] <= 7:
            return ValidationResult(
                name="Repeat Structure",
                status=ValidationStatus.WARNING,
                score=75.0,
                message=f"Minor tandem repeat detected ({suspicious_repeats[0]['count']}× {suspicious_repeats[0]['length']}bp)",
                details={"repeats": str(suspicious_repeats[0])},
                recommendations=["Verify repeat is biological (e.g., signal peptide)"]
            )
        else:
            return ValidationResult(
                name="Repeat Structure",
                status=ValidationStatus.FAIL,
                score=40.0,
                message=f"Multiple suspicious repeat structures detected ({len(suspicious_repeats)} regions)",
                details={"repeat_count": len(suspicious_repeats), "examples": str(suspicious_repeats[:3])},
                recommendations=[
                    "Possible polymerase slippage or assembly artifact",
                    "Verify with independent sequencing",
                    "Check reference for known repeat regions"
                ]
            )

    # ─────────────────────────────────────────────────────────────────────
    # 2. Assembly artifact detection
    # ─────────────────────────────────────────────────────────────────────

    def validate_coverage_uniformity(self) -> ValidationResult:
        """
        Validate that coverage is uniform across the gene.
        Sudden drops or discontinuities suggest chimeras or misassembly.
        """
        if not self.per_pos:
            return ValidationResult(
                name="Coverage Uniformity",
                status=ValidationStatus.UNKNOWN,
                score=50.0,
                message="No coverage data available"
            )

        depths = [p["depth"] for p in self.per_pos]

        if not depths or all(d == 0 for d in depths):
            return ValidationResult(
                name="Coverage Uniformity",
                status=ValidationStatus.FAIL,
                score=0.0,
                message="No coverage detected"
            )

        # Calculate statistics
        non_zero = [d for d in depths if d > 0]
        mean_depth = statistics.mean(non_zero) if non_zero else 0
        median_depth = statistics.median(non_zero) if non_zero else 0

        if len(non_zero) > 1:
            stdev_depth = statistics.stdev(non_zero)
            cv = (stdev_depth / mean_depth * 100) if mean_depth > 0 else 999
        else:
            stdev_depth = 0
            cv = 0

        # Detect discontinuities (sudden drops)
        discontinuities = []
        window_size = max(10, len(depths) // 20)

        for i in range(0, len(depths) - window_size, window_size):
            window1 = depths[i:i + window_size]
            window2 = depths[i + window_size:i + 2 * window_size]

            avg1 = statistics.mean(window1) if window1 else 0
            avg2 = statistics.mean(window2) if window2 else 0

            if avg1 > 5 and avg2 > 5:  # Both have decent coverage
                ratio = min(avg1, avg2) / max(avg1, avg2) if max(avg1, avg2) > 0 else 1
                if ratio < 0.3:  # 70%+ drop
                    discontinuities.append({
                        "position": i + window_size,
                        "before": avg1,
                        "after": avg2,
                        "drop_pct": (1 - ratio) * 100
                    })

        # Scoring
        if cv < 50 and len(discontinuities) == 0:
            status = ValidationStatus.PASS
            score = 100.0
            message = f"Coverage is uniform (CV={cv:.1f}%)"
        elif cv < 80 and len(discontinuities) <= 1:
            status = ValidationStatus.WARNING
            score = 75.0
            message = f"Coverage moderately variable (CV={cv:.1f}%)"
        elif cv < 120:
            status = ValidationStatus.WARNING
            score = 50.0
            message = f"Coverage highly variable (CV={cv:.1f}%)"
        else:
            status = ValidationStatus.FAIL
            score = 25.0
            message = f"Coverage extremely variable (CV={cv:.1f}%) - possible artifact"

        if discontinuities:
            status = ValidationStatus.FAIL
            score = min(score, 40.0)
            message += f" with {len(discontinuities)} discontinuities"

        recommendations = []
        if cv > 80:
            recommendations.append("High coverage variability suggests uneven amplification or chimera")
        if discontinuities:
            recommendations.append(f"Coverage discontinuities at positions: {[d['position'] for d in discontinuities]}")
            recommendations.append("Check for chimeric assembly or misalignment")

        return ValidationResult(
            name="Coverage Uniformity",
            status=status,
            score=score,
            message=message,
            details={
                "mean_depth": f"{mean_depth:.1f}×",
                "median_depth": f"{median_depth:.1f}×",
                "cv_percent": f"{cv:.1f}%",
                "discontinuities": len(discontinuities)
            },
            recommendations=recommendations
        )

    def validate_strand_bias(self) -> ValidationResult:
        """
        Check for severe strand bias which may indicate artifacts.
        Real genes should have reads from both strands.
        """
        if not self.alignments:
            return ValidationResult(
                name="Strand Bias",
                status=ValidationStatus.UNKNOWN,
                score=50.0,
                message="No alignment data available"
            )

        forward_reads = 0
        reverse_reads = 0

        for aln in self.alignments:
            flag = aln.get("flag", 0)
            if flag & 0x10:  # Reverse strand
                reverse_reads += 1
            else:
                forward_reads += 1

        total = forward_reads + reverse_reads
        if total == 0:
            return ValidationResult(
                name="Strand Bias",
                status=ValidationStatus.FAIL,
                score=0.0,
                message="No reads detected"
            )

        forward_pct = (forward_reads / total) * 100
        reverse_pct = (reverse_reads / total) * 100
        bias_ratio = max(forward_pct, reverse_pct) / min(forward_pct, reverse_pct) if min(forward_pct, reverse_pct) > 0 else 999

        # Scoring
        if 40 <= forward_pct <= 60:  # Near 50/50
            status = ValidationStatus.PASS
            score = 100.0
            message = f"No strand bias detected ({forward_pct:.0f}%/{reverse_pct:.0f}%)"
        elif 30 <= forward_pct <= 70:
            status = ValidationStatus.WARNING
            score = 75.0
            message = f"Minor strand bias ({forward_pct:.0f}%/{reverse_pct:.0f}%)"
        elif 20 <= forward_pct <= 80:
            status = ValidationStatus.WARNING
            score = 50.0
            message = f"Moderate strand bias ({forward_pct:.0f}%/{reverse_pct:.0f}%)"
        else:
            status = ValidationStatus.FAIL
            score = 25.0
            message = f"Severe strand bias ({forward_pct:.0f}%/{reverse_pct:.0f}%) - possible artifact"

        recommendations = []
        if bias_ratio > 3:
            recommendations.append("Severe strand bias may indicate:")
            recommendations.append("  - Strand-specific library prep")
            recommendations.append("  - Systematic sequencing errors")
            recommendations.append("  - Reference sequence on wrong strand")

        return ValidationResult(
            name="Strand Bias",
            status=status,
            score=score,
            message=message,
            details={
                "forward_reads": forward_reads,
                "reverse_reads": reverse_reads,
                "forward_percent": f"{forward_pct:.1f}%",
                "reverse_percent": f"{reverse_pct:.1f}%"
            },
            recommendations=recommendations
        )

    # ─────────────────────────────────────────────────────────────────────
    # 3. Reference consistency
    # ─────────────────────────────────────────────────────────────────────

    def validate_reference_identity(self) -> ValidationResult:
        """
        Validate similarity to reference sequence.
        Too much divergence may indicate wrong gene or contamination.
        """
        if not self.reference:
            return ValidationResult(
                name="Reference Identity",
                status=ValidationStatus.UNKNOWN,
                score=50.0,
                message="No reference sequence provided"
            )

        # Calculate identity
        min_len = min(len(self.consensus), len(self.reference))
        matches = sum(
            1 for i in range(min_len)
            if self.consensus[i] == self.reference[i] and self.consensus[i] != "N"
        )
        valid_pos = sum(
            1 for i in range(min_len)
            if self.consensus[i] != "N"
        )

        identity_pct = (matches / valid_pos * 100) if valid_pos > 0 else 0

        # Score based on identity
        if identity_pct >= 95:
            status = ValidationStatus.PASS
            score = 100.0
            message = f"High identity to reference ({identity_pct:.1f}%)"
        elif identity_pct >= 85:
            status = ValidationStatus.PASS
            score = 90.0
            message = f"Good identity to reference ({identity_pct:.1f}%)"
        elif identity_pct >= 70:
            status = ValidationStatus.WARNING
            score = 70.0
            message = f"Moderate identity to reference ({identity_pct:.1f}%) - verify gene identity"
        elif identity_pct >= 50:
            status = ValidationStatus.WARNING
            score = 50.0
            message = f"Low identity to reference ({identity_pct:.1f}%) - possible different allele"
        else:
            status = ValidationStatus.FAIL
            score = 25.0
            message = f"Very low identity to reference ({identity_pct:.1f}%) - wrong gene or contamination"

        recommendations = []
        if identity_pct < 85:
            recommendations.append("Verify gene identity with BLAST against database")
        if identity_pct < 70:
            recommendations.append("Consider possibility of:")
            recommendations.append("  - Different species/strain")
            recommendations.append("  - Horizontal gene transfer variant")
            recommendations.append("  - Wrong reference gene")
            recommendations.append("  - Contamination")

        return ValidationResult(
            name="Reference Identity",
            status=status,
            score=score,
            message=message,
            details={
                "identity_percent": f"{identity_pct:.1f}%",
                "matches": matches,
                "valid_positions": valid_pos
            },
            recommendations=recommendations
        )

    # ─────────────────────────────────────────────────────────────────────
    # 4. Variant quality
    # ─────────────────────────────────────────────────────────────────────

    def validate_variant_quality(self) -> ValidationResult:
        """
        Validate that called variants are high quality and not artifacts.
        """
        if not self.variants:
            return ValidationResult(
                name="Variant Quality",
                status=ValidationStatus.PASS,
                score=100.0,
                message="No variants called (identical to reference)"
            )

        high_qual_variants = 0
        low_qual_variants = 0
        suspicious_variants = []

        for var in self.variants:
            depth = var.get("depth", 0)
            freq = var.get("freq_pct", 0)
            var_type = var.get("type", "unknown")

            # High quality: depth >= 10, frequency >= 70%
            if depth >= 10 and freq >= 70:
                high_qual_variants += 1
            # Low quality: depth < 5 or frequency < 40%
            elif depth < 5 or freq < 40:
                low_qual_variants += 1
                suspicious_variants.append({
                    "position": var.get("pos"),
                    "type": var_type,
                    "depth": depth,
                    "frequency": freq
                })

        total_variants = len(self.variants)
        high_qual_pct = (high_qual_variants / total_variants * 100) if total_variants > 0 else 0

        if high_qual_pct >= 90 and low_qual_variants == 0:
            status = ValidationStatus.PASS
            score = 100.0
            message = f"All {total_variants} variants are high quality"
        elif high_qual_pct >= 75:
            status = ValidationStatus.WARNING
            score = 80.0
            message = f"{high_qual_variants}/{total_variants} variants high quality"
        elif high_qual_pct >= 50:
            status = ValidationStatus.WARNING
            score = 60.0
            message = f"Only {high_qual_pct:.0f}% variants high quality"
        else:
            status = ValidationStatus.FAIL
            score = 40.0
            message = f"Many low-quality variants ({low_qual_variants}/{total_variants})"

        recommendations = []
        if low_qual_variants > 0:
            recommendations.append(f"{low_qual_variants} low-quality variants detected")
            recommendations.append("Low-quality variants may be sequencing errors")
            recommendations.append("Consider increasing -min_depth threshold")

        return ValidationResult(
            name="Variant Quality",
            status=status,
            score=score,
            message=message,
            details={
                "total_variants": total_variants,
                "high_quality": high_qual_variants,
                "low_quality": low_qual_variants,
                "high_qual_percent": f"{high_qual_pct:.1f}%"
            },
            recommendations=recommendations
        )

    # ─────────────────────────────────────────────────────────────────────
    # Utility methods
    # ─────────────────────────────────────────────────────────────────────

    def _find_orf_in_frame(self, seq: str) -> Tuple[int, int]:
        """
        Scan *seq* in the given reading frame (already offset externally) and
        return ``(longest_segment_nt, internal_stop_count)``.

        ``longest_segment_nt`` is the length of the longest uninterrupted
        coding stretch (i.e. the stretch with the fewest stops on either side).

        ``internal_stop_count`` is the number of *unexpected* stop codons:
        every observed stop minus one (the single legitimate terminal stop).
        For a clean ORF ``(ATG … TAA)`` this returns 0; one premature stop
        returns 1; and so on.
        """
        if len(seq) < 3:
            return 0, 0

        stop_codons = {"TAA", "TAG", "TGA"}
        longest = 0
        current = 0
        total_stops = 0

        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if codon in stop_codons:
                if current > longest:
                    longest = current
                total_stops += 1
                current = 0
            else:
                current += 3

        # Final trailing stretch (no terminal stop — possibly truncated)
        if current > longest:
            longest = current

        # One stop is the expected ORF terminator; every extra stop is internal.
        internal_stops = max(0, total_stops - 1)
        return longest, internal_stops

    def _reverse_complement(self, seq: str) -> str:
        """Return reverse complement of DNA sequence"""
        complement = str.maketrans("ACGT", "TGCA")
        return seq.translate(complement)[::-1]

    # ─────────────────────────────────────────────────────────────────────
    # Main entry point
    # ─────────────────────────────────────────────────────────────────────

    def validate_all(self) -> "ValidationReport":
        """Run all checks and return the complete :class:`ValidationReport`."""

        # Biological plausibility
        self.report.add_result(self.validate_orf_integrity())
        self.report.add_result(self.validate_gc_content())
        self.report.add_result(self.validate_repeat_structure())
        # Assembly artifact detection
        self.report.add_result(self.validate_coverage_uniformity())
        self.report.add_result(self.validate_strand_bias())
        # Reference consistency
        self.report.add_result(self.validate_reference_identity())
        # Variant quality
        self.report.add_result(self.validate_variant_quality())
        return self.report


def validate_reconstruction(
    consensus: str,
    reference: Optional[str],
    per_pos: List[dict],
    alignments: List[dict],
    variants: List[dict],
    gene_name: str,
    gene_type: str = "coding",
) -> "ValidationReport":
    """
    Convenience wrapper: construct a :class:`ReconstructionValidator`, run all
    checks, and return the :class:`ValidationReport`.

    Parameters
    ----------
    consensus : str
        Reconstructed consensus sequence.
    reference : str or None
        Reference sequence (from database or MD-tag reconstruction).
    per_pos : list of dict
        Per-position coverage / base-frequency data.
    alignments : list of dict
        Read-alignment records for the gene.
    variants : list of dict
        Called variant records.
    gene_name : str
        Gene identifier used in the report header.
    gene_type : str
        ``"coding"`` (default) activates ORF-integrity checking.

    Returns
    -------
    ValidationReport
        Complete report with per-check results and overall status / score.
    """
    validator = ReconstructionValidator(
        consensus=consensus,
        reference=reference,
        per_pos=per_pos,
        alignments=alignments,
        variants=variants,
        gene_name=gene_name,
        gene_type=gene_type,
    )

    return validator.validate_all()

