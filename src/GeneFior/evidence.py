from dataclasses import dataclass, field
import re
from typing import Dict, Iterable, List, Tuple


EXACT_ALLELE_DETECTED = "EXACT_ALLELE_DETECTED"
CANDIDATE_ALLELE_DETECTED = "CANDIDATE_ALLELE_DETECTED"
ALLELE_LIKE = "ALLELE_LIKE"
FAMILY_DETECTED = "FAMILY_DETECTED"
PARTIAL_OR_DIVERGENT = "PARTIAL_OR_DIVERGENT"
MIXED_OR_MOSAIC = "MIXED_OR_MOSAIC"
PROFILE_DETECTED = "PROFILE_DETECTED"
NOT_DETECTED = "NOT_DETECTED"
MUST_FLAG_REVIEW = "MUST_FLAG_REVIEW"
DETECTED = "DETECTED"
# Backwards-compatible import alias. Legacy-facing reports use DETECTED.
LEGACY_RELAXED_DETECTED = DETECTED

DETECTION_SYSTEM_QUALIFIED = "qualified"
DETECTION_SYSTEM_LEGACY_RELAXED = "legacy-relaxed"
DETECTION_SYSTEMS = {
    DETECTION_SYSTEM_QUALIFIED,
    DETECTION_SYSTEM_LEGACY_RELAXED,
}

POSITIVE_EXACT_STATUSES = {EXACT_ALLELE_DETECTED, PROFILE_DETECTED}
EVIDENCE_PRESENT_STATUSES = {
    EXACT_ALLELE_DETECTED,
    CANDIDATE_ALLELE_DETECTED,
    ALLELE_LIKE,
    FAMILY_DETECTED,
    MIXED_OR_MOSAIC,
    PROFILE_DETECTED,
    LEGACY_RELAXED_DETECTED,
}


@dataclass
class EvidenceConfig:
    corroborating_depth: int = 2
    exact_identity_min: float = 100.0
    candidate_depth: int = 3
    candidate_identity_min: float = 98.0
    max_internal_gap_bp: int = 15
    max_internal_gap_fraction: float = 0.02
    min_unique_best_reads: int = 2
    min_unique_best_fraction: float = 0.10
    ambiguity_warning_fraction: float = 0.50
    score_tie_absolute: float = 1.0
    score_tie_fraction: float = 0.005
    partial_min_coverage: float = 20.0
    uneven_coverage_cv: float = 1.50


@dataclass
class EvidenceCall:
    gene: str
    family: str
    tool: str
    status: str
    exact_detected: bool = False
    evidence_present: bool = False
    warnings: List[str] = field(default_factory=list)
    competing_alleles: List[str] = field(default_factory=list)
    best_allele: str = ""
    rationale: str = ""

    def warning_text(self) -> str:
        return ";".join(self.warnings)

    @property
    def exact_allele_detected(self) -> bool:
        return self.status == EXACT_ALLELE_DETECTED

    @property
    def candidate_allele_detected(self) -> bool:
        return self.status == CANDIDATE_ALLELE_DETECTED

    @property
    def profile_detected(self) -> bool:
        return self.status == PROFILE_DETECTED


def infer_gene_family(gene: str) -> str:
    """Return a conservative family label from common AMR allele identifiers."""
    family = gene.split("|")[-1].strip()
    family = re.sub(r"_[0-9]+_[A-Z][A-Z0-9.]*$", "", family)
    family = re.sub(r"_[A-Z]{1,6}[0-9]+(?:\.[0-9]+)?$", "", family)
    if family.lower().startswith("bla"):
        family = re.sub(r"-[0-9]+$", "", family)
    return family or gene


def tool_modality(tool_name: str) -> str:
    name = tool_name.lower()
    if "hmmer" in name:
        return "profile"
    if any(token in name for token in ("blastx", "blastp", "diamond")):
        return "protein"
    return "nucleotide"


def scores_tied(best_score: float, candidate_score: float,
                config: EvidenceConfig) -> bool:
    tolerance = max(
        config.score_tie_absolute,
        abs(best_score) * config.score_tie_fraction,
    )
    return abs(best_score - candidate_score) <= tolerance


def normalise_detection_system(value: str) -> str:
    """Normalise public and backwards-compatible detection-system names."""
    name = str(value or DETECTION_SYSTEM_QUALIFIED).strip().lower()
    aliases = {
        "new": DETECTION_SYSTEM_QUALIFIED,
        "evidence": DETECTION_SYSTEM_QUALIFIED,
        "legacy": DETECTION_SYSTEM_LEGACY_RELAXED,
        "relaxed": DETECTION_SYSTEM_LEGACY_RELAXED,
    }
    name = aliases.get(name, name)
    if name not in DETECTION_SYSTEMS:
        raise ValueError(
            f"Unsupported detection system '{value}'. "
            f"Choose from: {', '.join(sorted(DETECTION_SYSTEMS))}"
        )
    return name


def legacy_thresholds_pass(stats, detection_min_coverage: float,
                           detection_min_identity: float,
                           detection_min_base_depth: float,
                           detection_min_num_reads: int,
                           profile: bool = False) -> bool:
    """Apply the original direct gene-level threshold rule."""
    identity_pass = profile or stats.avg_identity >= detection_min_identity
    return (
        stats.gene_coverage >= detection_min_coverage
        and identity_pass
        and stats.base_depth_hit >= detection_min_base_depth
        and stats.num_sequences >= detection_min_num_reads
    )


def _qualified_thresholds_pass(stats, detection_min_coverage: float,
                               detection_min_identity: float,
                               detection_min_base_depth: float,
                               detection_min_num_reads: int,
                               profile: bool = False) -> bool:
    """Apply user thresholds using distinct passing-read support where known."""
    identity_pass = profile or stats.avg_identity >= detection_min_identity
    sequence_support = (
        stats.passing_read_support
        if stats.passing_read_support > 0
        else stats.num_sequences
    )
    return (
        stats.gene_coverage >= detection_min_coverage
        and identity_pass
        and stats.base_depth_hit >= detection_min_base_depth
        and sequence_support >= detection_min_num_reads
    )


def classify_gene_legacy(
        gene: str,
        tool_name: str,
        stats,
        detection_min_coverage: float,
        detection_min_identity: float,
        detection_min_base_depth: float,
        detection_min_num_reads: int) -> EvidenceCall:
    """Reproduce the original relaxed binary detector without allele claims."""
    modality = tool_modality(tool_name)
    detected = legacy_thresholds_pass(
        stats,
        detection_min_coverage,
        detection_min_identity,
        detection_min_base_depth,
        detection_min_num_reads,
        profile=modality == "profile",
    )
    return EvidenceCall(
        gene=gene,
        family=infer_gene_family(gene),
        tool=tool_name,
        status=LEGACY_RELAXED_DETECTED if detected else NOT_DETECTED,
        exact_detected=False,
        evidence_present=detected,
        warnings=["LEGACY_RELAXED_RULES"],
        best_allele=gene,
        rationale=(
            f"legacy coverage={stats.gene_coverage:.2f}%; "
            f"identity={stats.avg_identity:.2f}%; "
            f"covered_depth={stats.base_depth_hit:.2f}; "
            f"sequences={stats.num_sequences}"
        ),
    )


def classify_gene_evidence(
        gene: str,
        tool_name: str,
        stats,
        config: EvidenceConfig,
        detection_min_coverage: float,
        detection_min_identity: float,
        detection_min_base_depth: float,
        detection_min_num_reads: int,
        reads_mode: bool = True,
        must_flag_override: bool = False) -> EvidenceCall:
    """Classify one gene without yet resolving competition inside its family."""
    family = infer_gene_family(gene)
    modality = tool_modality(tool_name)
    is_profile = modality == "profile"
    warnings: List[str] = []

    base_pass = _qualified_thresholds_pass(
        stats,
        detection_min_coverage,
        detection_min_identity,
        detection_min_base_depth,
        detection_min_num_reads,
        profile=is_profile,
    )

    corroborating_coverage = stats.coverage_at_depth(config.corroborating_depth)
    robust_pass = (
        not reads_mode
        or corroborating_coverage >= detection_min_coverage
    )
    max_allowed_gap = max(
        config.max_internal_gap_bp,
        int(round(stats.gene_length * config.max_internal_gap_fraction)),
    )
    discontinuous = (
        stats.longest_internal_gap > max_allowed_gap
        or stats.num_internal_gaps > 1
    )

    passing_support = max(
        stats.passing_read_support,
        stats.num_sequences,
    )
    unique_support = stats.unique_best_read_support
    ambiguous_support = stats.ambiguous_best_read_support
    best_support = max(stats.best_read_support, unique_support + ambiguous_support)
    ambiguity_fraction = (
        ambiguous_support / best_support if best_support > 0 else 0.0
    )
    unique_fraction = (
        unique_support / passing_support if passing_support > 0 else 0.0
    )
    clear_best_support = (
        not reads_mode
        or (
            unique_support >= config.min_unique_best_reads
            and unique_fraction >= config.min_unique_best_fraction
        )
    )

    if stats.gene_coverage < detection_min_coverage:
        warnings.append("PARTIAL_COVERAGE")
    if stats.avg_identity < detection_min_identity and not is_profile:
        warnings.append("LOW_MEAN_IDENTITY")
    if discontinuous:
        warnings.append("DISCONTINUOUS_COVERAGE")
    if reads_mode and not robust_pass:
        warnings.append("LOW_CORROBORATED_COVERAGE")
        if stats.gene_coverage >= detection_min_coverage:
            warnings.append("THRESHOLD_SENSITIVE")
            warnings.append("SINGLE_READ_BRIDGE_RISK")
    if stats.depth_cv > config.uneven_coverage_cv:
        warnings.append("UNEVEN_COVERAGE")
    if ambiguity_fraction >= config.ambiguity_warning_fraction:
        warnings.append("MULTI_ALLELE_READS")
    if reads_mode and modality == "nucleotide":
        if unique_support < config.min_unique_best_reads:
            warnings.append("INSUFFICIENT_UNIQUE_SUPPORT")
        if unique_fraction < config.min_unique_best_fraction:
            warnings.append("LOW_UNIQUE_SUPPORT_FRACTION")
    if modality == "protein":
        warnings.append("PROTEIN_ONLY_ALLELE_AMBIGUITY")

    # Full-length gene FASTA inputs are not read-depth evidence. A single
    # full-length nucleotide query can support a candidate/exact database
    # allele, while read inputs keep the configured per-base depth requirement.
    candidate_depth = 1 if not reads_mode else max(1, int(config.candidate_depth))
    candidate_coverage = stats.coverage_at_depth(candidate_depth)
    candidate_coverage_pass = candidate_coverage >= 100.0 - 1e-9
    candidate_identity_pass = (
        stats.avg_identity >= config.candidate_identity_min
    )
    candidate_pass = (
        modality == "nucleotide"
        and base_pass
        and candidate_coverage_pass
        and candidate_identity_pass
        and clear_best_support
        and not discontinuous
    )

    exact_coverage_pass = stats.gene_coverage >= 100.0 - 1e-9
    exact_identity_pass = stats.avg_identity >= 100.0 - 1e-9
    exact_corroboration_pass = (
        candidate_coverage_pass
        and (
            not reads_mode
            or stats.coverage_at_depth(
                max(config.corroborating_depth, candidate_depth)
            ) >= 100.0 - 1e-9
        )
    )
    if modality == "nucleotide" and base_pass:
        if not candidate_coverage_pass:
            warnings.append(
                "INCOMPLETE_CANDIDATE_ALLELE_COVERAGE"
                if not reads_mode
                else "INCOMPLETE_CANDIDATE_ALLELE_3X_COVERAGE"
            )
        if not candidate_identity_pass:
            warnings.append("LOW_CANDIDATE_ALLELE_IDENTITY")
        if not exact_coverage_pass:
            warnings.append("INCOMPLETE_EXACT_ALLELE_COVERAGE")
        if not exact_identity_pass:
            warnings.append("NON_IDENTICAL_ALLELE_SEQUENCE")
        if not exact_corroboration_pass:
            warnings.append("INCOMPLETE_EXACT_CORROBORATION")

    if must_flag_override and not base_pass:
        status = MUST_FLAG_REVIEW
        warnings.append("MUST_FLAG_THRESHOLD_OVERRIDE")
    elif not base_pass:
        if stats.num_sequences > 0 and stats.gene_coverage >= config.partial_min_coverage:
            status = PARTIAL_OR_DIVERGENT
        else:
            status = NOT_DETECTED
    elif discontinuous or (reads_mode and not robust_pass):
        status = MIXED_OR_MOSAIC
    elif modality == "profile":
        status = PROFILE_DETECTED
    elif modality == "protein":
        if ambiguity_fraction >= config.ambiguity_warning_fraction:
            status = FAMILY_DETECTED
        else:
            status = ALLELE_LIKE
    elif candidate_pass and (
            exact_coverage_pass
            and exact_identity_pass
            and exact_corroboration_pass):
        status = EXACT_ALLELE_DETECTED
    elif candidate_pass:
        status = CANDIDATE_ALLELE_DETECTED
    elif reads_mode and not clear_best_support:
        status = FAMILY_DETECTED
    else:
        status = ALLELE_LIKE

    exact_detected = status in POSITIVE_EXACT_STATUSES
    # Evidence is the original high-confidence detection tier: all configured
    # user thresholds must pass. Diagnostic partial or must-flag review states
    # remain visible, but are not positive evidence.
    evidence_present = base_pass
    rationale = (
        f"coverage={stats.gene_coverage:.2f}%; "
        f"coverage_{config.corroborating_depth}x={corroborating_coverage:.2f}%; "
        f"coverage_{candidate_depth}x={candidate_coverage:.2f}%; "
        f"identity={stats.avg_identity:.2f}%; "
        f"unique_best_reads={unique_support}; "
        f"ambiguous_best_reads={ambiguous_support}"
    )
    return EvidenceCall(
        gene=gene,
        family=family,
        tool=tool_name,
        status=status,
        exact_detected=exact_detected,
        evidence_present=evidence_present,
        warnings=sorted(set(warnings)),
        best_allele=gene,
        rationale=rationale,
    )


def resolve_family_calls(
        calls: Dict[str, EvidenceCall],
        stats_by_gene: Dict[str, object],
        config: EvidenceConfig) -> Dict[str, dict]:
    """Resolve competing alleles and return one summary row per family."""
    grouped: Dict[str, List[str]] = {}
    for gene, call in calls.items():
        if call.evidence_present:
            grouped.setdefault(call.family, []).append(gene)

    summaries: Dict[str, dict] = {}
    for family, genes in grouped.items():
        ranked = sorted(
            genes,
            key=lambda gene: (
                stats_by_gene[gene].unique_best_read_support,
                stats_by_gene[gene].best_read_support,
                stats_by_gene[gene].coverage_at_depth(config.corroborating_depth),
                stats_by_gene[gene].gene_coverage,
                stats_by_gene[gene].avg_identity,
            ),
            reverse=True,
        )
        best_gene = ranked[0]
        competitors = ranked[1:]
        best_stats = stats_by_gene[best_gene]

        unresolved = False
        if competitors:
            second_stats = stats_by_gene[competitors[0]]
            unresolved = (
                best_stats.unique_best_read_support
                <= second_stats.unique_best_read_support
                or best_stats.unique_best_read_support < config.min_unique_best_reads
            )

        if unresolved:
            for gene in ranked:
                call = calls[gene]
                call.competing_alleles = [other for other in ranked if other != gene]
                call.best_allele = best_gene
                if call.status in (
                        EXACT_ALLELE_DETECTED,
                        CANDIDATE_ALLELE_DETECTED):
                    call.status = FAMILY_DETECTED
                    call.exact_detected = False
                if "COMPETING_ALLELES" not in call.warnings:
                    call.warnings.append("COMPETING_ALLELES")
                    call.warnings.sort()

        best_call = calls[best_gene]
        summaries[family] = {
            "family": family,
            "status": best_call.status,
            "best_allele": best_gene,
            "competing_alleles": competitors,
            "exact_allele_resolved": (
                best_call.exact_allele_detected and not unresolved
            ),
            "candidate_allele_resolved": (
                best_call.candidate_allele_detected and not unresolved
            ),
            "warnings": sorted(set(best_call.warnings)),
        }
    return summaries


def status_rank(status: str) -> int:
    order = {
        EXACT_ALLELE_DETECTED: 8,
        PROFILE_DETECTED: 7,
        LEGACY_RELAXED_DETECTED: 7,
        CANDIDATE_ALLELE_DETECTED: 6,
        ALLELE_LIKE: 5,
        FAMILY_DETECTED: 4,
        MIXED_OR_MOSAIC: 3,
        PARTIAL_OR_DIVERGENT: 2,
        MUST_FLAG_REVIEW: 1,
        NOT_DETECTED: 0,
    }
    return order.get(status, 0)


def best_status(statuses: Iterable[str]) -> str:
    values = list(statuses)
    return max(values, key=status_rank) if values else NOT_DETECTED
