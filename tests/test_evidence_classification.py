from GeneFior.evidence import (
    ALLELE_LIKE,
    CANDIDATE_ALLELE_DETECTED,
    EXACT_ALLELE_DETECTED,
    EXACT_PROTEIN_DETECTED,
    WHOLE_GENOME_MAPPED,
    FAMILY_DETECTED,
    LEGACY_RELAXED_DETECTED,
    MIXED_OR_MOSAIC,
    EvidenceConfig,
    classify_gene_evidence,
    classify_gene_legacy,
    infer_gene_family,
    normalise_detection_system,
)
from GeneFior.gene_stats import GeneStats


def _classify(stats, tool="BLASTn", **overrides):
    config = EvidenceConfig(**overrides)
    return classify_gene_evidence(
        gene="blaCTX-M-1_1_DQ915955",
        tool_name=tool,
        stats=stats,
        config=config,
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_num_reads=1,
        reads_mode=True,
    )


def test_family_name_collapses_common_resfinder_alleles():
    assert infer_gene_family("blaCTX-M-1_1_DQ915955") == "blaCTX-M"
    assert infer_gene_family("blaCTX-M-36_1_AB177384") == "blaCTX-M"
    assert infer_gene_family("tet(A)_1_AJ313332") == "tet(A)"
    assert infer_gene_family("SHV-52|ARO:3001109") == "SHV"
    assert infer_gene_family("aadA1_1_X12870") == "aadA"


def test_single_read_bridge_is_not_an_exact_detection():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    stats.add_hit(1, 77, 100.0, 292)
    stats.add_hit(90, 139, 80.0, 292)
    stats.add_hit(144, 292, 100.0, 292)
    stats.add_read_support(
        mapped=True,
        passing=True,
        best=True,
        ambiguous_best=True,
    )
    stats.finalise()

    call = _classify(stats, tool="DIAMOND-BLASTX")

    assert call.status == MIXED_OR_MOSAIC
    assert not call.evidence_present
    assert not call.exact_detected
    assert "LOW_CORROBORATED_COVERAGE" in call.warnings
    assert "PROTEIN_ONLY_ALLELE_AMBIGUITY" in call.warnings


def test_legacy_relaxed_mode_does_not_accept_single_read_bridge():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    stats.add_hit(1, 77, 100.0, 292)
    stats.add_hit(90, 139, 80.0, 292)
    stats.add_hit(144, 292, 100.0, 292)
    stats.add_read_support(
        mapped=True,
        passing=True,
        best=True,
        ambiguous_best=True,
    )
    stats.finalise()

    qualified = _classify(stats, tool="DIAMOND-BLASTX")
    legacy = classify_gene_legacy(
        gene="blaCTX-M-1_1_DQ915955",
        tool_name="DIAMOND-BLASTX",
        stats=stats,
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_depth=1,
        detection_min_num_reads=1,
    )

    assert qualified.status == MIXED_OR_MOSAIC
    assert legacy.status != LEGACY_RELAXED_DETECTED
    assert not legacy.evidence_present
    assert not legacy.exact_allele_detected
    assert legacy.warnings == ["LEGACY_RELAXED_RULES"]


def test_legacy_depth_gate_enforces_three_x_floor():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    stats.add_hit(1, 876, 99.0, 876)
    stats.finalise()
    three_read_stats = GeneStats("blaCTX-M-1_1_DQ915955")
    for _ in range(3):
        three_read_stats.add_hit(1, 876, 99.0, 876)
    three_read_stats.finalise()

    requested_one_x = classify_gene_legacy(
        gene="blaCTX-M-1_1_DQ915955",
        tool_name="BLASTn",
        stats=stats,
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_depth=1,
        detection_min_num_reads=1,
    )
    actual_three_x = classify_gene_legacy(
        gene="blaCTX-M-1_1_DQ915955",
        tool_name="BLASTn",
        stats=three_read_stats,
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_depth=1,
        detection_min_num_reads=1,
    )
    requested_three_x = classify_gene_legacy(
        gene="blaCTX-M-1_1_DQ915955",
        tool_name="BLASTn",
        stats=stats,
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_depth=3,
        detection_min_num_reads=1,
    )

    assert not requested_one_x.evidence_present
    assert actual_three_x.evidence_present
    assert not requested_three_x.evidence_present


def test_detection_system_aliases_are_normalised():
    assert normalise_detection_system("new") == "qualified"
    assert normalise_detection_system("legacy") == "legacy-relaxed"


def test_legacy_mode_preserves_historical_hsp_counting_for_min_reads():
    stats = GeneStats("geneA")
    stats.add_hit(1, 100, 99.0, 100)
    stats.add_hit(1, 100, 99.0, 100)
    stats.add_hit(1, 100, 99.0, 100)
    stats.add_read_support(
        mapped=True,
        passing=True,
        best=True,
        unique_best=True,
    )
    stats.finalise()

    qualified = classify_gene_evidence(
        gene="geneA",
        tool_name="BLASTn",
        stats=stats,
        config=EvidenceConfig(),
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_num_reads=3,
        reads_mode=True,
    )
    legacy = classify_gene_legacy(
        gene="geneA",
        tool_name="BLASTn",
        stats=stats,
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_depth=1,
        detection_min_num_reads=3,
    )

    assert not qualified.evidence_present
    assert legacy.evidence_present


def test_full_protein_coverage_is_exact_at_protein_level_not_nucleotide_level():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    for _ in range(5):
        stats.add_hit(1, 292, 100.0, 292)
        stats.add_read_support(
            mapped=True,
            passing=True,
            best=True,
            unique_best=True,
        )
    stats.finalise()

    call = _classify(stats, tool="DIAMOND-BLASTX")

    assert call.status == EXACT_PROTEIN_DETECTED
    assert call.evidence_present
    assert not call.exact_detected
    assert call.exact_protein_detected
    assert not call.exact_allele_detected


def test_exact_nucleotide_call_requires_correlated_unique_support():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    for _ in range(3):
        stats.add_hit(1, 876, 100.0, 876)
        stats.add_read_support(
            mapped=True,
            passing=True,
            best=True,
            unique_best=True,
            high_confidence=True,
        )
    stats.finalise()

    call = _classify(stats, tool="BWA")

    assert call.status == EXACT_ALLELE_DETECTED
    assert call.evidence_present
    assert call.exact_detected


def test_near_exact_identity_is_not_labelled_literal_exact():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    for _ in range(3):
        stats.add_hit(1, 876, 99.5, 876)
        stats.add_read_support(
            mapped=True,
            passing=True,
            best=True,
            unique_best=True,
            high_confidence=True,
        )
    stats.finalise()

    call = _classify(stats, tool="BWA")

    assert call.status == CANDIDATE_ALLELE_DETECTED
    assert call.candidate_allele_detected
    assert not call.exact_detected
    assert "NON_IDENTICAL_ALLELE_SEQUENCE" in call.warnings


def test_candidate_allele_requires_full_three_x_depth():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    for _ in range(2):
        stats.add_hit(1, 876, 99.5, 876)
        stats.add_read_support(
            mapped=True,
            passing=True,
            best=True,
            unique_best=True,
            high_confidence=True,
        )
    stats.finalise()

    call = _classify(stats, tool="BWA")

    assert call.status == MIXED_OR_MOSAIC
    assert not call.candidate_allele_detected
    assert "INCOMPLETE_CANDIDATE_ALLELE_3X_COVERAGE" in call.warnings
    assert not call.evidence_present


def test_candidate_allele_requires_minimum_identity():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    for _ in range(3):
        stats.add_hit(1, 876, 97.9, 876)
        stats.add_read_support(
            mapped=True,
            passing=True,
            best=True,
            unique_best=True,
            high_confidence=True,
        )
    stats.finalise()

    call = _classify(stats, tool="BWA")

    assert call.status == ALLELE_LIKE
    assert not call.candidate_allele_detected
    assert "LOW_CANDIDATE_ALLELE_IDENTITY" in call.warnings


def test_ambiguous_nucleotide_support_is_family_level():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    for _ in range(4):
        stats.add_hit(1, 876, 99.5, 876)
        stats.add_read_support(
            mapped=True,
            passing=True,
            best=True,
            ambiguous_best=True,
        )
    stats.finalise()

    call = _classify(stats, tool="BLASTn")

    assert call.status == FAMILY_DETECTED
    assert call.evidence_present
    assert not call.exact_detected
    assert "MULTI_ALLELE_READS" in call.warnings


def test_full_length_dna_gene_can_be_candidate_at_one_x_coverage():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    stats.add_hit(1, 876, 99.5, 876)
    stats.add_read_support(mapped=True, passing=True, best=True)
    stats.finalise()

    call = classify_gene_evidence(
        gene="blaCTX-M-1_1_DQ915955",
        tool_name="BLASTn",
        stats=stats,
        config=EvidenceConfig(),
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_num_reads=1,
        reads_mode=False,
    )

    assert call.status == CANDIDATE_ALLELE_DETECTED
    assert call.candidate_allele_detected
    assert not call.exact_allele_detected


def test_full_length_dna_gene_can_be_exact_at_one_x_coverage():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    stats.add_hit(1, 876, 100.0, 876)
    stats.add_read_support(mapped=True, passing=True, best=True)
    stats.finalise()

    call = classify_gene_evidence(
        gene="blaCTX-M-1_1_DQ915955",
        tool_name="BLASTn",
        stats=stats,
        config=EvidenceConfig(),
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_num_reads=1,
        reads_mode=False,
    )

    assert call.status == EXACT_ALLELE_DETECTED
    assert call.exact_allele_detected


def test_full_length_protein_gene_is_exact_protein_not_nucleotide_allele():
    stats = GeneStats("blaCTX-M-1_1_DQ915955")
    stats.add_hit(1, 292, 100.0, 292)
    stats.add_read_support(mapped=True, passing=True, best=True)
    stats.finalise()

    call = classify_gene_evidence(
        gene="blaCTX-M-1_1_DQ915955",
        tool_name="DIAMOND-BLASTP",
        stats=stats,
        config=EvidenceConfig(),
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_num_reads=1,
        reads_mode=False,
    )

    assert call.status == EXACT_PROTEIN_DETECTED
    assert call.evidence_present
    assert not call.candidate_allele_detected
    assert not call.exact_allele_detected
    assert call.exact_protein_detected
    assert "PROTEIN_ONLY_ALLELE_AMBIGUITY" in call.warnings


def test_whole_genome_mapping_does_not_claim_gene_or_allele_resolution():
    stats = GeneStats("chromosome_1")
    for _ in range(3):
        stats.add_hit(1, 1000, 99.5, 1000)
        stats.add_read_support(mapped=True, passing=True, best=True)
    stats.finalise()

    call = classify_gene_evidence(
        gene="chromosome_1",
        tool_name="BWA",
        stats=stats,
        config=EvidenceConfig(),
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_num_reads=1,
        reads_mode=True,
        whole_genome=True,
    )

    assert call.status == WHOLE_GENOME_MAPPED
    assert call.evidence_present
    assert not call.exact_allele_detected
    assert not call.candidate_allele_detected


def test_partial_high_identity_mapping_is_not_positive_evidence():
    stats = GeneStats("aadE-Cc_1_CP013733")
    stats.add_hit(1, 50, 100.0, 300)
    stats.add_hit(82, 130, 100.0, 300)
    stats.add_read_support(
        mapped=True,
        passing=True,
        best=True,
        unique_best=True,
    )
    stats.add_read_support(
        mapped=True,
        passing=True,
        best=True,
        unique_best=True,
    )
    stats.finalise()

    call = classify_gene_evidence(
        gene="aadE-Cc_1_CP013733",
        tool_name="DIAMOND-BLASTX",
        stats=stats,
        config=EvidenceConfig(),
        detection_min_coverage=80.0,
        detection_min_identity=80.0,
        detection_min_num_reads=1,
        reads_mode=True,
    )

    assert call.status == "PARTIAL_OR_DIVERGENT"
    assert not call.evidence_present
    assert not call.exact_detected
    assert "PARTIAL_COVERAGE" in call.warnings
