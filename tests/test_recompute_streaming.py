import json
import logging
from types import SimpleNamespace

from GeneFior.GeneFior_Recompute import (
    apply_source_run_defaults,
    detect_diamond_mode,
    parse_query_paths,
    validate_recompute_relaxation,
)


def _options(input_dir, **overrides):
    values = {
        "input": str(input_dir),
        "tools": ["blastx"],
        "query_min_coverage": 40.0,
        "query_min_identity": 80.0,
        "min_bitscore": None,
        "evalue": None,
        "allow_incomplete_relaxation": False,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_recompute_allows_equal_or_stricter_source_thresholds(tmp_path):
    (tmp_path / "run_parameters.json").write_text(json.dumps({
        "query_min_coverage": 40.0,
        "query_min_identity": 80.0,
        "min_bitscore": None,
        "evalue": None,
    }))
    options = _options(
        tmp_path,
        query_min_coverage=60.0,
        query_min_identity=90.0,
    )
    assert validate_recompute_relaxation(options, logging.getLogger("test"))


def test_recompute_rejects_incomplete_query_relaxation(tmp_path):
    (tmp_path / "run_parameters.json").write_text(json.dumps({
        "query_min_coverage": 40.0,
        "query_min_identity": 80.0,
        "min_bitscore": None,
        "evalue": None,
    }))
    options = _options(tmp_path, query_min_identity=70.0)
    assert not validate_recompute_relaxation(options, logging.getLogger("test"))


def test_recompute_can_explicitly_allow_incomplete_relaxation(tmp_path):
    (tmp_path / "run_parameters.json").write_text(json.dumps({
        "query_min_coverage": 40.0,
        "query_min_identity": 80.0,
        "min_bitscore": None,
        "evalue": None,
    }))
    options = _options(
        tmp_path,
        query_min_identity=70.0,
        allow_incomplete_relaxation=True,
    )
    assert validate_recompute_relaxation(options, logging.getLogger("test"))


def test_bam_only_recompute_is_not_blocked_by_blast_source_thresholds(tmp_path):
    (tmp_path / "run_parameters.json").write_text(json.dumps({
        "query_min_coverage": 40.0,
        "query_min_identity": 80.0,
    }))
    options = _options(
        tmp_path,
        tools=["bwa"],
        query_min_identity=20.0,
    )
    assert validate_recompute_relaxation(options, logging.getLogger("test"))


def test_recompute_inherits_source_sequence_and_evidence_metadata(tmp_path):
    (tmp_path / "run_parameters.json").write_text(json.dumps({
        "sequence_type": "Genes-FASTA",
        "genes_type": "protein",
        "query_min_coverage": 55.0,
        "query_min_identity": 91.0,
        "detection_min_coverage": 87.0,
        "evidence_corroborating_depth": 3,
        "hamronized_output": "both",
        "hamronized_min_call": "candidate",
        "sample_id": "source-sample",
        "database_versions": {"resfinder": "2026-01"},
    }))
    options = _options(
        tmp_path,
        tools=["all"],
        sequence_type=None,
        genes_type=None,
        query_min_coverage=None,
        query_min_identity=None,
        detection_min_coverage=None,
        detection_min_identity=None,
        detection_min_depth=None,
        detection_min_num_reads=None,
        evidence_corroborating_depth=None,
        evidence_exact_identity=None,
        evidence_max_internal_gap_bp=None,
        evidence_max_internal_gap_fraction=None,
        evidence_min_unique_reads=None,
        evidence_min_unique_fraction=None,
        evidence_ambiguity_fraction=None,
        evidence_score_tie=None,
        hamronized_output=None,
        hamronized_min_call=None,
        sample_id=None,
    )

    apply_source_run_defaults(options, logging.getLogger("test"))

    assert options.sequence_type == "Genes-FASTA"
    assert options.genes_type == "aa"
    assert options.query_min_coverage == 55.0
    assert options.query_min_identity == 91.0
    assert options.detection_min_coverage == 87.0
    assert options.detection_min_depth == 1
    assert options.evidence_corroborating_depth == 3
    assert options.hamronized_output == "both"
    assert options.hamronized_min_call == "candidate"
    assert options.sample_id == "source-sample"
    assert options.source_run_parameters["database_versions"] == {
        "resfinder": "2026-01"
    }
    assert options.tools == ["all"]


def test_recompute_maps_old_base_depth_manifest_to_detection_depth(tmp_path):
    (tmp_path / "run_parameters.json").write_text(json.dumps({
        "detection_min_base_depth": 4,
    }))
    options = _options(
        tmp_path,
        detection_min_depth=None,
        detection_min_coverage=None,
        detection_min_identity=None,
        detection_min_num_reads=None,
        query_min_coverage=None,
        query_min_identity=None,
        evidence_corroborating_depth=None,
        evidence_exact_identity=None,
        evidence_candidate_depth=None,
        evidence_candidate_identity=None,
        evidence_max_internal_gap_bp=None,
        evidence_max_internal_gap_fraction=None,
        evidence_min_unique_reads=None,
        evidence_min_unique_fraction=None,
        evidence_ambiguity_fraction=None,
        evidence_score_tie=None,
        sequence_type=None,
        genes_type=None,
    )

    apply_source_run_defaults(options, logging.getLogger("test"))

    assert options.detection_min_depth == 4


def test_diamond_mode_prefers_source_manifest_for_protein_genes(tmp_path):
    diamond_file = tmp_path / "db_diamond_results.tsv"
    diamond_file.write_text("")

    assert detect_diamond_mode(
        diamond_file,
        {"sequence_type": "Genes-FASTA", "genes_type": "protein"},
    ) == "DIAMOND-BLASTP"


def test_parse_query_paths_supports_paired_files(tmp_path):
    read1 = tmp_path / "reads_1.fastq.gz"
    read2 = tmp_path / "reads_2.fastq.gz"

    assert parse_query_paths(f"{read1},{read2}") == [read1, read2]
