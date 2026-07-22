import csv
import json
from pathlib import Path

import pytest

from GeneFior.evidence import infer_gene_family
from GeneFior.hamronization import (
    HAMRONIZED_FIELDS,
    export_hamronized_reports,
    infer_sample_id,
    parse_database_versions,
    validate_database_versions,
    validate_hamronization_request,
)


QUALIFIED_FIELDS = [
    "Gene",
    "Gene_Length",
    "Gene_Coverage",
    "Base_Coverage",
    "Detection_Depth_Coverage",
    "Avg_Identity",
    "Unique_Best_Read_Support",
    "Evidence_Status",
    "Family",
    "Top_Database_Candidate",
    "Evidence_Present",
    "Candidate_Allele_Detected",
    "Exact_Allele_Detected",
    "Detected",
]


def _write_stats(root: Path, database: str, tool: str, rows, fields=None):
    stats_dir = root / "tool_stats"
    stats_dir.mkdir(parents=True, exist_ok=True)
    path = stats_dir / f"{database}_{tool}_stats.tsv"
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fields or QUALIFIED_FIELDS,
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def _qualified_row(gene, status, family, **overrides):
    row = {
        "Gene": gene,
        "Gene_Length": "876",
        "Gene_Coverage": "100.00",
        "Base_Coverage": "42.50",
        "Detection_Depth_Coverage": "100.00",
        "Avg_Identity": "99.20",
        "Unique_Best_Read_Support": "12",
        "Evidence_Status": status,
        "Family": family,
        "Top_Database_Candidate": gene,
        "Evidence_Present": "1",
        "Candidate_Allele_Detected": "0",
        "Exact_Allele_Detected": "0",
        "Detected": "1",
    }
    row.update({key: str(value) for key, value in overrides.items()})
    return row


def test_database_version_parser_and_validation():
    assert parse_database_versions(["ResFinder=2026-01", "card=4.0.1"]) == {
        "resfinder": "2026-01",
        "card": "4.0.1",
    }
    with pytest.raises(ValueError, match="DATABASE=VERSION"):
        parse_database_versions(["resfinder"])
    with pytest.raises(ValueError, match="Conflicting versions"):
        parse_database_versions(["card=4.0", "CARD=4.1"])
    with pytest.raises(ValueError, match="Missing: card"):
        validate_database_versions("tsv", ["resfinder", "card"], {
            "resfinder": "2026-01"
        })
    validate_database_versions(None, ["resfinder"], {})
    with pytest.raises(ValueError, match="requires --detection-system qualified"):
        validate_hamronization_request("tsv", "candidate", "legacy-relaxed")
    validate_hamronization_request("tsv", "evidence", "legacy-relaxed")


def test_sample_id_inference_for_paired_and_single_inputs(tmp_path):
    paired = (
        "/data/232758E_18_1_trimmed.fastq.gz,"
        "/data/232758E_18_2_trimmed.fastq.gz"
    )
    assert infer_sample_id(paired, tmp_path) == "232758E_18"
    assert infer_sample_id("/data/isolate_A.fasta.gz", tmp_path) == "isolate_A"
    assert infer_sample_id(None, tmp_path) == tmp_path.name


def test_qualified_export_consolidates_tools_and_excludes_review_rows(tmp_path):
    _write_stats(tmp_path, "resfinder", "BLASTx", [
        _qualified_row(
            "blaCTX-M-1_1_DQ915955",
            "FAMILY_DETECTED",
            "blaCTX-M",
            Gene_Length=292,
            Avg_Identity=100,
        ),
        _qualified_row(
            "review_1_X00001",
            "MUST_FLAG_REVIEW",
            "review",
            Evidence_Present=0,
            Detected=0,
        ),
    ])
    _write_stats(tmp_path, "resfinder", "BLASTn", [
        _qualified_row(
            "blaCTX-M-1_1_DQ915955",
            "EXACT_ALLELE_DETECTED",
            "blaCTX-M",
            Avg_Identity=100,
            Candidate_Allele_Detected=1,
            Exact_Allele_Detected=1,
        ),
    ])

    count = export_hamronized_reports(
        tmp_path,
        output_format="both",
        minimum_call="evidence",
        database_versions={"resfinder": "2026-01"},
        analysis_software_name="GeneFior",
        analysis_software_version="0.10.3",
        sample_id="sample-1",
    )

    assert count == 1
    with open(tmp_path / "hamronized_report.tsv", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == HAMRONIZED_FIELDS
    assert rows[0]["input_file_name"] == "sample-1"
    assert rows[0]["gene_symbol"] == "blaCTX-M-1"
    assert rows[0]["reference_accession"] == "DQ915955"
    assert rows[0]["reference_database_name"] == "resfinder"
    assert rows[0]["reference_database_version"] == "2026-01"
    assert rows[0]["reference_gene_length"] == "876"
    assert rows[0]["reference_protein_length"] == ""
    assert rows[0]["genetic_variation_type"] == "gene_presence_detected"

    json_rows = json.loads((tmp_path / "hamronized_report.json").read_text())
    assert len(json_rows) == 1
    assert json_rows[0]["coverage_percentage"] == "100.0"
    assert json_rows[0]["analysis_software_version"] == "0.10.3"


def test_protein_family_export_uses_family_and_protein_length(tmp_path):
    _write_stats(tmp_path, "resfinder", "DIAMOND-BLASTX", [
        _qualified_row(
            "blaCTX-M-1_1_DQ915955",
            "FAMILY_DETECTED",
            "blaCTX-M",
            Gene_Length=292,
            Gene_Coverage=93.84,
            Avg_Identity=80,
        )
    ])
    export_hamronized_reports(
        tmp_path,
        output_format="tsv",
        minimum_call="evidence",
        database_versions={"resfinder": "2026-01"},
        analysis_software_name="AMRfior",
        analysis_software_version="0.10.3",
        sample_id="sample-2",
    )
    with open(tmp_path / "hamronized_report.tsv", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["gene_symbol"] == "blaCTX-M"
    assert row["reference_accession"] == "DQ915955"
    assert row["reference_gene_length"] == ""
    assert row["reference_protein_length"] == "292"
    assert row["sequence_identity"] == "80.0"


def test_candidate_and_exact_minimum_call_tiers(tmp_path):
    _write_stats(tmp_path, "resfinder", "BLASTn", [
        _qualified_row("family_1_AA00001", "FAMILY_DETECTED", "family"),
        _qualified_row(
            "candidate_1_AA00002",
            "CANDIDATE_ALLELE_DETECTED",
            "candidate",
            Candidate_Allele_Detected=1,
        ),
        _qualified_row(
            "exact_1_AA00003",
            "EXACT_ALLELE_DETECTED",
            "exact",
            Candidate_Allele_Detected=1,
            Exact_Allele_Detected=1,
        ),
    ])
    common = {
        "output_format": "tsv",
        "database_versions": {"resfinder": "2026-01"},
        "analysis_software_name": "GeneFior",
        "analysis_software_version": "0.10.3",
        "sample_id": "sample-3",
    }
    assert export_hamronized_reports(
        tmp_path, minimum_call="candidate", **common
    ) == 2
    assert export_hamronized_reports(
        tmp_path, minimum_call="exact", **common
    ) == 1


def test_legacy_export_is_binary_and_rejects_allele_tiers(tmp_path):
    fields = [
        "Gene", "Gene_Length", "Gene_Coverage", "Base_Coverage",
        "Avg_Identity", "Detected",
    ]
    _write_stats(tmp_path, "resfinder", "BLASTn", [{
        "Gene": "tet(A)_1_AJ313332",
        "Gene_Length": "1200",
        "Gene_Coverage": "95",
        "Base_Coverage": "8",
        "Avg_Identity": "97",
        "Detected": "1",
    }], fields=fields)
    kwargs = {
        "output_format": "tsv",
        "database_versions": {"resfinder": "2026-01"},
        "analysis_software_name": "GeneFior-Recompute",
        "analysis_software_version": "0.10.3",
        "sample_id": "sample-4",
    }
    assert export_hamronized_reports(
        tmp_path, minimum_call="evidence", **kwargs
    ) == 1
    with pytest.raises(ValueError, match="unavailable for legacy-relaxed"):
        export_hamronized_reports(
            tmp_path, minimum_call="candidate", **kwargs
        )


def test_card_family_uses_gene_symbol_not_aro_accession():
    assert infer_gene_family("AAC(6')-Ib7|ARO:3002578") == "AAC(6')-Ib"
    assert infer_gene_family("SHV-52|ARO:3001109") == "SHV"
    assert infer_gene_family("CmlA9_1_JN406319") == "CmlA"


def test_card_alias_names_consolidate_to_one_family(tmp_path):
    _write_stats(tmp_path, "card", "BLASTn", [
        _qualified_row(
            "aadA5|ARO:3002602",
            "FAMILY_DETECTED",
            "ARO:3002602",
        ),
        _qualified_row(
            "aadA6/aadA|ARO:3002603",
            "FAMILY_DETECTED",
            "aadA6/aadA",
        ),
    ])
    assert export_hamronized_reports(
        tmp_path,
        output_format="tsv",
        minimum_call="evidence",
        database_versions={"card": "4.0.1"},
        analysis_software_name="AMRfior",
        analysis_software_version="0.10.3",
        sample_id="sample-5",
    ) == 1
