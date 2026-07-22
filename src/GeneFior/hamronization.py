"""Export finalised GeneFior calls in the PHA4GE hAMRonization schema."""

from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

try:
    from .evidence import infer_gene_family
except (ImportError, ModuleNotFoundError):
    from evidence import infer_gene_family


HAMRONIZED_FIELDS = [
    "input_file_name",
    "gene_symbol",
    "gene_name",
    "reference_database_name",
    "reference_database_version",
    "reference_accession",
    "analysis_software_name",
    "analysis_software_version",
    "genetic_variation_type",
    "antimicrobial_agent",
    "coverage_percentage",
    "coverage_depth",
    "coverage_ratio",
    "drug_class",
    "input_gene_length",
    "input_gene_start",
    "input_gene_stop",
    "input_protein_length",
    "input_protein_start",
    "input_protein_stop",
    "input_sequence_id",
    "nucleotide_mutation",
    "nucleotide_mutation_interpretation",
    "predicted_phenotype",
    "predicted_phenotype_confidence_level",
    "amino_acid_mutation",
    "amino_acid_mutation_interpretation",
    "reference_gene_length",
    "reference_gene_start",
    "reference_gene_stop",
    "reference_protein_length",
    "reference_protein_start",
    "reference_protein_stop",
    "resistance_mechanism",
    "strand_orientation",
    "sequence_identity",
]

_CALL_RANK = {
    "EXACT_ALLELE_DETECTED": 7,
    "CANDIDATE_ALLELE_DETECTED": 6,
    "ALLELE_LIKE": 5,
    "PROFILE_DETECTED": 4,
    "FAMILY_DETECTED": 3,
    "MIXED_OR_MOSAIC": 2,
    "DETECTED": 1,
}


def add_hamronization_arguments(parser, inherit_defaults: bool = False):
    """Add the shared hAMRonization output options to an argparse parser."""
    minimum_help = (
        "source run, otherwise evidence"
        if inherit_defaults else "evidence"
    )
    group = parser.add_argument_group("hAMRonization Output")
    group.add_argument(
        "--hamronized-output",
        choices=["tsv", "json", "both"],
        default=None,
        help="Write positive calls using the PHA4GE hAMRonization schema.",
    )
    group.add_argument(
        "--hamronized-min-call",
        choices=["evidence", "candidate", "exact"],
        default=None if inherit_defaults else "evidence",
        help=(
            "Minimum call tier included in the harmonized report "
            f"(default: {minimum_help})."
        ),
    )
    group.add_argument(
        "--sample-id",
        default=None,
        help=(
            "Sample identifier used as input_file_name. Single runs infer it "
            "from the input when omitted; batch runs use each sample name."
        ),
    )
    group.add_argument(
        "--database-version",
        dest="database_version_specs",
        action="append",
        default=None,
        metavar="DATABASE=VERSION",
        help=(
            "Reference database version required by hAMRonization. Repeat "
            "once per selected database."
        ),
    )
    return group


def parse_database_versions(values: Optional[Iterable[str]]) -> Dict[str, str]:
    """Parse repeated DATABASE=VERSION values into a normalized mapping."""
    versions: Dict[str, str] = {}
    for raw_value in values or []:
        value = str(raw_value).strip()
        if "=" not in value:
            raise ValueError(
                f"Invalid database version '{raw_value}'; expected DATABASE=VERSION"
            )
        database, version = (part.strip() for part in value.split("=", 1))
        if not database or not version:
            raise ValueError(
                f"Invalid database version '{raw_value}'; database and version "
                "must both be non-empty"
            )
        key = database.casefold()
        if key in versions and versions[key] != version:
            raise ValueError(
                f"Conflicting versions supplied for database '{database}'"
            )
        versions[key] = version
    return versions


def validate_database_versions(
        output_format: Optional[str], database_names: Iterable[str],
        database_versions: Mapping[str, str]) -> None:
    """Require explicit versions before starting a harmonized analysis."""
    if not output_format:
        return
    normalized = {str(key).casefold(): value for key, value in database_versions.items()}
    missing = sorted(
        str(database) for database in database_names
        if not normalized.get(str(database).casefold())
    )
    if missing:
        options = " ".join(
            f"--database-version {database}=VERSION" for database in missing
        )
        raise ValueError(
            "hAMRonization requires a reference database version for every "
            f"selected database. Missing: {', '.join(missing)}. Supply {options}"
        )


def validate_hamronization_request(
        output_format: Optional[str], minimum_call: Optional[str],
        detection_system: str) -> None:
    """Reject allele-tier exports from the binary legacy detector."""
    if (
            output_format
            and minimum_call in {"candidate", "exact"}
            and str(detection_system).casefold() == "legacy-relaxed"):
        raise ValueError(
            f"--hamronized-min-call {minimum_call} requires "
            "--detection-system qualified; legacy-relaxed can export only "
            "the evidence tier"
        )


def infer_sample_id(input_spec=None, output_dir=None) -> str:
    """Infer a stable sample identifier from single or paired sequence inputs."""
    if isinstance(input_spec, (tuple, list)):
        raw_paths = [str(value) for value in input_spec if value]
    else:
        raw_paths = [
            part.strip() for part in str(input_spec or "").split(",")
            if part.strip()
        ]

    names = [_strip_sequence_suffixes(Path(path).name) for path in raw_paths]
    if len(names) > 1:
        prefix = os.path.commonprefix(names).rstrip("._- ")
        if prefix:
            return prefix
    if names:
        name = re.sub(
            r"(?:[._-](?:R?[12]))(?:[._-](?:trimmed|cleaned|filtered))?$",
            "",
            names[0],
            flags=re.IGNORECASE,
        ).rstrip("._- ")
        if name:
            return name
    if output_dir:
        return Path(output_dir).resolve().name
    return "sample"


def _strip_sequence_suffixes(name: str) -> str:
    result = name
    for suffix in (".gz", ".gzip", ".bz2", ".xz"):
        if result.lower().endswith(suffix):
            result = result[:-len(suffix)]
            break
    for suffix in (".fastq", ".fq", ".fasta", ".fna", ".faa", ".fsa", ".fa"):
        if result.lower().endswith(suffix):
            result = result[:-len(suffix)]
            break
    return result


def _as_bool(value) -> bool:
    return str(value or "").strip().casefold() in {"1", "true", "yes", "y"}


def _as_float(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _as_int(value) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _sequence_basis(tool: str) -> str:
    name = str(tool).upper()
    if "HMMER-DNA" in name or "NHMMER" in name:
        return "nucleotide"
    if any(token in name for token in ("BLASTX", "BLASTP", "DIAMOND")):
        return "protein"
    if "HMMER" in name:
        return "protein"
    return "nucleotide"


def _display_gene_name(gene: str) -> str:
    value = str(gene or "").strip()
    if "|" in value:
        parts = [part.strip() for part in value.split("|") if part.strip()]
        non_accessions = [
            part for part in parts
            if not re.fullmatch(r"ARO:\d+", part, flags=re.IGNORECASE)
        ]
        if non_accessions:
            value = non_accessions[0]
    value = re.sub(
        r"_[0-9]+_[A-Z]{1,6}_?[A-Z0-9]*(?:\.[0-9]+)?$",
        "",
        value,
    )
    return value or str(gene)


def _reference_accession(gene: str) -> str:
    value = str(gene or "").strip()
    aro = re.search(r"ARO:\d+", value, flags=re.IGNORECASE)
    if aro:
        return aro.group(0).upper()
    accessions = re.findall(
        r"(?<![A-Za-z0-9])[A-Z]{1,6}_?\d{3,}(?:\.\d+)?(?![A-Za-z0-9])",
        value,
    )
    if accessions:
        return accessions[-1]
    suffix = re.search(
        r"_([A-Z]{1,6}_?[A-Z0-9]*\d[A-Z0-9]*(?:\.\d+)?)$",
        value,
    )
    if suffix:
        return suffix.group(1)
    return value


def _normalized_family(gene: str, reported_family: str) -> str:
    family = str(reported_family or "").strip()
    if not family or re.fullmatch(r"ARO:\d+", family, re.IGNORECASE):
        family = gene
    aliases = [part.strip() for part in family.split("/") if part.strip()]
    normalized_aliases = [infer_gene_family(alias) for alias in aliases]
    if normalized_aliases and len({
            alias.casefold() for alias in normalized_aliases
    }) == 1:
        return normalized_aliases[0]
    return infer_gene_family(family)


def _stats_file_identity(path: Path) -> Tuple[str, str]:
    name = path.name
    if not name.endswith("_stats.tsv"):
        raise ValueError(f"Not a GeneFior tool-stat path: {path}")
    stem = name[:-len("_stats.tsv")]
    parts = stem.rsplit("_", 1)
    if len(parts) != 2 or not all(parts):
        raise ValueError(f"Cannot identify database and tool from {path}")
    return parts[0], parts[1]


def _row_is_selected(row: Mapping[str, str], minimum_call: str) -> bool:
    qualified = "Evidence_Present" in row
    if qualified:
        if minimum_call == "evidence":
            return _as_bool(row.get("Evidence_Present"))
        if minimum_call == "candidate":
            return (
                _as_bool(row.get("Candidate_Allele_Detected"))
                or _as_bool(row.get("Exact_Allele_Detected"))
            )
        return _as_bool(row.get("Exact_Allele_Detected"))
    if minimum_call != "evidence":
        raise ValueError(
            f"The '{minimum_call}' hAMRonization tier is unavailable for "
            "legacy-relaxed results; use --hamronized-min-call evidence or "
            "recompute with --detection-system qualified"
        )
    return _as_bool(row.get("Detected"))


def _candidate_rank(row: Mapping[str, str], tool: str) -> Tuple:
    status = row.get("Evidence_Status") or (
        "DETECTED" if _as_bool(row.get("Detected")) else "NOT_DETECTED"
    )
    basis_rank = 2 if _sequence_basis(tool) == "nucleotide" else 1
    return (
        _CALL_RANK.get(status, 0),
        basis_rank,
        _as_float(row.get("Detection_Depth_Coverage")),
        _as_float(row.get("Gene_Coverage")),
        _as_float(row.get("Avg_Identity")),
        _as_float(row.get("Base_Coverage")),
        _as_int(row.get("Unique_Best_Read_Support")),
    )


def _collect_candidates(stats_dir: Path, minimum_call: str) -> List[dict]:
    selected: Dict[Tuple[str, str], dict] = {}
    for stats_path in sorted(stats_dir.glob("*_stats.tsv")):
        database, tool = _stats_file_identity(stats_path)
        with open(stats_path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if not row.get("Gene") or not _row_is_selected(row, minimum_call):
                    continue
                gene = row["Gene"].strip()
                family = _normalized_family(gene, row.get("Family"))
                record = {
                    "database": database,
                    "tool": tool,
                    "gene": gene,
                    "family": family,
                    "row": row,
                }
                key = (database.casefold(), family.casefold())
                previous = selected.get(key)
                if previous is None or _candidate_rank(row, tool) > _candidate_rank(
                        previous["row"], previous["tool"]):
                    selected[key] = record
    return sorted(
        selected.values(),
        key=lambda item: (
            item["database"].casefold(), item["family"].casefold(), item["gene"]
        ),
    )


def _hamronized_row(
        candidate: Mapping[str, object], sample_id: str,
        database_versions: Mapping[str, str], analysis_software_name: str,
        analysis_software_version: str) -> dict:
    database = str(candidate["database"])
    tool = str(candidate["tool"])
    source = candidate["row"]
    gene = str(candidate["gene"])
    family = str(candidate["family"])
    status = source.get("Evidence_Status", "DETECTED")
    top_candidate = source.get("Top_Database_Candidate") or gene
    allele_resolved = status in {
        "EXACT_ALLELE_DETECTED", "CANDIDATE_ALLELE_DETECTED"
    }
    symbol = _display_gene_name(top_candidate) if allele_resolved else _display_gene_name(family)
    basis = _sequence_basis(tool)
    gene_length = _as_int(source.get("Gene_Length"))
    version = database_versions[database.casefold()]

    result = {field: None for field in HAMRONIZED_FIELDS}
    result.update({
        "input_file_name": sample_id,
        "gene_symbol": symbol,
        "gene_name": symbol,
        "reference_database_name": database,
        "reference_database_version": version,
        "reference_accession": _reference_accession(top_candidate),
        "analysis_software_name": analysis_software_name,
        "analysis_software_version": analysis_software_version,
        "genetic_variation_type": "gene_presence_detected",
        "coverage_percentage": _as_float(source.get("Gene_Coverage")),
        "coverage_depth": _as_float(source.get("Base_Coverage")),
        "sequence_identity": _as_float(source.get("Avg_Identity")),
    })
    if basis == "nucleotide":
        result["reference_gene_length"] = gene_length
    else:
        result["reference_protein_length"] = gene_length
    return result


def export_hamronized_reports(
        output_dir, output_format: str, minimum_call: str,
        database_versions: Mapping[str, str], analysis_software_name: str,
        analysis_software_version: str, sample_id: Optional[str] = None,
        input_spec=None, logger=None) -> int:
    """Write consolidated hAMRonization TSV/JSON reports from tool statistics."""
    root = Path(output_dir)
    stats_dir = root / "tool_stats"
    if not stats_dir.is_dir():
        raise ValueError(f"GeneFior tool-stat directory was not found: {stats_dir}")
    normalized_versions = {
        str(key).casefold(): str(value)
        for key, value in database_versions.items()
        if str(value).strip()
    }
    candidates = _collect_candidates(stats_dir, minimum_call)
    missing = sorted({
        item["database"] for item in candidates
        if item["database"].casefold() not in normalized_versions
    })
    if missing:
        raise ValueError(
            "Cannot write hAMRonization output without database version(s): "
            + ", ".join(missing)
        )

    effective_sample_id = sample_id or infer_sample_id(input_spec, root)
    rows = [
        _hamronized_row(
            candidate,
            effective_sample_id,
            normalized_versions,
            analysis_software_name,
            analysis_software_version,
        )
        for candidate in candidates
    ]

    formats = ("tsv", "json") if output_format == "both" else (output_format,)
    for format_name in formats:
        output_path = root / f"hamronized_report.{format_name}"
        if format_name == "tsv":
            with open(output_path, "w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    delimiter="\t",
                    fieldnames=HAMRONIZED_FIELDS,
                    lineterminator=os.linesep,
                )
                writer.writeheader()
                writer.writerows(rows)
        elif format_name == "json":
            clean_rows = [
                {
                    key: str(value) if value else ""
                    for key, value in row.items()
                }
                for row in rows
            ]
            with open(output_path, "w") as handle:
                json.dump(clean_rows, handle)
                handle.write("\n")
        else:
            raise ValueError(f"Unsupported hAMRonization output format: {format_name}")
        if logger:
            logger.info(
                f"Wrote {len(rows)} hAMRonized result(s) to: {output_path}"
            )
    return len(rows)


def export_from_options(
        options, output_dir, analysis_software_name: str,
        analysis_software_version: str, logger=None,
        sample_id: Optional[str] = None, input_spec=None) -> int:
    """Export when requested using the normalized attributes on CLI options."""
    output_format = getattr(options, "hamronized_output", None)
    if not output_format:
        return 0
    return export_hamronized_reports(
        output_dir=output_dir,
        output_format=output_format,
        minimum_call=getattr(options, "hamronized_min_call", None) or "evidence",
        database_versions=getattr(options, "database_versions", {}) or {},
        analysis_software_name=analysis_software_name,
        analysis_software_version=analysis_software_version,
        sample_id=sample_id or getattr(options, "sample_id", None),
        input_spec=input_spec,
        logger=logger,
    )
