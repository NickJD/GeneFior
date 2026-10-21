"""DIAMOND translated-read mapping and RPKM reporting for protein gene sets."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


DEFAULT_IDENTITY = 50.0
SPECIAL_IDENTITY = 60.0
DEFAULT_QUERY_COVERAGE = 80.0
SPECIAL_MARKERS = (
    "coxl", "mmox", "bmox", "tmoa", "isoa", "dmpn", "etnc", "prma", "zmoa", "fefe",
)
FIELDS = (
    "qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "qlen", "slen", "qcovhsp", "scovhsp", "stitle",
)


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix.lower() == ".gz" else path.open("rt")


def read_protein_fasta(path: Path) -> dict[str, tuple[int, str]]:
    """Return sequence-id -> (AA length, full title), validating a protein FASTA."""
    records: dict[str, tuple[int, str]] = {}
    title = None
    length = 0
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if title is not None:
                    seqid = title.split()[0]
                    records[seqid] = (length, title)
                title, length = line[1:], 0
            elif title is not None:
                length += len(line)
    if title is not None:
        records[title.split()[0]] = (length, title)
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def count_reads(path: Path) -> int:
    """Count FASTQ or FASTA records without loading sequences into memory."""
    suffixes = [s.lower() for s in path.suffixes]
    is_fastq = any(s in {".fastq", ".fq"} for s in suffixes)
    if not is_fastq and not any(s in {".fasta", ".fa", ".fna", ".fas"} for s in suffixes):
        raise ValueError(f"Cannot infer FASTA/FASTQ format from {path.name}")
    count = 0
    with open_text(path) as handle:
        if is_fastq:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                if not (sequence and plus and quality):
                    raise ValueError(f"Truncated four-line FASTQ record in {path}")
                if not header.startswith("@") or not plus.startswith("+"):
                    raise ValueError(f"Malformed FASTQ record {count + 1} in {path}")
                if len(sequence.strip()) != len(quality.strip()):
                    raise ValueError(f"Sequence/quality length mismatch in {path}, record {count + 1}")
                count += 1
        else:
            count = sum(line.startswith(">") for line in handle)
    return count


def default_identity_rules() -> list[dict]:
    """Method-specific defaults; the same rule format is accepted from a file."""
    rules = [{"pattern": marker, "minimum_identity": SPECIAL_IDENTITY,
              "match_type": "contains", "group_name": marker, "source": "built_in"}
             for marker in SPECIAL_MARKERS]
    rules.append({"pattern": "group4 nife hydrogenase", "minimum_identity": SPECIAL_IDENTITY,
                  "match_type": "all_terms", "group_name": "group 4 [NiFe] hydrogenases",
                  "source": "built_in"})
    return rules


def load_identity_rules(path: Path | None) -> list[dict]:
    """Load generic identity rules from comma- or tab-delimited text."""
    if path is None:
        return []
    with path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        except csv.Error:
            dialect = csv.excel_tab if path.suffix.lower() in {".tsv", ".tab"} else csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError("Identity-threshold file needs a header row")
        names = {name.strip().lower(): name for name in reader.fieldnames if name}
        pattern_col = next((names[n] for n in ("pattern", "gene_id", "gene_group") if n in names), None)
        identity_col = next((names[n] for n in ("minimum_identity", "identity_threshold", "identity_threshold_pct") if n in names), None)
        if not pattern_col or not identity_col:
            raise ValueError("Identity-threshold file needs pattern and minimum_identity columns")
        rules = []
        for line_no, row in enumerate(reader, start=2):
            pattern = (row.get(pattern_col) or "").strip()
            if not pattern:
                continue
            try:
                threshold = float(row.get(identity_col, ""))
            except (TypeError, ValueError):
                raise ValueError(f"Invalid minimum_identity on line {line_no}") from None
            if not 0 <= threshold <= 100:
                raise ValueError(f"minimum_identity on line {line_no} must be in [0,100]")
            match_type = (row.get(names.get("match_type", ""), "exact") or "exact").strip().lower()
            if match_type not in {"exact", "contains", "regex", "all_terms"}:
                raise ValueError(f"Unsupported match_type {match_type!r} on line {line_no}")
            if match_type == "regex":
                try:
                    re.compile(pattern, re.IGNORECASE)
                except re.error as exc:
                    raise ValueError(f"Invalid regex on line {line_no}: {exc}") from None
            group_name = (row.get(names.get("group_name", ""), "") or pattern).strip()
            rules.append({"pattern": pattern, "minimum_identity": threshold,
                          "match_type": match_type, "group_name": group_name,
                          "source": "custom", "line": line_no})
    return rules


def _rule_matches(rule: dict, seqid: str, title: str) -> bool:
    identifier = seqid.strip().lower()
    description = title.strip().lower()
    combined = identifier + " " + description
    pattern = rule["pattern"]
    kind = rule["match_type"]
    if kind == "exact":
        return pattern.lower() == identifier
    if kind == "contains":
        return re.sub(r"[^a-z0-9]+", "", pattern.lower()) in re.sub(r"[^a-z0-9]+", "", combined)
    if kind == "regex":
        return re.search(pattern, combined, re.IGNORECASE) is not None
    compact = re.sub(r"[^a-z0-9]+", "", combined)
    terms = [re.sub(r"[^a-z0-9]+", "", term.lower()) for term in pattern.split()]
    return all(term and term in compact for term in terms)


def identity_threshold(seqid: str, title: str = "", ordinary: float = DEFAULT_IDENTITY,
                       rules: list[dict] | None = None,
                       use_default_exceptions: bool = True) -> tuple[float, str]:
    """Return the first matching configured group threshold, else the baseline."""
    custom = [rule for rule in (rules or []) if rule.get("source") != "built_in"]
    defaults = default_identity_rules() if use_default_exceptions else []
    # Exact IDs outrank broad patterns; then preserve user row order.
    ordered = sorted(custom, key=lambda r: (r["match_type"] != "exact", r.get("line", 0))) + defaults
    for rule in ordered:
        if _rule_matches(rule, seqid, title):
            return float(rule["minimum_identity"]), str(rule["group_name"])
    return ordinary, "default"


def _parse_hit(line: str) -> dict | None:
    fields = line.rstrip("\n").split("\t", len(FIELDS) - 1)
    if len(fields) != len(FIELDS):
        return None
    try:
        row = dict(zip(FIELDS, fields))
        for key in ("pident", "evalue", "bitscore", "qlen", "slen", "qcovhsp", "scovhsp"):
            row[key] = float(row[key])
        row["qstart"], row["qend"] = int(row["qstart"]), int(row["qend"])
        row["sstart"], row["send"] = int(row["sstart"]), int(row["send"])
        return row
    except (ValueError, TypeError):
        return None


def _select_query_hit(query: str, candidates: list[dict]) -> dict | None:
    if not candidates:
        return None
    candidates.sort(key=lambda r: (-r["bitscore"], r["evalue"], -r["pident"], r["sseqid"]))
    top = candidates[0]
    co_best = sorted({r["sseqid"] for r in candidates
                      if math.isclose(r["bitscore"], top["bitscore"], abs_tol=1e-9)
                      and math.isclose(r["evalue"], top["evalue"], rel_tol=1e-9, abs_tol=1e-300)
                      and math.isclose(r["pident"], top["pident"], abs_tol=1e-9)})
    top["qseqid"] = query
    top["co_best_targets"] = ";".join(co_best)
    top["assignment_tie"] = len(co_best) > 1
    return top


def iter_best_hits(raw_tsv: Path, min_query_coverage: float,
                   annotations: dict[str, dict] | None = None,
                   ordinary_identity: float = DEFAULT_IDENTITY,
                   query_prefix: str = "", identity_rules: list[dict] | None = None,
                   use_default_exceptions: bool = True):
    """Stream one best qualifying target per query from DIAMOND's query-sorted output."""
    current_query = None
    candidates: list[dict] = []
    with raw_tsv.open() as handle:
        for line in handle:
            row = _parse_hit(line)
            if row is None:
                continue
            query = row["qseqid"]
            if current_query is not None and query != current_query:
                selected = _select_query_hit(query_prefix + current_query, candidates)
                if selected:
                    yield selected
                candidates = []
            current_query = query
            annotation = (annotations or {}).get(row["sseqid"], {})
            threshold, threshold_group = identity_threshold(
                row["sseqid"], row["stitle"] + " " + " ".join(str(v) for v in annotation.values()),
                ordinary_identity, identity_rules, use_default_exceptions)
            # DIAMOND's --query-cover is an early filter; repeat strict checks
            # here so the documented query-cover > threshold is explicit.
            if row["qcovhsp"] > min_query_coverage and row["pident"] >= threshold:
                row["identity_threshold"] = threshold
                row["identity_threshold_group"] = threshold_group
                candidates.append(row)
    if current_query is not None:
        selected = _select_query_hit(query_prefix + current_query, candidates)
        if selected:
            yield selected


def load_annotations(path: Path | None) -> dict[str, dict]:
    if path is None:
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "gene_id" not in reader.fieldnames:
            raise ValueError("Annotation TSV must have a 'gene_id' column")
        return {row["gene_id"]: row for row in reader if row.get("gene_id")}


def load_gene_lengths(path: Path | None) -> dict[str, int]:
    if path is None:
        return {}
    lengths = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"gene_id", "length_bp"} <= set(reader.fieldnames):
            raise ValueError("Gene-length TSV must have 'gene_id' and 'length_bp' columns")
        for row in reader:
            if row.get("gene_id"):
                length = int(row["length_bp"])
                if length < 1:
                    raise ValueError(f"Invalid length_bp for {row['gene_id']}: {length}")
                lengths[row["gene_id"]] = length
    return lengths


def gene_abundance(best_hits: Iterable[dict], total_reads: int,
                   reference_records: dict[str, tuple[int, str]],
                   annotations: dict[str, dict], include_zeroes: bool,
                   ordinary_identity: float = DEFAULT_IDENTITY,
                   nucleotide_lengths: dict[str, int] | None = None,
                   identity_rules: list[dict] | None = None,
                   use_default_exceptions: bool = True) -> tuple[list[dict], int]:
    counts: dict[str, int] = defaultdict(int)
    aa_lengths: dict[str, int] = {}
    titles: dict[str, str] = {}
    intervals: dict[str, bytearray] = {}
    identity_sums: dict[str, float] = defaultdict(float)
    query_cov_sums: dict[str, float] = defaultdict(float)
    ties: dict[str, int] = defaultdict(int)
    hit_count = 0
    for hit in best_hits:
        hit_count += 1
        gene = hit["sseqid"]
        counts[gene] += 1
        aa_length = int(reference_records.get(gene, (int(hit["slen"]), hit["stitle"]))[0])
        aa_lengths[gene] = aa_length
        titles[gene] = reference_records.get(gene, (aa_length, hit["stitle"]))[1]
        if gene not in intervals:
            intervals[gene] = bytearray(aa_length)
        start = max(0, min(hit["sstart"], hit["send"]) - 1)
        end = min(aa_length, max(hit["sstart"], hit["send"]))
        if end > start:
            intervals[gene][start:end] = b"\x01" * (end - start)
        identity_sums[gene] += hit["pident"]
        query_cov_sums[gene] += hit["qcovhsp"]
        ties[gene] += int(hit["assignment_tie"])

    genes = set(reference_records) if include_zeroes and reference_records else set(counts)
    rows = []
    for gene in sorted(genes):
        nreads = counts.get(gene, 0)
        aa_length = aa_lengths.get(gene, reference_records.get(gene, (0, ""))[0])
        length_bp = (nucleotide_lengths or {}).get(gene, aa_length * 3)
        length_source = "provided_nucleotide_length" if gene in (nucleotide_lengths or {}) else "protein_length_times_3"
        length_kb = length_bp / 1000.0
        rpkm = nreads / ((total_reads / 1_000_000.0) * length_kb) if total_reads and length_kb else 0.0
        breadth = (sum(intervals[gene]) / aa_length * 100
                   if aa_length and gene in intervals else 0.0)
        annotation = annotations.get(gene, {})
        threshold, threshold_group = identity_threshold(
            gene, titles.get(gene, "") + " " + " ".join(str(v) for v in annotation.values()),
            ordinary_identity, identity_rules, use_default_exceptions)
        rows.append({
            "gene_id": gene,
            "gene_title": titles.get(gene, reference_records.get(gene, (0, ""))[1]),
            "classification": annotation.get("classification", annotation.get("family", titles.get(gene, gene))),
            "subgroup": annotation.get("subgroup", ""),
            "family": annotation.get("family", ""),
            "hyddb_classification": annotation.get("hyddb_classification", ""),
            "gene_length_aa": aa_length,
            "gene_length_bp_used": length_bp,
            "gene_length_kb_used": f"{length_kb:.6f}",
            "gene_length_source": length_source,
            "read_count": nreads,
            "total_input_reads": total_reads,
            "identity_threshold_pct": f"{threshold:.1f}",
            "identity_threshold_group": threshold_group,
            "mean_read_identity_pct": f"{identity_sums.get(gene, 0.0) / nreads:.3f}" if nreads else "0",
            "mean_query_coverage_pct": f"{query_cov_sums.get(gene, 0.0) / nreads:.3f}" if nreads else "0",
            "gene_breadth_coverage_pct": f"{breadth:.3f}",
            "rpkm": f"{rpkm:.6f}",
            "tied_best_hit_reads": ties.get(gene, 0),
        })
    return rows, hit_count


def write_tsv(path: Path, rows: list[dict], default_fields: Iterable[str]):
    fields = list(rows[0]) if rows else list(default_fields)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Search metagenomic reads against protein gene references with DIAMOND blastx and report RPKM."
    )
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="One read FASTA/FASTQ, or paired-end R1 R2 FASTA/FASTQ files (optionally .gz).")
    parser.add_argument("--database", "--db", required=True, type=Path,
                        help="Prebuilt DIAMOND protein database (.dmnd or database prefix).")
    parser.add_argument("--reference-fasta", type=Path,
                        help="Protein FASTA used to build the database; supplies lengths and zero-hit rows.")
    parser.add_argument("--annotations", type=Path,
                        help="Optional TSV keyed by gene_id with classification, subgroup, family, hyddb_classification columns.")
    parser.add_argument("--gene-lengths", type=Path,
                        help="Optional TSV with gene_id and length_bp columns; overrides amino-acid length × 3 for RPKM.")
    parser.add_argument("-o", "--output-dir", required=True, type=Path)
    parser.add_argument("--sample-id", default="sample")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--query-coverage", type=float, default=DEFAULT_QUERY_COVERAGE,
                        help="Strict per-read DIAMOND query coverage threshold; only hits > this value pass (default: 80).")
    parser.add_argument("--minimum-identity", type=float, default=DEFAULT_IDENTITY,
                        help="Baseline identity for targets without a matching rule (default: 50).")
    parser.add_argument("--identity-thresholds", type=Path,
                        help="Optional CSV/TSV rules with pattern, minimum_identity, and optional match_type/group_name columns.")
    parser.add_argument("--no-default-identity-exceptions", action="store_true",
                        help="Disable the built-in paper-method 60%% group exceptions; custom rules remain active.")
    parser.add_argument("--minimum-gene-coverage", type=float, default=0.0,
                        help="Optional minimum union coverage of the protein target by passing reads (default: 0).")
    parser.add_argument("--minimum-reads", type=int, default=1,
                        help="Minimum assigned reads for a gene to be marked detected (default: 1).")
    parser.add_argument("--sensitivity", choices=("default", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"),
                        default="sensitive")
    parser.add_argument("--evalue", type=float, default=0.001)
    parser.add_argument("--max-target-seqs", type=int, default=25)
    parser.add_argument("--require-diamond-version", default=None,
                        help="Require the installed DIAMOND version string to contain this value, e.g. 2.1.6.")
    parser.add_argument("--include-zeroes", action="store_true",
                        help="Report every reference sequence at zero count/RPKM; requires --reference-fasta.")
    args = parser.parse_args(argv)

    if len(args.input) not in (1, 2):
        parser.error("--input takes one file or exactly two paired-end files")
    reads = [Path(value).expanduser() for value in args.input]
    if any(not path.is_file() for path in reads):
        parser.error("Every input read file must exist")
    if len(reads) == 2:
        counts_per_mate = [count_reads(path) for path in reads]
        if counts_per_mate[0] != counts_per_mate[1]:
            parser.error("Paired-end input files must contain the same number of read records")
        total_reads = sum(counts_per_mate)
    else:
        total_reads = count_reads(reads[0])
    if total_reads <= 0:
        parser.error("No input reads found")
    if args.threads < 1 or args.max_target_seqs < 1:
        parser.error("--threads and --max-target-seqs must be positive")
    if not 0 <= args.query_coverage < 100 or not 0 <= args.minimum_identity <= 100:
        parser.error("Coverage must be in [0,100) and identity in [0,100]")
    if args.minimum_gene_coverage < 0 or args.minimum_gene_coverage > 100:
        parser.error("--minimum-gene-coverage must be between 0 and 100")
    if args.minimum_reads < 1:
        parser.error("--minimum-reads must be at least 1")
    if args.include_zeroes and args.reference_fasta is None:
        parser.error("--include-zeroes requires --reference-fasta")
    if args.reference_fasta is not None and not args.reference_fasta.is_file():
        parser.error("--reference-fasta does not exist")
    if args.annotations is not None and not args.annotations.is_file():
        parser.error("--annotations does not exist")
    if args.gene_lengths is not None and not args.gene_lengths.is_file():
        parser.error("--gene-lengths does not exist")
    if args.identity_thresholds is not None and not args.identity_thresholds.is_file():
        parser.error("--identity-thresholds does not exist")
    if not args.database.is_file() and not Path(str(args.database) + ".dmnd").is_file():
        parser.error(f"DIAMOND database not found: {args.database}")
    diamond = shutil.which("diamond")
    if not diamond:
        parser.error("DIAMOND executable not found on PATH")
    try:
        version_result = subprocess.run([diamond, "version"], check=True, capture_output=True, text=True)
        diamond_version = (version_result.stdout or version_result.stderr).strip()
    except (OSError, subprocess.CalledProcessError) as exc:
        parser.error(f"Unable to determine DIAMOND version: {exc}")
    if args.require_diamond_version and args.require_diamond_version not in diamond_version:
        parser.error(f"DIAMOND {args.require_diamond_version} required; found: {diamond_version}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    refs = read_protein_fasta(args.reference_fasta) if args.reference_fasta else {}
    annotations = load_annotations(args.annotations)
    nucleotide_lengths = load_gene_lengths(args.gene_lengths)
    try:
        custom_identity_rules = load_identity_rules(args.identity_thresholds)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    default_rules = default_identity_rules() if not args.no_default_identity_exceptions else []
    identity_rules = custom_identity_rules + default_rules
    search_identity = min([args.minimum_identity] +
                          [rule["minimum_identity"] for rule in identity_rules])
    raw_paths = []
    commands = []
    for mate_index, read_path in enumerate(reads, start=1):
        suffix = f".R{mate_index}" if len(reads) == 2 else ""
        raw_path = args.output_dir / f"diamond_hits{suffix}.raw.tsv"
        cmd = [diamond, "blastx", "--query", str(read_path), "--db", str(args.database),
               "--out", str(raw_path), "--outfmt", "6", *FIELDS,
               "--id", str(search_identity),
               "--query-cover", str(args.query_coverage), "--max-target-seqs", str(args.max_target_seqs),
               "--evalue", str(args.evalue), "--threads", str(args.threads)]
        if args.sensitivity != "default":
            cmd.append("--" + args.sensitivity)
        commands.append(cmd)
        try:
            subprocess.run(cmd, check=True)
        except (OSError, subprocess.CalledProcessError) as exc:
            parser.error(f"DIAMOND search failed for {read_path}: {exc}")
        raw_paths.append(raw_path)

    assignment_path = args.output_dir / "best_read_assignments.tsv"
    assignment_fields = ("query_id", "gene_id", "gene_title", "identity_pct", "identity_threshold_pct",
                         "identity_threshold_group",
                         "query_coverage_pct", "subject_coverage_pct", "bitscore", "evalue",
                         "co_best_targets", "assignment_tie")

    def qualifying_hits():
        with assignment_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=assignment_fields, delimiter="\t")
            writer.writeheader()
            for mate_index, raw_path in enumerate(raw_paths, start=1):
                prefix = f"R{mate_index}:" if len(reads) == 2 else ""
                for hit in iter_best_hits(raw_path, args.query_coverage, annotations,
                                          args.minimum_identity, prefix, identity_rules,
                                          not args.no_default_identity_exceptions):
                    writer.writerow({"query_id": hit["qseqid"], "gene_id": hit["sseqid"],
                                     "gene_title": hit["stitle"], "identity_pct": f"{hit['pident']:.3f}",
                                     "identity_threshold_pct": f"{hit['identity_threshold']:.1f}",
                                     "identity_threshold_group": hit["identity_threshold_group"],
                                     "query_coverage_pct": f"{hit['qcovhsp']:.3f}",
                                     "subject_coverage_pct": f"{hit['scovhsp']:.3f}",
                                     "bitscore": f"{hit['bitscore']:.3f}", "evalue": hit["evalue"],
                                     "co_best_targets": hit["co_best_targets"],
                                     "assignment_tie": int(hit["assignment_tie"])})
                    yield hit

    abundance, passing_read_count = gene_abundance(qualifying_hits(), total_reads, refs,
                                                    annotations, args.include_zeroes,
                                                    args.minimum_identity, nucleotide_lengths,
                                                    identity_rules, not args.no_default_identity_exceptions)
    for row in abundance:
        row["sample_id"] = args.sample_id
        row["detected"] = int(int(row["read_count"]) >= args.minimum_reads
                              and float(row["gene_breadth_coverage_pct"]) >= args.minimum_gene_coverage)
        row["minimum_reads_threshold"] = args.minimum_reads
        row["minimum_gene_coverage_pct_threshold"] = f"{args.minimum_gene_coverage:.2f}"
    write_tsv(args.output_dir / "gene_abundance_rpkm.tsv", abundance,
              ("sample_id", "gene_id", "gene_title", "classification", "subgroup", "family",
               "hyddb_classification", "gene_length_aa", "gene_length_bp_used", "gene_length_kb_used",
               "gene_length_source",
               "read_count", "total_input_reads", "identity_threshold_pct", "identity_threshold_group",
               "mean_read_identity_pct",
               "mean_query_coverage_pct", "gene_breadth_coverage_pct", "rpkm", "tied_best_hit_reads",
               "detected", "minimum_reads_threshold", "minimum_gene_coverage_pct_threshold"))
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sample_id": args.sample_id,
        "diamond_version": diamond_version,
        "search_mode": "blastx",
        "search_strategy": "Each read mate searched independently; one best threshold-passing target assigned per read.",
        "inputs": [str(path.resolve()) for path in reads],
        "read_counts_per_input": [count_reads(path) for path in reads],
        "total_input_reads": total_reads,
        "database": str(args.database.resolve()),
        "reference_fasta": str(args.reference_fasta.resolve()) if args.reference_fasta else None,
        "gene_lengths_tsv": str(args.gene_lengths.resolve()) if args.gene_lengths else None,
        "query_coverage_threshold_pct_strictly_greater_than": args.query_coverage,
        "baseline_identity_threshold_pct": args.minimum_identity,
        "identity_thresholds_file": str(args.identity_thresholds.resolve()) if args.identity_thresholds else None,
        "identity_threshold_rules": identity_rules,
        "default_identity_exceptions_enabled": not args.no_default_identity_exceptions,
        "diamond_prefilter_identity_pct": search_identity,
        "minimum_gene_breadth_coverage_pct_for_detected": args.minimum_gene_coverage,
        "minimum_assigned_reads_for_detected": args.minimum_reads,
        "rpkm_formula": "read_count / (total_input_reads / 1,000,000) / gene_length_kb",
        "protein_length_to_coding_bp_fallback": "amino_acid_length * 3",
        "commands": commands,
        "outputs": ["diamond_hits*.raw.tsv", "best_read_assignments.tsv", "gene_abundance_rpkm.tsv"],
    }
    (args.output_dir / "metagenomic_run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Sample: {args.sample_id}")
    print(f"Input reads counted: {total_reads}")
    print(f"Reads with a threshold-passing best target: {passing_read_count}")
    print("Raw DIAMOND alignments: " + ", ".join(map(str, raw_paths)))
    print(f"Best read assignments: {args.output_dir / 'best_read_assignments.tsv'}")
    print(f"RPKM abundance table: {args.output_dir / 'gene_abundance_rpkm.tsv'}")
    print(f"Run metadata: {args.output_dir / 'metagenomic_run_manifest.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
