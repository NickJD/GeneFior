"""Create SNP, allele-evidence, and trimming figures for QC-passing genes."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
try:
    from Bio.Align import PairwiseAligner
except ImportError:  # pragma: no cover - reported by CLI dependencies
    PairwiseAligner = None

if __package__:
    from .GeneFior_Architect import align_to_contig, find_assemblies, make_aligner, read_fasta
else:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from GeneFior.GeneFior_Architect import align_to_contig, find_assemblies, make_aligner, read_fasta


ACCEPTED = {"PASS", "LIKELY_PASS"}
BASE_ORDER = {"A": 0, "C": 1, "G": 2, "T": 3}
BASE_COLORS = {"A": "#2471A3", "C": "#239B56", "G": "#D4AC0D", "T": "#C0392B"}


def read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_references(path: Path) -> dict[str, str]:
    records = read_fasta(path)
    safe = lambda name: re.sub(r"[^\w\-.]", "_", name)
    return {safe(name): seq for name, seq in records}


def qc_pass_reconstructions(report: Path) -> list[dict]:
    rows = read_tsv(report)
    groups = defaultdict(list)
    for row in rows:
        if row.get("record", "").endswith("sample_with_insertions"):
            continue
        mapper = Path(row["reconstruction_file"]).parent.name.lower()
        row["mapper"] = mapper
        key = (row["sample_id"], row["gene"], mapper, row["reconstruction_file"], row["record"])
        groups[key].append(row)
    selected = []
    for group in groups.values():
        by_assembly = {row["assembly_kind"]: row for row in group}
        if {"long_read", "enhanced_contig"} <= by_assembly.keys() and all(
            by_assembly[k]["status"] in ACCEPTED for k in ("long_read", "enhanced_contig")
        ):
            selected.append(by_assembly["long_read"])
            selected.append(by_assembly["enhanced_contig"])
    return selected


def _extract_truth_sequence(row: dict, assembly_path: Path) -> Optional[str]:
    """Extract the reference-locus span from an assembly in reference orientation."""
    contig = row.get("truth_contig", "")
    start = row.get("truth_reference_locus_start_1based", "")
    end = row.get("truth_reference_locus_end_1based", "")
    if not contig or not start or not end:
        return None
    sequence = dict(read_fasta(assembly_path)).get(contig)
    if sequence is None:
        return None
    start0, end0 = int(start) - 1, int(end)
    if start0 < 0 or end0 > len(sequence) or end0 <= start0:
        return None
    result = sequence[start0:end0]
    if row.get("truth_reference_strand") == "-":
        result = result.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    return result


def _alignment_differences(reference: str, query: str, aligner) -> list[dict]:
    """Project global alignment differences to 1-based reference coordinates."""
    alignment = aligner.align(reference, query)[0]
    coords = alignment.coordinates
    differences = []
    ref_pos = 0
    query_pos = 0
    for i in range(coords.shape[1] - 1):
        r0, r1 = int(coords[0, i]), int(coords[0, i + 1])
        q0, q1 = int(coords[1, i]), int(coords[1, i + 1])
        dr, dq = r1 - r0, q1 - q0
        if dr and dq:
            ref_chunk, query_chunk = reference[r0:r1], query[q0:q1]
            common = min(dr, dq)
            for j in range(common):
                if ref_chunk[j] != query_chunk[j]:
                    differences.append({"position_1based": r0 + j + 1, "reference": ref_chunk[j],
                                        "observed": query_chunk[j],
                                        "type": "SNP" if ref_chunk[j] in "ACGT" and query_chunk[j] in "ACGT" else "AMBIGUITY"})
            if dr > common:
                differences.append({"position_1based": r0 + common + 1,
                                    "reference": reference[r0 + common:r1], "observed": "-",
                                    "type": "DELETION"})
            elif dq > common:
                differences.append({"position_1based": r0 + common,
                                    "reference": "-", "observed": query[q0 + common:q1],
                                    "type": "INSERTION"})
        elif dr:
            differences.append({"position_1based": r0 + 1, "reference": reference[r0:r1],
                                "observed": "-", "type": "DELETION"})
        elif dq:
            differences.append({"position_1based": r0, "reference": "-",
                                "observed": query[q0:q1], "type": "INSERTION"})
        ref_pos, query_pos = r1, q1
    return differences


def plot_three_way_comparison(reference: str, sequences: dict[str, str], out: Path,
                              title: str, metadata: dict) -> list[dict]:
    """Write a three-source variation map and its tabular site evidence."""
    if PairwiseAligner is None:
        raise RuntimeError("Three-way alignment figures require Biopython.")
    aligner = PairwiseAligner(mode="global")
    aligner.match_score, aligner.mismatch_score = 2, -3
    aligner.open_gap_score, aligner.extend_gap_score = -5, -1
    tracks = {name: _alignment_differences(reference, seq, aligner)
              for name, seq in sequences.items()}
    colors = {"reconstruction": "#C0392B", "long_read_ground_truth": "#2471A3",
              "enhanced_contig_ground_truth": "#239B56"}
    fig, ax = plt.subplots(figsize=(15, 4.8))
    lanes = {name: i + 1 for i, name in enumerate(sequences)}
    ax.hlines(0, 1, len(reference), color="#555555", linewidth=2, label="Reference")
    details = []
    for name, diffs in tracks.items():
        y = lanes[name]
        ax.hlines(y, 1, len(reference), color="#dddddd", linewidth=1)
        previous_pos = None
        label_stack = 0
        for item in diffs:
            pos = max(1, min(len(reference), int(item["position_1based"])))
            label_stack = label_stack + 1 if previous_pos is not None and pos - previous_pos <= 8 else 0
            previous_pos = pos
            color = colors[name]
            marker = {"SNP": "o", "DELETION": "v", "INSERTION": "^", "AMBIGUITY": "D"}[item["type"]]
            ax.scatter(pos, y, marker=marker, s=38, color=color, edgecolor="white",
                       linewidth=0.45, zorder=3)
            if item["type"] in {"SNP", "AMBIGUITY"}:
                label = f"{item['reference']}>{item['observed']}"
            elif item["type"] == "INSERTION":
                label = f"+{item['observed']}"
            else:
                label = f"-{item['reference']}"
            ax.annotate(label, (pos, y), xytext=(2, 5 + label_stack * 8), textcoords="offset points",
                        fontsize=6.5, rotation=0, ha="left", va="bottom", color=color)
            details.append({"source": name, **item})
        ax.text(len(reference) * 1.09, y,
                f"{len(diffs)} differences", va="center", fontsize=9, color=colors[name])
    ax.set_xlim(0, len(reference) * 1.22)
    ax.set_ylim(-0.6, len(lanes) + 0.65)
    ax.set_yticks([0, *lanes.values()], ["Reference", *lanes.keys()])
    ax.set_xlabel("Reference gene coordinate (1-based; indels are anchored to the nearest reference base)")
    ax.set_title(title, fontsize=12, pad=13)
    ax.grid(axis="x", color="#e6e6e6", linewidth=0.55)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(0.01, -0.27,
            f"Reconstruction length: {metadata['reconstruction_length']} bp · "
            f"Long-read locus: {metadata['long_read_locus']} · "
            f"Enhanced-contig locus: {metadata['enhanced_contig_locus']} · "
            "circles = substitutions, diamonds = ambiguous calls, triangles = indels",
            transform=ax.transAxes, fontsize=8.5)
    fig.subplots_adjust(left=0.20, right=0.84, top=0.88, bottom=0.25)
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return details


def generate_three_way_visuals(report: Path, references: dict[str, str], output_dir: Path) -> list[dict]:
    """Create one reference/reconstruction/assembly figure per high-quality consensus."""
    rows = read_tsv(report)
    grouped = defaultdict(dict)
    for row in rows:
        if row.get("record", "").endswith("sample_with_insertions"):
            continue
        row["mapper"] = Path(row["reconstruction_file"]).parent.name.lower()
        key = (row["sample_id"], row["gene"], row["mapper"], row["record"], row["reconstruction_file"])
        grouped[key][row["assembly_kind"]] = row
    manifest, summary = [], []
    for (sample, gene, mapper, record, recon_path), assembly_rows in sorted(grouped.items()):
        if not {"long_read", "enhanced_contig"} <= assembly_rows.keys():
            continue
        lr, enhanced = assembly_rows["long_read"], assembly_rows["enhanced_contig"]
        if any(r.get("status") not in ACCEPTED or float(r.get("n_pct") or 100) >= 1.0
               for r in (lr, enhanced)):
            continue
        reference = reference_for_gene(gene, references)
        if not reference:
            continue
        recon_records = read_fasta(Path(recon_path))
        if not recon_records:
            continue
        reconstruction = recon_records[0][1]
        assemblies = {}
        for row, key, label in ((lr, "long_read_ground_truth", "long_read"),
                                (enhanced, "enhanced_contig_ground_truth", "enhanced_contig")):
            seq = _extract_truth_sequence(row, Path(row["assembly_file"]))
            if seq is None:
                assemblies = {}
                break
            assemblies[key] = seq
        if len(assemblies) != 2:
            continue
        sequences = {"reconstruction": reconstruction, **assemblies}
        case_dir = output_dir / sample / gene
        case_dir.mkdir(parents=True, exist_ok=True)
        stem = f"{mapper}_{Path(recon_path).stem}_three_way"
        png = case_dir / f"{stem}.png"
        metadata = {"reconstruction_length": len(reconstruction),
                    "long_read_locus": f"{lr['truth_contig']}:{lr['truth_reference_locus_start_1based']}-{lr['truth_reference_locus_end_1based']}",
                    "enhanced_contig_locus": f"{enhanced['truth_contig']}:{enhanced['truth_reference_locus_start_1based']}-{enhanced['truth_reference_locus_end_1based']}"}
        diffs = plot_three_way_comparison(reference, sequences, png,
                                          f"{sample} · {gene} · {mapper} · three-way variation",
                                          metadata)
        tsv = case_dir / f"{stem}_differences.tsv"
        write_tsv(diffs, tsv)
        summary.append({"sample_id": sample, "gene": gene, "mapper": mapper,
                        "reconstruction_file": recon_path, "figure": str(png),
                        "differences_tsv": str(tsv), "reconstruction_length": len(reconstruction),
                        "reconstruction_n_pct": lr.get("n_pct", ""),
                        "long_read_status": lr["status"], "enhanced_contig_status": enhanced["status"],
                        "reconstruction_snp_count": sum(x["source"] == "reconstruction" and x["type"] == "SNP" for x in diffs),
                        "long_read_snp_count": sum(x["source"] == "long_read_ground_truth" and x["type"] == "SNP" for x in diffs),
                        "enhanced_contig_snp_count": sum(x["source"] == "enhanced_contig_ground_truth" and x["type"] == "SNP" for x in diffs)})
        manifest.extend([png, tsv])
    summary_path = output_dir / "three_way_visualisation_summary.tsv"
    write_tsv(summary, summary_path)
    return summary


def reference_for_gene(gene: str, refs: dict[str, str]) -> Optional[str]:
    if gene in refs:
        return refs[gene]
    return next((seq for name, seq in refs.items() if name == gene or name.startswith(gene)), None)


def plot_snp_profiles(records: list[dict], references: dict[str, str], out: Path):
    groups = defaultdict(dict)
    for row in records:
        groups[(row["sample_id"], row["gene"])][row["mapper"]] = Path(row["reconstruction_file"])
    genes = sorted(groups)
    ncols = 2
    nrows = math.ceil(len(genes) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, max(5, nrows * 2.35)), squeeze=False)
    fig.suptitle("QC-passing SNP profiles against reference genes", fontsize=16, y=0.995)
    for ax, (sample, gene) in zip(axes.flat, genes):
        ref = reference_for_gene(gene, references)
        gene_len = len(ref) if ref else 1
        mapper_paths = groups[(sample, gene)]
        y_positions = {mapper: i for i, mapper in enumerate(sorted(mapper_paths))}
        n_snps = 0
        for mapper, reconstruction in sorted(mapper_paths.items()):
            variant_path = reconstruction.parent / f"{reconstruction.stem.replace('_consensus', '')}_variants.tsv"
            variants = [v for v in read_tsv(variant_path) if v.get("Type") == "SNP"]
            y = y_positions[mapper]
            for v in variants:
                pos = int(v["Position"])
                alt = v["Sample_Base"]
                ref_base = v["Ref_Base"]
                depth = max(1, int(float(v.get("Depth", 1) or 1)))
                color = BASE_COLORS.get(alt, "#555555")
                ax.vlines(pos, y, y + 0.30, color=color, linewidth=1.3, alpha=0.85)
                ax.scatter(pos, y + 0.30, s=min(95, 24 + depth * 0.20), color=color,
                           edgecolor="black", linewidth=0.25, zorder=3)
                freq = float(v.get("Freq_pct", 0) or 0)
                ax.annotate(f"{ref_base}>{alt} {freq:.0f}%", (pos, y + 0.30), xytext=(2, 4),
                            textcoords="offset points", fontsize=7, rotation=35, ha="left")
                n_snps += 1
        ax.set_xlim(0, gene_len)
        ax.set_ylim(-0.30, max(0.9, len(y_positions) - 0.1))
        ax.set_yticks(list(y_positions.values()), [m.capitalize() for m in y_positions])
        ax.set_xlabel("Reference coordinate (bp)")
        ax.set_title(f"{sample} · {gene} · {n_snps} mapper-SNP calls", fontsize=10)
        ax.grid(axis="x", color="#d8d8d8", linewidth=0.55)
        ax.spines[["top", "right"]].set_visible(False)
        if n_snps == 0:
            ax.text(0.5, 0.5, "No callable SNPs", transform=ax.transAxes,
                    ha="center", va="center", color="#555555")
    for ax in axes.flat[len(genes):]:
        ax.remove()
    fig.text(0.5, 0.01, "Each stem is a consensus substitution; label shows alternate-allele fraction and marker size scales with depth. Different bases at the same coordinate across mappers are flagged in the TSV.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.035, 1, 0.97))
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


def count_reference_loci(reference: str, assembly_path: Path, aligner,
                         min_identity: float = 90.0, min_coverage: float = 80.0) -> list[dict]:
    loci = []
    for contig_name, sequence in read_fasta(assembly_path):
        masked = list(sequence)
        for _ in range(12):
            hit = align_to_contig(reference, "".join(masked), aligner)
            if (hit is None or hit.identity_pct < min_identity
                    or hit.query_coverage_pct < min_coverage
                    or hit.target_end <= hit.target_start):
                break
            loci.append({"contig": contig_name, "identity_pct": hit.identity_pct,
                         "coverage_pct": hit.query_coverage_pct,
                         "start": hit.target_start + 1, "end": hit.target_end,
                         "strand": hit.strand})
            pad = max(25, len(reference) // 20)
            left = max(0, hit.target_start - pad)
            right = min(len(masked), hit.target_end + pad)
            masked[left:right] = ["N"] * (right - left)
    return loci


def read_validation(validation_path: Path) -> tuple[str, int, int]:
    text = validation_path.read_text(errors="replace") if validation_path.is_file() else ""
    state = re.search(r"Classification:\s*(\w+)", text)
    sites = re.search(r"Informative sites:\s*(\d+)", text)
    hap_fasta = validation_path.with_name(validation_path.name.replace("_validation.txt", "_haplotypes.fasta"))
    haplotypes = sum(1 for name, _ in read_fasta(hap_fasta)) if hap_fasta.is_file() else 0
    return state.group(1).lower() if state else "unknown", int(sites.group(1)) if sites else 0, haplotypes


def family_for_candidate(gene: str, rows: list[dict]) -> set[str]:
    families = set()
    safe = lambda name: re.sub(r"[^\w\-.]", "_", name)
    for row in rows:
        top = row.get("Top_Database_Candidate", "")
        competing = set(filter(None, row.get("Competing_Alleles", "").split(";")))
        if safe(gene) == safe(top) or safe(gene) in {safe(item) for item in competing}:
            families.add(row.get("Family", ""))
    return families


def build_allele_summary(passing: list[dict], amrfior_root: Path,
                         truth_root: Path, references: dict[str, str], aligner):
    unique = {}
    for row in passing:
        unique[(row["sample_id"], row["gene"], row["mapper"], row["reconstruction_file"])] = row
    group_keys = sorted({(r["sample_id"], r["gene"]) for r in unique.values()})
    summary, locus_cache = [], {}
    for sample, gene in group_keys:
        evidence_rows = read_tsv(amrfior_root / sample / "resfinder_allele_resolution.tsv")
        families = family_for_candidate(gene, evidence_rows)
        relevant = [r for r in evidence_rows if r.get("Family") in families]
        direct_status = {"CANDIDATE_ALLELE_DETECTED", "EXACT_ALLELE_DETECTED", "EXACT_PROTEIN_DETECTED"}
        top_candidates = sorted({r.get("Top_Database_Candidate", "") for r in relevant
                                 if r.get("Top_Database_Candidate") and r.get("Status") in direct_status})
        competing = sorted({allele for r in relevant for allele in r.get("Competing_Alleles", "").split(";") if allele})
        mapper_rows = [r for r in unique.values() if r["sample_id"] == sample and r["gene"] == gene]
        mv = [read_validation(Path(r["reconstruction_file"]).with_name(
            Path(r["reconstruction_file"]).name.replace("_consensus.fasta", "_validation.txt"))) for r in mapper_rows]
        mv_state = max((x[0] for x in mv), key=lambda x: {"single": 0, "uncertain": 1, "multi": 2}.get(x, -1), default="unknown")
        informative = max((x[1] for x in mv), default=0)
        hap_count = max((x[2] for x in mv), default=0)
        ref = reference_for_gene(gene, references)
        asm_hits = {}
        for kind, path in find_assemblies(truth_root, sample):
            # The same assembly is queried for several distinct gene references.
            key = (str(path), ref)
            if key not in locus_cache:
                locus_cache[key] = count_reference_loci(ref, path, aligner) if ref else []
            asm_hits[kind] = locus_cache[key]
        long_loci = asm_hits.get("long_read", [])
        enhanced_loci = asm_hits.get("enhanced_contig", [])
        summary.append({
            "sample_id": sample, "gene": gene,
            "family": ";".join(sorted(families)),
            "amrfior_top_candidate_labels": ";".join(top_candidates),
            "n_distinct_top_candidate_labels": len(top_candidates),
            "competing_reference_labels": ";".join(competing),
            "n_competing_reference_labels": len(competing),
            "reconstruction_multi_version": mv_state,
            "informative_read_sites_max": informative,
            "haplotype_sequences_emitted_max": hap_count,
            "copy_number_estimate": "not_inferred",
            "allele_evidence_interpretation": (
                "Read haplotypes are evidence for sequence versions; assembly hits are candidate loci. "
                "Neither alone establishes allele dosage."
            ),
            "long_read_reference_like_loci": len(long_loci),
            "enhanced_contig_reference_like_loci": len(enhanced_loci),
            "long_read_loci": ";".join(f"{x['contig']}:{x['start']}-{x['end']}({x['identity_pct']:.1f}%)" for x in long_loci),
            "enhanced_contig_loci": ";".join(f"{x['contig']}:{x['start']}-{x['end']}({x['identity_pct']:.1f}%)" for x in enhanced_loci),
        })
    return summary


def plot_allele_evidence(summary: list[dict], out: Path):
    cols = [
        ("n_distinct_top_candidate_labels", "AMRFior\ntop allele IDs"),
        ("n_competing_reference_labels", "Competing\nreference IDs"),
        ("informative_read_sites_max", "Read informative\nsites"),
        ("long_read_reference_like_loci", "Long-read\nloci"),
        ("enhanced_contig_reference_like_loci", "Enhanced\nloci"),
    ]
    labels = [f"{r['sample_id']} · {r['gene']}" for r in summary]
    data = np.array([[int(r[key]) for key, _ in cols] for r in summary], dtype=float)
    fig, ax = plt.subplots(figsize=(12, max(5, len(summary) * 0.47)))
    im = ax.imshow(data, aspect="auto", cmap="Blues", interpolation="nearest")
    ax.set_xticks(range(len(cols)), [name for _, name in cols])
    ax.set_yticks(range(len(labels)), labels)
    ax.set_title("Allele and genome-locus evidence for QC-passing reconstructions", pad=14)
    ax.set_xlabel("Counts (evidence types are distinct; they are not interchangeable copy-number estimates)")
    ax.tick_params(axis="x", labelrotation=0)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = int(data[i, j])
            ax.text(j, i, str(val), ha="center", va="center",
                    color="white" if val >= max(1, data.max() * 0.55) else "#222222", fontsize=9)
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("Count")
    ax.set_xticks(np.arange(-0.5, data.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, data.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.4)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_boundaries(records: list[dict], references: dict[str, str], out: Path):
    by_gene = defaultdict(dict)
    for row in records:
        key = (row["sample_id"], row["gene"])
        by_gene[key][(row["mapper"], row["assembly_kind"])] = row
    genes = sorted(by_gene)
    ncols = 2
    nrows = math.ceil(len(genes) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, max(5, nrows * 2.3)), squeeze=False)
    fig.suptitle("Matched reconstruction boundaries relative to the reference gene", fontsize=15, y=0.995)
    assembly_colors = {"long_read": "#2874A6", "enhanced_contig": "#D35400"}
    for ax, (sample, gene) in zip(axes.flat, genes):
        reference = reference_for_gene(gene, references)
        length = len(reference) if reference else 1
        entries = by_gene[(sample, gene)]
        mappers = sorted({k[0] for k in entries})
        lane = 0
        ylabels = []
        for mapper in mappers:
            for assembly_kind in ("long_read", "enhanced_contig"):
                row = entries.get((mapper, assembly_kind))
                if not row:
                    continue
                d5 = row.get("gene_5prime_extension_bp")
                d3 = row.get("gene_3prime_extension_bp")
                if d5 in (None, "") or d3 in (None, "") or row.get("boundary_strands_match") != "True":
                    continue
                x0 = -float(d5) / length * 100
                x1 = 100 + float(d3) / length * 100
                color = assembly_colors[assembly_kind]
                ax.barh(lane, 100, left=0, height=0.48, color="#e4e4e4", zorder=1)
                ax.barh(lane, x1 - x0, left=x0, height=0.27, color=color, alpha=0.88, zorder=2)
                ax.plot([x0, x1], [lane, lane], "|", color="#222222", markersize=7, zorder=3)
                ax.text(x0, lane - 0.23, f"{int(d5):+d}", ha="center", va="top", fontsize=7, color=color)
                ax.text(x1, lane - 0.23, f"{int(d3):+d}", ha="center", va="top", fontsize=7, color=color)
                ylabels.append(f"{mapper} · {'LR' if assembly_kind == 'long_read' else 'enhanced'}")
                lane += 1
        ax.axvline(0, color="#555555", linewidth=0.75)
        ax.axvline(100, color="#555555", linewidth=0.75)
        ax.axvline(-20, color="#888888", linewidth=0.7, linestyle="--")
        ax.axvline(120, color="#888888", linewidth=0.7, linestyle="--")
        ax.set_xlim(-35, 135)
        ax.set_yticks(range(len(ylabels)), ylabels)
        ax.invert_yaxis()
        ax.set_xlabel("Position relative to reference gene (%)")
        ax.set_title(f"{sample} · {gene}", fontsize=10)
        ax.grid(axis="x", color="#dddddd", linewidth=0.55)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        if not ylabels:
            ax.text(0.5, 0.5, "No matching boundary coordinates", transform=ax.transAxes,
                    ha="center", va="center")
    for ax in axes.flat[len(genes):]:
        ax.remove()
    fig.text(0.5, 0.01,
             "Gray bar = reference span (0–100%); colored bar = reconstruction span matched to the assembly. Dashed lines = current ±20% extension search limit. Endpoint labels are signed bp relative to the reference boundaries.",
             ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot SNP, allele evidence, and boundary behavior for assembly-QC-passing reconstructions.")
    parser.add_argument("--architect-report", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--reconstruction-dir", required=True, type=Path)
    parser.add_argument("--amrfior-root", required=True, type=Path)
    parser.add_argument("--truth-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    if not args.architect_report.is_file():
        parser.error(f"Architect report not found: {args.architect_report}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    passing = qc_pass_reconstructions(args.architect_report)
    if not passing:
        parser.error("No consensus reconstruction passes against both long-read and enhanced-contig assemblies.")
    references = load_references(args.reference_fasta)
    aligner = make_aligner()
    stem = args.output_dir
    three_way = generate_three_way_visuals(args.architect_report, references, stem / "three_way")
    plot_snp_profiles(passing, references, stem / "reconstruction-snp-profiles.png")
    snp_evidence = collect_snp_evidence(passing)
    write_tsv(snp_evidence, stem / "reconstruction-snp-evidence.tsv")
    evidence = build_allele_summary(passing, args.amrfior_root, args.truth_root, references, aligner)
    plot_allele_evidence(evidence, stem / "reconstruction-allele-evidence.png")
    plot_boundaries(passing, references, stem / "reconstruction-boundaries.png")
    if evidence:
        with (stem / "reconstruction-allele-evidence.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(evidence[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(evidence)
    write_manifest(args, passing, evidence, stem, three_way)
    print(f"QC-passing mapper/assembly comparisons: {len(passing)}")
    print(f"QC-passing sample/gene groups: {len(evidence)}")
    print(f"High-quality three-way reconstruction figures: {len(three_way)}")
    print(f"Figures written under: {stem}")
    return 0


def collect_snp_evidence(records: list[dict]) -> list[dict]:
    """Summarize each callable SNP position across passing mapper consensuses."""
    by_site = defaultdict(dict)
    for row in records:
        path = Path(row["reconstruction_file"])
        variant_path = path.parent / f"{path.stem.replace('_consensus', '')}_variants.tsv"
        for variant in read_tsv(variant_path):
            if variant.get("Type") != "SNP":
                continue
            key = (row["sample_id"], row["gene"], int(variant["Position"]))
            by_site[key][row["mapper"]] = variant
    output = []
    for (sample, gene, position), calls in sorted(by_site.items()):
        alts = sorted({v.get("Sample_Base", "") for v in calls.values()})
        output.append({
            "sample_id": sample, "gene": gene, "reference_position_1based": position,
            "mapper_calls": ";".join(f"{m}:{v.get('Ref_Base','')}>{v.get('Sample_Base','')}"
                                      for m, v in sorted(calls.items())),
            "mapper_support": len(calls), "n_distinct_alternate_bases": len(alts),
            "mapper_disagreement": len(alts) > 1,
            "depths": ";".join(f"{m}:{v.get('Depth','')}" for m, v in sorted(calls.items())),
            "alternate_fractions_pct": ";".join(f"{m}:{v.get('Freq_pct','')}" for m, v in sorted(calls.items())),
        })
    return output


def write_tsv(rows: list[dict], path: Path):
    fields = list(rows[0]) if rows else ["sample_id", "gene", "reference_position_1based"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(args, passing: list[dict], evidence: list[dict], output_dir: Path,
                   three_way: list[dict]):
    sources = {"architect_report": args.architect_report, "reference_fasta": args.reference_fasta}
    for row in passing:
        reconstruction = Path(row["reconstruction_file"])
        sources[str(reconstruction)] = reconstruction
        prefix = reconstruction.stem.replace("_consensus", "")
        for suffix in ("_variants.tsv", "_validation.txt", "_depth.tsv"):
            candidate = reconstruction.with_name(prefix + suffix)
            if candidate.is_file():
                sources[str(candidate)] = candidate
    for sample in sorted({r["sample_id"] for r in passing}):
        table = args.amrfior_root / sample / "resfinder_allele_resolution.tsv"
        if table.is_file():
            sources[str(table)] = table
        for _, assembly in find_assemblies(args.truth_root, sample):
            sources[str(assembly)] = assembly
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "consensus reconstructions passing against both assembly sources",
        "qc_passing_mapper_assembly_comparisons": len(passing),
        "qc_passing_sample_gene_groups": len(evidence),
        "three_way_high_quality_reconstructions": len(three_way),
        "inputs": [{"path": str(path), "size_bytes": path.stat().st_size,
                    "modified_time_ns": path.stat().st_mtime_ns}
                   for path in sorted(set(sources.values())) if path.is_file()],
        "outputs": ["reconstruction-snp-profiles.png", "reconstruction-snp-evidence.tsv",
                    "reconstruction-allele-evidence.png", "reconstruction-allele-evidence.tsv",
                    "reconstruction-boundaries.png", "three_way/three_way_visualisation_summary.tsv",
                    "per-reconstruction PNG and differences TSV files under three_way/sample/gene/"],
    }
    (output_dir / "reconstruction-visualisation-manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )


if __name__ == "__main__":
    sys.exit(main())
