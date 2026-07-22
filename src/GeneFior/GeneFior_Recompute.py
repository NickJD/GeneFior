import argparse
import sys
import os
from pathlib import Path
import logging
import json
import tempfile
from datetime import datetime
from collections import defaultdict

try:
    from .constants import GENEFIOR_VERSION
    from .read_store import ReadMatchStore
    from .workflow import Workflow
    from .utils import load_gene_id_file
    from .hamronization import (
        add_hamronization_arguments,
        export_from_options,
        infer_sample_id,
        parse_database_versions,
        validate_database_versions,
        validate_hamronization_request,
    )
except (ModuleNotFoundError, ImportError):
    from constants import GENEFIOR_VERSION
    from read_store import ReadMatchStore
    from workflow import Workflow
    from utils import load_gene_id_file
    from hamronization import (
        add_hamronization_arguments,
        export_from_options,
        infer_sample_id,
        parse_database_versions,
        validate_database_versions,
        validate_hamronization_request,
    )



def discover_files(options, logger, databases_found, tools_found):

    #databases_found.add(database)
    #tools_found[database].add(tool)

    # Discover all alignment output files in the input directory.
    logger.info("Discovering alignment output files...")

    if str(options.input).endswith('raw_outputs'):
        raw_dir = Path(options.input)
    else:
        raw_dir = Path(options.input) / "raw_outputs"
    if not raw_dir.exists():
        logger.error(f"Raw outputs directory not found: {raw_dir}")
        return False

    # Pattern matching for different file types
    all_patterns = {
        'blastn': '*_blastn_results.tsv',
        'blastp': '*_blastp_results.tsv',
        'blastx': '*_blastx_results.tsv',
        'diamond': '*_diamond_results.tsv',
        'bowtie2': '*_bowtie2_results_sorted.bam',
        'bwa': '*_bwa_results_sorted.bam',
        'minimap2': '*_minimap2_results_sorted.bam'
    }
    if 'all' in options.tools: # User selection
        patterns = all_patterns
    else:
        patterns = {tool: all_patterns[tool] for tool in options.tools if tool in all_patterns}

    found_files = defaultdict(list)

    for tool, pattern in patterns.items():
        files = list(raw_dir.glob(pattern))
        for file in files:
            # Extract database name from filename
            # Format: {database}_{tool}_results...
            filename = file.stem.replace('_results', '').replace('_sorted', '')

            if tool in ['bowtie2', 'bwa', 'minimap2']:
                database = filename.replace(f'_{tool}', '')
            elif tool == 'blastn':
                database = filename.replace('_blastn', '')
            elif tool == 'blastx':
                database = filename.replace('_blastx', '')
            elif tool == 'blastp':
                database = filename.replace('_blastp', '')
            elif tool == 'diamond':
                database = filename.replace('_diamond', '')

            found_files[database].append((tool, file))
            databases_found.add(database)
            tools_found[database].add(tool)

    if not found_files:
        logger.error("No alignment output files found!")
        return False

    logger.info(f"Found {len(databases_found)} databases: {', '.join(databases_found)}")
    for db in databases_found:
        logger.info(f"  {db}: {', '.join(sorted(tools_found[db]))}")


    return found_files, databases_found, tools_found


def detect_diamond_mode(file_path: Path, source_parameters=None) -> str:
    """Try to infer whether a DIAMOND output was produced using blastx or blastp
    by scanning the header/comment lines of the output file. Returns one of
    'DIAMOND-BLASTX', 'DIAMOND-BLASTP' or 'DIAMOND' (unknown)."""
    source_parameters = source_parameters or {}
    source_sequence_type = source_parameters.get('sequence_type')
    source_genes_type = source_parameters.get('genes_type')
    if source_sequence_type == 'Genes-FASTA':
        if source_genes_type in ('aa', 'protein'):
            return 'DIAMOND-BLASTP'
        if source_genes_type == 'dna':
            return 'DIAMOND-BLASTX'

    try:
        with open(file_path, 'r') as fh:
            for _ in range(200):
                line = fh.readline()
                if not line:
                    break
                l = line.lower()
                if 'blastx' in l:
                    return 'DIAMOND-BLASTX'
                if 'blastp' in l:
                    return 'DIAMOND-BLASTP'
                # sometimes the DIAMOND header includes the full command
                if 'diamond' in l and ('blastx' in l or 'blastp' in l):
                    if 'blastp' in l:
                        return 'DIAMOND-BLASTP'
                    if 'blastx' in l:
                        return 'DIAMOND-BLASTX'
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) >= 14:
                    try:
                        alignment_length = float(fields[3])
                        query_span = abs(int(fields[7]) - int(fields[6])) + 1
                        if alignment_length > 0:
                            ratio = query_span / alignment_length
                            if 2.5 <= ratio <= 3.5:
                                return 'DIAMOND-BLASTX'
                            if 0.75 <= ratio <= 1.25:
                                return 'DIAMOND-BLASTP'
                    except (TypeError, ValueError):
                        pass
    except Exception:
        pass
    return 'DIAMOND'


def load_original_run_parameters(input_path: str):
    input_dir = Path(input_path)
    root = input_dir.parent if input_dir.name == "raw_outputs" else input_dir
    params_path = root / "run_parameters.json"
    if not params_path.exists():
        return None, params_path
    try:
        with open(params_path, "r") as handle:
            return json.load(handle), params_path
    except Exception:
        return None, params_path


def parse_query_paths(query_fasta):
    """Return query sequence paths from a single or paired CLI value."""
    return [
        Path(part.strip())
        for part in str(query_fasta or '').split(',')
        if part.strip()
    ]


def apply_source_run_defaults(options, logger):
    """Fill unspecified recompute settings from the source run manifest."""
    original, params_path = load_original_run_parameters(options.input)
    original = original or {}
    if original:
        logger.info(f"Using source run metadata from: {params_path}")
    else:
        logger.warning(
            f"Source run metadata was not available at {params_path}; "
            "using documented recompute defaults."
        )

    defaults = {
        'query_min_coverage': 40.0,
        'query_min_identity': 80.0,
        'detection_min_coverage': 80.0,
        'detection_min_identity': 80.0,
        'detection_min_depth': 1,
        'detection_min_num_reads': 1,
        'evidence_corroborating_depth': 3,
        'evidence_exact_identity': 100.0,
        'evidence_candidate_depth': 3,
        'evidence_candidate_identity': 98.0,
        'evidence_max_internal_gap_bp': 15,
        'evidence_max_internal_gap_fraction': 0.02,
        'evidence_min_unique_reads': 2,
        'evidence_min_unique_fraction': 0.10,
        'evidence_ambiguity_fraction': 0.50,
        'evidence_score_tie': 1.0,
        'sequence_type': 'Single-FASTA',
        'genes_type': None,
        'evalue': None,
        'min_bitscore': None,
        'always_flag_genes': [],
        'hamronized_output': None,
        'hamronized_min_call': 'evidence',
        'sample_id': None,
    }
    for name, fallback in defaults.items():
        if getattr(options, name, None) is not None:
            continue
        if name == 'detection_min_depth':
            value = original.get(
                'detection_min_depth',
                original.get('detection_min_base_depth', fallback),
            )
        else:
            value = original.get(name, fallback)
        if value is None and fallback is not None:
            value = fallback
        if name == 'genes_type' and value == 'protein':
            value = 'aa'
        setattr(options, name, value)

    options.source_run_parameters = original
    return original


def validate_recompute_relaxation(options, logger) -> bool:
    """Reject threshold relaxation that cannot be recovered from filtered output."""
    if 'all' not in options.tools and not any(
            tool in options.tools for tool in ['blastn', 'blastx', 'blastp', 'diamond']):
        return True
    original, params_path = load_original_run_parameters(options.input)
    if original is None:
        logger.warning(
            f"Original run parameters were not available at {params_path}; "
            "cannot verify whether the requested thresholds are complete."
        )
        return True

    incomplete = []
    original_qcov = original.get("query_min_coverage")
    original_qid = original.get("query_min_identity")
    original_bitscore = original.get("min_bitscore")
    original_evalue = original.get("evalue")

    if original_qcov is not None and options.query_min_coverage < float(original_qcov):
        incomplete.append(
            f"query coverage {options.query_min_coverage}% < source {original_qcov}%"
        )
    if original_qid is not None and options.query_min_identity < float(original_qid):
        incomplete.append(
            f"query identity {options.query_min_identity}% < source {original_qid}%"
        )
    if original_bitscore is not None:
        if options.min_bitscore is None or options.min_bitscore < float(original_bitscore):
            incomplete.append(
                f"minimum bitscore {options.min_bitscore} is less restrictive than "
                f"source {original_bitscore}"
            )
    if original_evalue is not None:
        if options.evalue is None or options.evalue > float(original_evalue):
            incomplete.append(
                f"e-value {options.evalue} is less restrictive than source {original_evalue}"
            )

    if not incomplete:
        return True

    logger.error(
        "The requested recompute thresholds are more permissive than the source run: "
        + "; ".join(incomplete)
    )
    logger.error(
        "The stored BLAST/DIAMOND tables contain only hits retained by the source "
        "query-level filters, so discarded hits cannot be recovered."
    )
    if options.allow_incomplete_relaxation:
        logger.warning(
            "Proceeding because --allow-incomplete-relaxation was supplied. Results may "
            "underestimate coverage, depth, identity support, and detections."
        )
        return True
    logger.error(
        "Re-run the source search with permissive thresholds, or use "
        "--allow-incomplete-relaxation to explicitly accept incomplete results."
    )
    return False




def run(options, workflow, logger):

    databases_found = set()
    tools_found = defaultdict(set)

    logger.info("=" * 70)
    logger.info(f"Genefíor-Recompute {GENEFIOR_VERSION}")
    logger.info("=" * 70)
    logger.info(f"Input directory: {options.input}")
    logger.info(f"Output directory: {options.output}")
    logger.info(f"Detection thresholds:")
    logger.info(f"  Query min coverage: {options.query_min_coverage}%")
    logger.info(f"  Query min id: {options.query_min_identity}%")
    logger.info(f"  Gene min coverage: {options.detection_min_coverage}%")
    logger.info(f"  Min identity: {options.detection_min_identity}%")
    logger.info(f"  Legacy detection depth: {options.detection_min_depth}×")
    logger.info(f"  Qualified evidence depth: {options.evidence_corroborating_depth}×")
    logger.info(f"  Min num reads: {options.detection_min_num_reads}")
    logger.info(f"  Detection system: {options.detection_system}")
    logger.info("=" * 70)

    # Discover files
    discovered = discover_files(options, logger, databases_found, tools_found)
    if not discovered:
        return False
    found_files, databases_found, tools_found = discovered
    try:
        validate_database_versions(
            options.hamronized_output,
            found_files.keys(),
            options.database_versions,
        )
    except ValueError as exc:
        logger.error(str(exc))
        return False

    #if not discover_files():
    #    return False

    # Process each database
    for database, tools in sorted(found_files.items()):
        logger.info(f"\n### Processing {database.upper()} ###")
        workflow.databases = database
        #files = found_files[database]
        for tool_info in tools:
            tool = tool_info[0]
            file_path = tool_info[1]
            read_store = None
            wants_read_outputs = bool(options.report_read_names or workflow.report_fasta)
            if wants_read_outputs:
                fd, temp_name = tempfile.mkstemp(
                    prefix="genefior_recompute_reads_",
                    suffix=".sqlite",
                    dir=str(workflow.output_dir),
                )
                os.close(fd)
                os.unlink(temp_name)
                read_store = ReadMatchStore(Path(temp_name))

            try:
                if tool in ['blastn', 'blastx', 'blastp', 'diamond']:
                    # Map short tool keys to the tool_name strings expected by Workflow.parse_blast_results
                    if tool == 'blastn':
                        tool_name = 'BLASTn'
                    elif tool == 'blastx':
                        tool_name = 'BLASTx'
                    elif tool == 'blastp':
                        tool_name = 'BLASTp'
                    else:  # diamond
                        tool_name = detect_diamond_mode(
                            file_path,
                            getattr(options, 'source_run_parameters', None),
                        )

                    detected, gene_reads = workflow.parse_blast_results(
                        file_path, database, tool_name,
                        count_only=True,
                        read_store=read_store,
                    )
                    workflow.write_tool_stats(database, tool_name, gene_reads)
                    workflow.gene_stats.get(database, {}).pop(tool_name, None)
                elif tool in ['bowtie2', 'bwa', 'minimap2']:
                    tool_name = {
                        'bowtie2': 'Bowtie2',
                        'bwa': 'BWA',
                        'minimap2': 'Minimap2',
                    }[tool]
                    detected, gene_reads = workflow.parse_bam_results(
                        file_path, database, tool_name,
                        count_only=True,
                        read_store=read_store,
                    )
                    workflow.write_tool_stats(database, tool_name, gene_reads)
                    workflow.gene_stats.get(database, {}).pop(tool_name, None)
                else:
                    continue

                if read_store is not None:
                    read_store.flush()
                    evidence_genes = {
                        gene
                        for gene, calls in workflow.evidence_calls
                        .get(database, {}).items()
                        if tool_name in calls
                        and calls[tool_name].evidence_present
                    }
                    exact_genes = {
                        gene
                        for gene, calls in workflow.evidence_calls
                        .get(database, {}).items()
                        if tool_name in calls
                        and calls[tool_name].exact_allele_detected
                    }
                    candidate_genes = {
                        gene
                        for gene, calls in workflow.evidence_calls
                        .get(database, {}).items()
                        if tool_name in calls
                        and calls[tool_name].candidate_allele_detected
                    }
                    if workflow.report_fasta and tool in ['blastn', 'blastx', 'blastp', 'diamond']:
                        matched_records = 0
                        updated_rows = 0
                        query_paths = parse_query_paths(options.query_fasta)
                        for index, query_path in enumerate(query_paths):
                            suffix = (
                                f"/{index + 1}"
                                if len(query_paths) == 2
                                else None
                            )
                            matched, updated = read_store.populate_sequences(
                                query_path,
                                suffix=suffix,
                            )
                            matched_records += matched
                            updated_rows += updated
                        logger.info(
                            f"Matched {matched_records} query records to "
                            f"{updated_rows} stored {database}/{tool_name} read-gene matches."
                        )

                    if options.report_read_names:
                        read_name_mode = options.report_read_names
                        selected_genes = evidence_genes
                        if read_name_mode in ('evidence', 'evidence-all'):
                            selected_genes = evidence_genes
                            read_name_mode = (
                                'detected'
                                if read_name_mode == 'evidence'
                                else 'detected-all'
                            )
                        elif read_name_mode in ('exact', 'exact-all'):
                            selected_genes = exact_genes
                            read_name_mode = (
                                'detected'
                                if read_name_mode == 'exact'
                                else 'detected-all'
                            )
                        elif read_name_mode in ('candidate', 'candidate-all'):
                            selected_genes = candidate_genes
                            read_name_mode = (
                                'detected'
                                if read_name_mode == 'candidate'
                                else 'detected-all'
                            )
                        exported = read_store.export_names(
                            workflow.output_dir / "read_names",
                            read_name_mode,
                            detected_genes=selected_genes,
                        )
                        logger.info(
                            f"Exported {exported} read names for {database}/{tool_name}."
                        )

                    if workflow.report_fasta:
                        fasta_mode = workflow.report_fasta
                        selected_genes = evidence_genes
                        if fasta_mode in ('evidence', 'evidence-all'):
                            selected_genes = evidence_genes
                            fasta_mode = (
                                'detected'
                                if fasta_mode == 'evidence'
                                else 'detected-all'
                            )
                        elif fasta_mode in ('exact', 'exact-all'):
                            selected_genes = exact_genes
                            fasta_mode = (
                                'detected'
                                if fasta_mode == 'exact'
                                else 'detected-all'
                            )
                        elif fasta_mode in ('candidate', 'candidate-all'):
                            selected_genes = candidate_genes
                            fasta_mode = (
                                'detected'
                                if fasta_mode == 'candidate'
                                else 'detected-all'
                            )
                        exported, missing = read_store.export_fasta(
                            workflow.fasta_dir,
                            fasta_mode,
                            detected_genes=selected_genes,
                        )
                        logger.info(
                            f"Exported {exported} FASTA records for {database}/{tool_name}; "
                            f"{missing} selected records had no sequence available."
                        )
            finally:
                if read_store is not None:
                    read_store.close(delete=True)

        # Generate detection matrix
        workflow.generate_detection_matrix(database)

    export_from_options(
        options,
        options.output,
        "GeneFior-Recompute",
        GENEFIOR_VERSION.lstrip("v"),
        logger=logger,
        input_spec=(options.query_fasta or options.input),
    )

    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("RECOMPUTATION COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Recomputed statistics saved to: {workflow.stats_dir}")
    logger.info(f"Detection matrices saved to: {options.output}")
    logger.info("=" * 70)

    return True


def main():
    parser = argparse.ArgumentParser(
        description='GeneFíor ' + GENEFIOR_VERSION + ' - GeneFíor-Recompute: Recalculate detection statistics from existing sequence search outputs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Recompute with different thresholds
  GeneFior-recompute -i original_results/ -o recomputed_90_90/ \\
    --d-min-cov 90 --d-min-id 90

  # More stringent depth requirement
  GeneFior-recompute -i original_results/ -o high_depth/ \\
    --d-min-depth 5 --d-min-reads 10

        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input directory containing Genefíor results (with raw_outputs/ subdirectory)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for recomputed results')
    parser.add_argument('--tools', nargs='+',
                            choices=['blastn', 'blastx', 'blastp', 'diamond', 'bowtie2', 'bwa', 'minimap2', 'all'], #, 'hmmer_dna', 'hmmer_protein'],
                            default=['all'], #, 'hmmer_dna','hmmer_protein'],
                            help='Specify which tools to recompute - "all" will recompute for all detected tools (default: all)')

    query_threshold_group = parser.add_argument_group('Query threshold Parameters')
    query_threshold_group.add_argument('--q-min-cov', '--query-min-coverage', type=float, default=None,
                                      dest='query_min_coverage',
                                      help='Minimum coverage threshold in percent (default: source run, otherwise 40.0)')
    query_threshold_group.add_argument('--q-min-id', '--query-min-identity', type=float, default=None,
                                      dest='query_min_identity',
                                      help='Minimum identity threshold in percent (default: source run, otherwise 80.0)')

    gene_detection_group = parser.add_argument_group('Gene Detection Parameters')
    gene_detection_group.add_argument('--d-min-cov', '--detection-min-coverage', type=float, default=None,
                              dest='detection_min_coverage',
                              help='Minimum coverage threshold in percent (default: 80.0)')
    gene_detection_group.add_argument('--d-min-id', '--detection-min-identity', type=float, default=None,
                              dest='detection_min_identity',
                              help='Minimum identity threshold in percent (default: 80.0)')
    gene_detection_group.add_argument('--d-min-depth', '--detection-min-depth',
                              '--d-min-base-depth', '--detection-min-base-depth',
                              type=int, default=None,
                              dest='detection_min_depth',
                              help='Minimum per-base depth required across the '
                                   'configured detection coverage in legacy-relaxed '
                                   'mode. The old base-depth option names are '
                                   'accepted as deprecated aliases (default: source '
                                   'run, otherwise 3).')
    gene_detection_group.add_argument('--d-min-reads', '--detection-min-num-reads',
                              type=int, default=None,
                              dest='detection_min_num_reads',
                              help='Minimum number of reads required for detection (default: 1)')
    gene_detection_group.add_argument(
        '--detection-system', '--detection-mode',
        choices=['qualified', 'legacy-relaxed'],
        default='qualified',
        help='Detection interpretation: "qualified" uses evidence, family/mosaic, '
             'and exact-allele resolution (default); "legacy-relaxed" reproduces '
             'the original direct threshold-only detector.')
    gene_detection_group.add_argument('--evidence-corroborating-depth', type=int, default=None,
                              help='Per-base depth required across the configured detection coverage for qualified evidence calls (default: source run, otherwise 3)')
    gene_detection_group.add_argument('--evidence-exact-identity', type=float, default=None,
                              help='Deprecated identity setting retained for manifest compatibility; literal exact calls require 100%% identity and full candidate-depth coverage')
    gene_detection_group.add_argument('--evidence-candidate-depth', type=int, default=None,
                              help='Minimum per-base depth required across 100%% of a nucleotide allele candidate (default: 3)')
    gene_detection_group.add_argument('--evidence-candidate-identity', type=float, default=None,
                              help='Minimum mean identity for a nucleotide allele candidate (default: 98.0)')
    gene_detection_group.add_argument('--evidence-max-internal-gap-bp', type=int, default=None,
                              help='Largest unsupported internal gap allowed for an exact call (default: 15 bp)')
    gene_detection_group.add_argument('--evidence-max-internal-gap-fraction', type=float, default=None,
                              help='Gene-length fraction allowed as an internal gap (default: 0.02)')
    gene_detection_group.add_argument('--evidence-min-unique-reads', type=int, default=None,
                              help='Minimum competitively unique reads for an exact allele (default: 2)')
    gene_detection_group.add_argument('--evidence-min-unique-fraction', type=float, default=None,
                              help='Minimum fraction of passing reads uniquely supporting the allele (default: 0.10)')
    gene_detection_group.add_argument('--evidence-ambiguity-fraction', type=float, default=None,
                              help='Ambiguous-best read fraction that triggers an allele warning (default: 0.50)')
    gene_detection_group.add_argument('--evidence-score-tie', type=float, default=None,
                              help='Absolute bit-score difference treated as a competitive tie (default: 1.0)')
    gene_detection_group.add_argument(
        '--always-flag-genes', '--always-flag-gene-list',
        dest='always_flag_genes_file',
        default=None,
        help='CSV/TSV/text file of database gene IDs to always surface for '
             'review when any evidence exists. These genes are not counted as '
             'detected unless they pass the normal thresholds.')

    # Allow e-value and bitscore thresholds to be supplied for recompute
    parser.add_argument('--evalue', type=float, default=None,
                        help='Optional e-value threshold to apply when recomputing (passed to BLAST/DIAMOND)')
    parser.add_argument('--min-bitscore', type=float, default=None,
                        dest='min_bitscore',
                        help='Optional minimum bitscore to require for BLAST/DIAMOND hits when recomputing')
    parser.add_argument('--allow-incomplete-relaxation', action='store_true', default=False,
                        help='Allow thresholds that are more permissive than the source run, '
                             'acknowledging that previously discarded hits cannot be recovered.')

    # Allow user to indicate input is a Genes-FASTA and optionally its type
    parser.add_argument('--sequence-type', choices=['Single-FASTA', 'Paired-FASTQ', 'Genes-FASTA'], default=None,
                        help='Input sequence type (default: inherit from source run)')
    parser.add_argument('--genes-type', choices=['aa', 'dna'], default=None,
                        help='When using --sequence-type Genes-FASTA, optionally declare whether genes are amino-acid (aa) or nucleotide (dna)')

    # Output selection
    output_group = parser.add_argument_group('Output Parameters')
    output_group.add_argument('--report-fasta',
                            choices=[
                                'None', 'all', 'detected', 'detected-all',
                                'evidence', 'evidence-all',
                                'candidate', 'candidate-all',
                                'exact', 'exact-all',
                            ],
                            default=None,
                            dest='report_fasta',
                            help='Specify whether to output sequences that "mapped" to genes.'
                                 '"all" should only be used for deep investigation/debugging.'
                                 '"detected" will report the reads that passed detection thresholds for each detected gene.'
                                 '"detected-all" reports all mapped reads for each evidence gene.'
                                 '"evidence" and "evidence-all" are explicit aliases for the detected modes.'
                                 '"candidate" and "candidate-all" restrict output to candidate allele calls.'
                                 '"exact" and "exact-all" restrict output to exact nucleotide alleles.  (default: None)')
    output_group.add_argument('--report-read-names',
                            choices=[
                                'all', 'detected', 'detected-all',
                                'evidence', 'evidence-all',
                                'candidate', 'candidate-all',
                                'exact', 'exact-all',
                            ],
                            default=None,
                            dest='report_read_names',
                            help='Write deduplicated read-name files per database/tool/gene. '
                                 '"all" reports every mapped read, "detected" reports threshold-passing '
                                 'reads for detected genes, and "detected-all" reports every mapped read '
                                 'for evidence genes. "evidence" modes are aliases for detected modes; '
                                 '"candidate" modes restrict output to candidate allele calls; '
                                 '"exact" modes restrict output to exact nucleotide alleles. Uses disk-backed '
                                 'storage to avoid retaining names in memory.')
    output_group.add_argument('--query-fasta', dest='query_fasta',
                              help='Specify the original query FASTA/FASTQ file used for alignment (required for reporting '
                                   'mapped sequences for BLAST/DIAMOND).')

    misc_group = parser.add_argument_group('Miscellaneous Parameters')
    misc_group.add_argument('-v','--version', action='version',
                            version='GeneFíor ' + GENEFIOR_VERSION,
                            help='Show program version and exit')

    add_hamronization_arguments(parser, inherit_defaults=True)
    options = parser.parse_args()
    if options.report_fasta == 'None':
        options.report_fasta = None

    ## Setup logging
    start_time = datetime.now()
    from pathlib import Path
    log_file = Path(options.output) / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = logging.getLogger('GeneFíor-Recompute')
    logger.setLevel(logging.INFO)

    # Create stream handler for console output
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(stream_handler)

    log_file.parent.mkdir(parents=True, exist_ok=True)
    # Create file handler explicitly (use string path) and prepare formatter
    try:
        file_handler = logging.FileHandler(str(log_file), mode='a')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    except Exception as _err:
        file_handler = None  # fallback to console-only if file cannot be opened
    ################

    # Check input directory exists
    if not os.path.exists(options.input):
        logger.error(f"Error: Input directory '{options.input}' not found")
        sys.exit(1)

        # Tool selection
    apply_source_run_defaults(options, logger)
    try:
        inherited_versions = {
            str(key).casefold(): str(value)
            for key, value in (
                getattr(options, 'source_run_parameters', {})
                .get('database_versions', {}) or {}
            ).items()
        }
        supplied_versions = parse_database_versions(
            options.database_version_specs
        )
        options.database_versions = {
            **inherited_versions,
            **supplied_versions,
        }
        validate_hamronization_request(
            options.hamronized_output,
            options.hamronized_min_call,
            options.detection_system,
        )
    except ValueError as exc:
        logger.error(str(exc))
        sys.exit(2)
    try:
        inherited_flags = set(getattr(options, 'always_flag_genes', []) or [])
        supplied_flags = load_gene_id_file(
            getattr(options, 'always_flag_genes_file', None),
            logger,
        )
        options.always_flag_genes = inherited_flags | supplied_flags
    except Exception as exc:
        logger.error(f"Failed to load --always-flag-genes file: {exc}")
        sys.exit(1)
    if options.always_flag_genes:
        logger.info(
            f"Always-flag review genes: {len(options.always_flag_genes)}"
        )
    if not validate_recompute_relaxation(options, logger):
        sys.exit(2)
    # If user requested FASTA reporting for BLAST/DIAMOND outputs, require --query-fasta
    recomputes_tabular_searches = (
        'all' in options.tools
        or any(
            tool in options.tools
            for tool in ['blastx', 'blastn', 'blastp', 'diamond']
        )
    )
    if options.report_fasta is not None and recomputes_tabular_searches:
        if options.query_fasta is None:
            logger.error("Error: --query-fasta must be provided when --report-fasta is used with blast/diamond outputs")
            sys.exit(1)
        query_paths = parse_query_paths(options.query_fasta)
        if not query_paths or len(query_paths) > 2 or any(
                not path.is_file() for path in query_paths):
            logger.error(
                "Error: --query-fasta must identify one file, or two "
                "comma-separated paired files, and every file must exist: "
                f"{options.query_fasta}"
            )
            sys.exit(1)

    # Create output directory
    # options.output.mkdir(parents=True, exist_ok=True)
    # stats_dir = options.output / "recomputed_stats"
    # stats_dir.mkdir(exist_ok=True)

    workflow = Workflow(
        input_fasta=options.query_fasta,
        input_fastq=None,
        output_dir=options.output,
        databases={},  # to be set based on user input or discovered files
        threads=options.threads if hasattr(options, 'threads') else 4,
        tool_sensitivity_params={},  # to be set if user provides
        query_min_coverage=options.query_min_coverage,
        query_min_identity=options.query_min_identity,
        detection_min_coverage=options.detection_min_coverage,
        detection_min_identity=options.detection_min_identity,
        detection_min_depth=options.detection_min_depth,
        detection_min_num_reads=options.detection_min_num_reads,
        detection_system=options.detection_system,
        evalue=options.evalue,
        min_bitscore=options.min_bitscore,
        run_dna=True,
        run_protein=True,
        sequence_type=options.sequence_type,
        report_fasta=options.report_fasta,
        no_cleanup=False,
        verbose=False,
        logger = logger,
        genes_type=(
            'protein' if options.genes_type == 'aa'
            else ('dna' if options.genes_type == 'dna' else None)
        ),
        evidence_corroborating_depth=options.evidence_corroborating_depth,
        evidence_exact_identity=options.evidence_exact_identity,
        evidence_candidate_depth=options.evidence_candidate_depth,
        evidence_candidate_identity=options.evidence_candidate_identity,
        evidence_max_internal_gap_bp=options.evidence_max_internal_gap_bp,
        evidence_max_internal_gap_fraction=options.evidence_max_internal_gap_fraction,
        evidence_min_unique_reads=options.evidence_min_unique_reads,
        evidence_min_unique_fraction=options.evidence_min_unique_fraction,
        evidence_ambiguity_fraction=options.evidence_ambiguity_fraction,
        evidence_score_tie=options.evidence_score_tie,
        always_flag_genes=getattr(options, 'always_flag_genes', set()),
    )



    # write run parameters manifest to recompute output directory
    try:
        outdir = Path(options.output)
        outdir.mkdir(parents=True, exist_ok=True)
        params = {
            'query_min_coverage': options.query_min_coverage,
            'query_min_identity': options.query_min_identity,
            'detection_min_coverage': options.detection_min_coverage,
            'detection_min_identity': options.detection_min_identity,
            'detection_min_depth': options.detection_min_depth,
            'detection_min_num_reads': options.detection_min_num_reads,
            'detection_system': options.detection_system,
            'evidence_corroborating_depth': options.evidence_corroborating_depth,
            'evidence_exact_identity': options.evidence_exact_identity,
            'evidence_candidate_depth': options.evidence_candidate_depth,
            'evidence_candidate_identity': options.evidence_candidate_identity,
            'evidence_max_internal_gap_bp': options.evidence_max_internal_gap_bp,
            'evidence_max_internal_gap_fraction': options.evidence_max_internal_gap_fraction,
            'evidence_min_unique_reads': options.evidence_min_unique_reads,
            'evidence_min_unique_fraction': options.evidence_min_unique_fraction,
            'evidence_ambiguity_fraction': options.evidence_ambiguity_fraction,
            'evidence_score_tie': options.evidence_score_tie,
            'tools': options.tools,
            'report_fasta': options.report_fasta,
            'report_read_names': options.report_read_names,
            'allow_incomplete_relaxation': options.allow_incomplete_relaxation,
            'always_flag_genes_file': options.always_flag_genes_file,
            'always_flag_genes': sorted(getattr(options, 'always_flag_genes', set())),
            'hamronized_output': getattr(options, 'hamronized_output', None),
            'hamronized_min_call': getattr(options, 'hamronized_min_call', 'evidence'),
            'sample_id': (
                getattr(options, 'sample_id', None)
                or infer_sample_id(options.query_fasta or options.input, options.output)
            ),
            'database_versions': dict(getattr(options, 'database_versions', {})),
            'sequence_type': options.sequence_type,
            'genes_type': (
                'protein' if options.genes_type == 'aa'
                else options.genes_type
            ),
            'evalue': options.evalue,
            'min_bitscore': options.min_bitscore,
        }
        source_input = Path(options.input).resolve()
        source_raw_outputs = (
            source_input
            if source_input.name == "raw_outputs"
            else source_input / "raw_outputs"
        )
        params['source_input'] = str(source_input)
        if source_raw_outputs.is_dir():
            params['source_raw_outputs'] = str(source_raw_outputs)
        with open(outdir / 'run_parameters.json', 'w') as pf:
            json.dump(params, pf, indent=2)
    except Exception:
        logger.debug('Failed to write run parameters to recompute output: %s', options.output)

    success = run(options, workflow, logger)

    #shutil.rmtree(workflow.raw_dir)

    end_time = datetime.now()
    elapsed = end_time - start_time
    if success:
        logger.info(f"GeneFíor-Recompute completed in {elapsed}.")
    else:
        logger.error(f"GeneFíor-Recompute failed after {elapsed}.")

    for handler in logger.handlers:
        handler.flush()
    for handler in logger.handlers:
        handler.close()
    logging.shutdown()
    if not success:
        sys.exit(1)



if __name__ == "__main__":
    main()
