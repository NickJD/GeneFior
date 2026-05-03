import argparse
import sys
import os
from pathlib import Path
import logging
import json
from datetime import datetime
from collections import defaultdict

try:
    from .gene_stats import GeneStats
    from .constants import GENEFIOR_VERSION
    from .workflow import Workflow
except (ModuleNotFoundError, ImportError):
    from gene_stats import GeneStats
    from constants import GENEFIOR_VERSION
    from workflow import Workflow



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


def detect_diamond_mode(file_path: Path) -> str:
    """Try to infer whether a DIAMOND output was produced using blastx or blastp
    by scanning the header/comment lines of the output file. Returns one of
    'DIAMOND-BLASTX', 'DIAMOND-BLASTP' or 'DIAMOND' (unknown)."""
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
    except Exception:
        pass
    return 'DIAMOND'




def run(options, workflow, logger):

    # Storage structures
    gene_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(GeneStats)))
    detections = defaultdict(lambda: defaultdict(lambda: defaultdict(bool)))
    databases_found = set()
    tools_found = defaultdict(set)

    logger.info("=" * 70)
    logger.info(f"Genefíor-Recompute {GENEFIOR_VERSION}")
    logger.info("=" * 70)
    logger.info(f"Input directory: {options.input}")
    logger.info(f"Output directory: {options.output}")
    logger.info(f"Detection thresholds:")
    logger.info(f"  Query min coverage: {options.query_min_coverage}%")
    logger.info(f"  Query min id: {options.query_min_id}%")
    logger.info(f"  Gene min coverage: {options.detection_min_coverage}%")
    logger.info(f"  Min identity: {options.detection_min_identity}%")
    logger.info(f"  Min base depth: {options.detection_min_base_depth}×")
    logger.info(f"  Min num reads: {options.detection_min_num_reads}")
    logger.info("=" * 70)

    # Discover files
    found_files, databases_found, tools_found = discover_files(options, logger, databases_found, tools_found)

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
            if tool in ['blastn', 'blastx', 'blastp', 'diamond']:
                # Map short tool keys to the tool_name strings expected by Workflow.parse_blast_results
                if tool == 'blastn':
                    tool_name = 'BLASTn'
                elif tool == 'blastx':
                    tool_name = 'BLASTx'
                elif tool == 'blastp':
                    tool_name = 'BLASTp'
                else:  # diamond
                    tool_name = detect_diamond_mode(file_path)

                detected, gene_reads = workflow.parse_blast_results(file_path, database, tool_name)
                workflow.write_tool_stats(database, tool_name, gene_reads)
            elif tool in ['bowtie2', 'bwa', 'minimap2']:
                #tool_name = {'bowtie2': 'Bowtie2', 'bwa': 'BWA', 'minimap2': 'Minimap2'}[tool]
                detected, gene_reads = workflow.parse_bam_results(file_path, database, tool)
                workflow.write_tool_stats(database, tool, gene_reads)

        # Generate detection matrix
        workflow.generate_detection_matrix(database)

    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("RECOMPUTATION COMPLETE")
    logger.info("=" * 70)
    stats_dir = Path(options.output) / "recomputed_stats"
    logger.info(f"Recomputed statistics saved to: {stats_dir}")
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
    --d-min-base-depth 5.0 --d-min-reads 10

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
    query_threshold_group.add_argument('--q-min-cov', '--query-min-coverage', type=float, default=40.0,
                                      dest='query_min_coverage',
                                      help='Minimum coverage threshold in percent (default: 40.0)')
    query_threshold_group.add_argument('--q-min-id', '--query-min-identity', type=float, default=80.0,
                                      dest='query_min_identity',
                                      help='Minimum identity threshold in percent (HSP for blast/diamond) (default: 80.0)')

    gene_detection_group = parser.add_argument_group('Gene Detection Parameters')
    gene_detection_group.add_argument('--d-min-cov', '--detection-min-coverage', type=float, default=80.0,
                              dest='detection_min_coverage',
                              help='Minimum coverage threshold in percent (default: 80.0)')
    gene_detection_group.add_argument('--d-min-id', '--detection-min-identity', type=float, default=80.0,
                              dest='detection_min_identity',
                              help='Minimum identity threshold in percent (default: 80.0)')
    gene_detection_group.add_argument('--d-min-base-depth', '--detection-min-base-depth',
                              type=float, default=1.0,
                              dest='detection_min_base_depth',
                              help='Minimum average base depth for detection '
                                   '- calculated against regions of the detected gene with at least one read hit (default: 1.0)')
    gene_detection_group.add_argument('--d-min-reads', '--detection-min-num-reads',
                              type=int, default=1,
                              dest='detection_min_num_reads',
                              help='Minimum number of reads required for detection (default: 1)')

    # Allow e-value and bitscore thresholds to be supplied for recompute
    parser.add_argument('--evalue', type=float, default=None,
                        help='Optional e-value threshold to apply when recomputing (passed to BLAST/DIAMOND)')
    parser.add_argument('--min-bitscore', type=float, default=None,
                        dest='min_bitscore',
                        help='Optional minimum bitscore to require for BLAST/DIAMOND hits when recomputing')

    # Allow user to indicate input is a Genes-FASTA and optionally its type
    parser.add_argument('--sequence-type', choices=['Single-FASTA', 'Paired-FASTQ', 'Genes-FASTA'], default='Single-FASTA',
                        help='Indicate the type of input sequences. Set to "Genes-FASTA" to treat inputs as full-length gene FASTA(s)')
    parser.add_argument('--genes-type', choices=['aa', 'dna'], default=None,
                        help='When using --sequence-type Genes-FASTA, optionally declare whether genes are amino-acid (aa) or nucleotide (dna)')

    # Output selection
    output_group = parser.add_argument_group('Output Parameterts')
    output_group.add_argument('--report-fasta',
                            choices=['None', 'all', 'detected', 'detected-all'],
                            default=None,
                            dest='report_fasta',
                            help='Specify whether to output sequences that "mapped" to genes.'
                                 '"all" should only be used for deep investigation/debugging.'
                                 '"detected" will report the reads that passed detection thresholds for each detected gene.'
                                 '"detected-all" will report all reads for each detected gene.  (default: None)')
    output_group.add_argument('--query-fasta', dest='query_fasta',
                              help='Specify the original query FASTA/FASTQ file used for alignment (required for reporting '
                                   'mapped sequences for BLAST/DIAMOND).')

    misc_group = parser.add_argument_group('Miscellaneous Parameters')
    misc_group.add_argument('-v','--version', action='version',
                            version='GeneFíor ' + GENEFIOR_VERSION,
                            help='Show program version and exit')

    options = parser.parse_args()

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
    if options.tools == ['all']:
        # Expand "all" differently depending on sequence type. For recompute workflows
        # that operate on gene FASTA inputs we should avoid including read-mappers.
        if getattr(options, 'sequence_type', None) == 'Genes-FASTA':
            genes_type = getattr(options, 'genes_type', None)
            if genes_type == 'aa':
                options.tools = ['blastp', 'diamond']
            elif genes_type == 'dna':
                options.tools = ['blastn', 'diamond']
            else:
                # Unknown genes_type: include nucleotide and translated searches plus DIAMOND
                options.tools = ['blastn', 'blastx', 'diamond']
        else:
            # Non-Genes inputs: include full toolset
            options.tools = ['blastn', 'blastx', 'diamond', 'bowtie2', 'bwa', 'minimap2']  # , 'hmmer_dna','hmmer_protein']
    # If user requested FASTA reporting for BLAST/DIAMOND outputs, require --query-fasta
    if options.report_fasta is not None and any(tool in options.tools for tool in ['blastx', 'blastn', 'diamond']):
        if options.query_fasta is None:
            logger.error("Error: --query-fasta must be provided when --report-fasta is used with blast/diamond outputs")
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
        detection_min_base_depth=options.detection_min_base_depth,
        detection_min_num_reads=options.detection_min_num_reads,
        evalue=options.evalue,
        min_bitscore=options.min_bitscore,
        run_dna=True,
        run_protein=True,
        sequence_type=options.sequence_type,
        report_fasta=options.report_fasta,
        no_cleanup=False,
        verbose=False,
        logger = logger,
        genes_type=('aa' if options.genes_type == 'aa' else ('dna' if options.genes_type == 'dna' else None))
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
            'detection_min_base_depth': options.detection_min_base_depth,
            'detection_min_num_reads': options.detection_min_num_reads,
            'tools': options.tools
        }
        with open(outdir / 'run_parameters.json', 'w') as pf:
            json.dump(params, pf, indent=2)
    except Exception:
        logger.debug('Failed to write run parameters to recompute output: %s', options.output)

    run(options, workflow, logger)

    #shutil.rmtree(workflow.raw_dir)

    for handler in logger.handlers:
        handler.flush()
    for handler in logger.handlers:
        handler.close()
    logging.shutdown()

    end_time = datetime.now()
    elapsed = end_time - start_time
    logger.info(f"GeneFíor-Recompute completed in {elapsed}.")



if __name__ == "__main__":
    main()
