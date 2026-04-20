import argparse
import sys, os
import logging
import json
from pathlib import Path
from datetime import datetime

try:
    from .constants import *
    from .databases import RESFINDER_DATABASES, CARD_DATABASES, NCBI_DATABASES, gather_databases
    from .workflow import Workflow
    from .gene_stats import GeneStats
    from .utils import (
        handle_all_input_files,
        cleanup_all_temp_files,
        discover_samples_from_input_dir,
        discover_samples_from_subdirs,
        add_file_handler_for_path,
        remove_file_handler,
        get_max_fasta_chunk_bytes,
    )

except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    from databases import RESFINDER_DATABASES, CARD_DATABASES, NCBI_DATABASES, gather_databases
    from workflow import Workflow
    from gene_stats import GeneStats
    from utils import (
        handle_all_input_files,
        cleanup_all_temp_files,
        discover_samples_from_input_dir,
        discover_samples_from_subdirs,
        add_file_handler_for_path,
        remove_file_handler,
        get_max_fasta_chunk_bytes,
    )


def main():
    parser = argparse.ArgumentParser(
        description='GeneFíor ' + GENEFIOR_VERSION + ' - AMRfíor - The Multi-Tool AMR Gene Detection Module.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default tools (runs DNA & protein tools)
  AMRfíor -i reads.fasta -st Single-FASTA -o results

  # Select specific tools and output detected FASTA sequences
  AMRfíor -i reads.fasta -st Single-FASTA -o results \
    --tools diamond bowtie2 \
    --report_fasta detected

  # Custom thresholds, paired-fastq input, threads and dna-only mode
  AMRfíor -i reads_R1.fastq,reads_R2.fastq -st Paired-FASTQ -o results/ \
    -t 16 --d-min-cov 90 --d-min-id 85 \
    --dna-only
        """
    )

    # Required arguments
    required_group = parser.add_argument_group('Required selection')
    required_group.add_argument('-i', '--input', required=False,
                        help='Input FASTA/FASTAQ file(s) with sequences to analyse - Separate FASTQ R1 and R2 '
                             'with a comma for Paired-FASTQ or single file path for Single-FASTA - .gz files accepted')
    # Optional directory-based inputs
    required_group.add_argument('--input-dir', dest='input_dir', required=False, default=None,
                        help='Path to a directory containing multiple sample FASTA files (for Single-FASTA) or multiple paired FASTQ files (for Paired-FASTQ).')
    required_group.add_argument('--input-subdirs', dest='input_subdirs', required=False, default=None,
                        help='Path to a directory where each subdirectory contains one sample (a FASTA or paired FASTQ files).')
    required_group.add_argument('-st', '--sequence-type', required=True,
                        choices=['Single-FASTA', 'Paired-FASTQ', 'Genes-FASTA'],
                        help='Specify the input Sequence Type: Single-FASTA, Paired-FASTQ (R1+R2) or Genes-FASTA. '
                             'When Genes-FASTA is selected the pipeline treats the input as full-length gene FASTA(s) '
                             '(DNA or protein) and will skip read-mappers (bowtie2, bwa, minimap2).')
    # When using Genes-FASTA, require the user to state whether sequences are DNA or AA
    required_group.add_argument('--genes-type', choices=['dna', 'aa'], dest='genes_type', default=None,
                                   help='(Required when -st Genes-FASTA) Specify whether provided genes FASTA contains DNA (dna) or amino-acid sequences (aa)')
    required_group.add_argument('-o', '--output', required=True,
                        help='Output directory for results')

    # Output selection
    output_group = parser.add_argument_group('Output selection')
    output_group.add_argument('--report-fasta',
                            choices=['all', 'detected', 'detected-all'], #, 'hmmer_dna', 'hmmer_protein'],
                            default=None,
                            dest='report_fasta',
                            help='Specify whether to output sequences that "mapped" to genes.'
                                 '"all" should only be used for deep investigation/debugging.'
                                 '"detected" will report the reads that passed detection thresholds for each detected gene.'
                                 '"detected-all" will report all reads for each detected gene.  (default: None)')

    # Tool selection
    tool_group = parser.add_argument_group('Tool selection')
    tool_group.add_argument('--tools', nargs='+',
                            choices=['blastn', 'blastx', 'blastp', 'diamond', 'bowtie2', 'bwa', 'minimap2', 'all'], #, 'hmmer_dna', 'hmmer_protein'],
                            default=['diamond', 'bowtie2', 'bwa', 'minimap2'], #, 'hmmer_dna','hmmer_protein'],
                            help='Specify which tools to run - "all" will run all tools'
                                 ' (default: all except blastx/p as they can be slow!!)')

    # Database selection
    db_group = parser.add_argument_group('Database selection')
    db_group.add_argument('--databases', nargs='+',
                          choices=['resfinder', 'card', 'ncbi', 'user-provided'],
                          default=['resfinder', 'card'],
                          help='Specify which AMR gene databases to use (default: resfinder and card) -If "user-provided" is selected, '
                               'please ensure the path contains the appropriate databases set up as per the documentation '
                               'and specify the path with --user-db-path.')
    db_group.add_argument('--user-db-path', type=str,
                          help='Path to the directory containing user-provided databases (required if --databases includes "user-provided")')

    query_threshold_group = parser.add_argument_group('Query threshold Parameters')
    query_threshold_group.add_argument('--q-min-cov', '--query-min-coverage', type=float, default=40.0,
                                      dest='query_min_coverage',
                                      help='Minimum coverage threshold in percent (HSP for blastx/n) (default: 40.0)')
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

    # gene_detection_group.add_argument( '--max_target_seqs', dest='max_target_seqs', type=int, default=100,
    #                           help='Maximum number of "hits" to return per query sequence (default: 100)')


    # Mode selection
    mode_group = parser.add_argument_group('Mode Selection')
    mode_group.add_argument('--dna-only', action='store_true',
                            help='Run only DNA-based tools')
    mode_group.add_argument('--protein-only', action='store_true',
                            help='Run only protein-based tools')
    mode_group.add_argument('--sensitivity', type=str, default='default',
                            choices=['default', 'very-fast', 'fast', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'],
                            help='Preset sensitivity levels. "default" leaves tools unchanged. Other presets map to common '
                                 'DIAMOND and Bowtie2 presets (e.g. "very-sensitive" -> DIAMOND --very-sensitive / Bowtie2 --very-sensitive-local). '
                                 'Use --tool-param to override specific tool parameters if needed.')

    # Tool-specific parameters
    tool_params_group = parser.add_argument_group('Tool-Specific Parameters')
    tool_params_group.add_argument('--minimap2-preset', default='sr',
                                   choices=['sr', 'map-ont', 'map-pb', 'map-hifi'],
                                   help='Minimap2 preset: sr=short reads, map-ont=Oxford Nanopore, '
                                        'map-pb=PacBio, map-hifi=PacBio HiFi (default: sr)')
    tool_params_group.add_argument('--blastx-task', default='blastx-fast',
                                   choices=['blastx', 'blastx-fast'], dest='blastx_task',
                                   help='Run blastx with task blastx-fast (default: blastx-fast)')
    tool_params_group.add_argument('--tool-param', action='append', dest='tool_param', default=[],
                                   help="Specify extra parameters for a tool as TOOL=\"args\". Can be provided multiple times. Example: --tool-param 'diamond=--more-sensitive -e 1e-5'")
    # Tool-wide evalue: belongs with tool-specific parameters
    tool_params_group.add_argument('-e', '--evalue', type=float, default=None,
                                   help='E-value threshold to pass to BLASTn/x and DIAMOND (if not set, the tools\' default is used)')
    tool_params_group.add_argument('--min-bitscore', type=float, default=None,
                                   help='Minimum bitscore threshold to apply to BLAST/DIAMOND hits (if not set, no bitscore filtering is applied)')

    # Runtime parameters
    runtime_group = parser.add_argument_group('Runtime Parameters')
    runtime_group.add_argument('-t', '--threads', type=int, default=4,
                              help='Number of threads to use (default: 4)')
    runtime_group.add_argument('--chunk-jobs', type=int, default=None,
                               help='Number of concurrent BLAST chunk jobs to run when chunking is active. If unset the pipeline auto-derives concurrency from total threads or defaults to 1')
    runtime_group.add_argument('--chunk-threads-per-job', type=int, default=None,
                               help='If set, reserve this many threads per chunk job; otherwise total threads are divided evenly across concurrent chunk jobs')
    runtime_group.add_argument('--preserve-chunks', action='store_true', default=False,
                               help='Keep chunk files and per-chunk outputs after concatenation (useful for debugging)')
    runtime_group.add_argument('--max-fasta-chunk-mb', type=float, default=200.0,
                               help='Max FASTA chunk size in MiB (default: 200.0). Inputs larger than this will be chunked for per-chunk BLAST runs')
    runtime_group.add_argument('-tmp', '--temp-directory', type=str, default=None,
                               help='Path to temporary to place input FASTA/Q file(s) for faster IO during BLAST - '
                                    'Path will also be used for all temporary files (default: output directory)')
    runtime_group.add_argument('--force-modify-fastq', action='store_true',
                               help='Force addition of /1 and /2 suffixes to paired FASTQ read IDs even if they appear unique')
    runtime_group.add_argument('--no_cleanup',  action='store_true',)
    runtime_group.add_argument( '--verbose', action='store_true',)


    misc_group = parser.add_argument_group('Miscellaneous Parameters')
    misc_group.add_argument('-v','--version', action='version',
                            version='AMRfíor ' + GENEFIOR_VERSION,
                            help='Show program version and exit')

    options = parser.parse_args()

    # Ensure a temporary directory is defined. If the user did not supply
    # --temp-directory, default to the output directory (behaviour documented
    # in menu). This makes options.temp_directory available to
    # downstream code (e.g. utils.handle_all_input_files) without needing
    # checks everywhere.
    if getattr(options, 'temp_directory', None) is None:
        options.temp_directory = options.output

    # Enforce genes_type when sequence_type is Genes-FASTA
    if getattr(options, 'sequence_type', None) == 'Genes-FASTA' and getattr(options, 'genes_type', None) is None:
        print("Error: --genes-type (dna|aa) is required when -st/--sequence-type is Genes-FASTA", file=sys.stderr)
        sys.exit(1)

    # Ensure at least one input specification is provided: -i or --input-dir or --input-subdirs
    if (getattr(options, 'input', None) in (None, '') and
            getattr(options, 'input_dir', None) in (None, '') and
            getattr(options, 'input_subdirs', None) in (None, '')):
        print("Error: Please provide input using -i/--input or --input-dir or --input-subdirs", file=sys.stderr)
        sys.exit(1)

    ## Setup logging
    start_time = datetime.now()
    from pathlib import Path
    log_file = Path(options.output) / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = logging.getLogger('AMRfíor')
    # Only configure handlers if none exist to avoid duplicate logs when reimported
    if not getattr(logger, 'handlers', None):
        logger.setLevel(logging.INFO)
        # Create stream handler for console output
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(logging.INFO)
        stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(stream_handler)
        logger.propagate = False

    log_file.parent.mkdir(parents=True, exist_ok=True)
    # Create a single master file handler for all runs (including batch samples)
    # so all sample logs go into the root output directory as a single pipeline log.
    file_handler = add_file_handler_for_path(logger, log_file)
    ################

    # Normalize report_fasta: allow None or explicit string values
    if options.report_fasta in (None, 'None'):
        options.report_fasta = None


    try:
       # Tool selection
        if options.tools == ['all']:
            # For Genes-FASTA inputs, limit "all" to tools appropriate for gene FASTA (no read-mappers).
            if getattr(options, 'sequence_type', None) == 'Genes-FASTA':
                genes_type = getattr(options, 'genes_type', None)
                if genes_type == 'aa':
                    options.tools = ['blastp', 'diamond']
                elif genes_type == 'dna':
                    options.tools = ['blastn', 'diamond']
                else:
                    options.tools = ['blastn', 'blastx', 'diamond']
            else:
                options.tools = ['blastn', 'blastx', 'blastp', 'diamond', 'bowtie2', 'bwa', 'minimap2']

        # Initialise database variables to avoid UnboundLocalError
        resfinder_dbs = None
        card_dbs = None
        ncbi_dbs = None
        user_dbs = None
        # Load database paths from databases
        #logger.info("Database selected: " + ', '.join(options.databases))
        if 'resfinder' in options.databases:
            resfinder_dbs = {}
            for tool in options.tools or []:
                # Direct mapping if present
                val = RESFINDER_DATABASES.get(tool)
                # Allow 'blastp' to fall back to packaged 'blastx' DBs
                if val is None and tool == 'blastp':
                    val = RESFINDER_DATABASES.get('blastx')
                if val:
                    resfinder_dbs[tool] = val
        if 'card' in options.databases:
            card_dbs = {}
            for tool in options.tools or []:
                val = CARD_DATABASES.get(tool)
                if val is None and tool == 'blastp':
                    val = CARD_DATABASES.get('blastx')
                if val:
                    card_dbs[tool] = val
        if 'ncbi' in options.databases:
            ncbi_dbs = {}
            for tool in options.tools or []:
                val = NCBI_DATABASES.get(tool)
                if val is None and tool == 'blastp':
                    val = NCBI_DATABASES.get('blastx')
                if val:
                    ncbi_dbs[tool] = val
        if 'user-provided' in options.databases:
            if not hasattr(options, 'user_db_path') or not os.path.isdir(options.user_db_path):
                print("Error: Please provide a valid directory path for user-provided databases using --user-db-path",
                      file=sys.stderr)
                sys.exit(1)
            user_dbs = gather_databases(options.user_db_path, options.tools)

        databases = {
            'resfinder': resfinder_dbs,
            'card': card_dbs,
            'ncbi': ncbi_dbs,
            'user-provided-db': user_dbs
        }
        # Filter out None values
        databases = {key: value for key, value in databases.items() if value}
        if not databases:
            logger.error("Error: At least one database must be specified in databases.py or provided by the user")
            sys.exit(1)

        # Preflight: log which database tools (blastn/blastx/blastp/diamond) have database files available
        try:
            logger.info("Database tool availability (per selected database):")
            for db_name, db_map in databases.items():
                try:
                    # db_map is expected to be a dict mapping tool -> db_path (or similar truthy value)
                    # Some packaged DB constants expose 'blastx' for protein DBs but not 'blastp'.
                    # If the user explicitly requested 'blastp' we should expose a fallback mapping
                    # from an available 'blastx' entry so BLASTp can be used for protein Genes-FASTA.
                    if isinstance(db_map, dict):
                        if 'blastp' not in db_map and 'blastx' in db_map:
                            db_map['blastp'] = db_map.get('blastx')

                    bn = 'Yes' if db_map.get('blastn') else 'No'
                    bx = 'Yes' if db_map.get('blastx') else 'No'
                    bp = 'Yes' if db_map.get('blastp') else 'No'
                    dm = 'Yes' if db_map.get('diamond') else 'No'
                    logger.info(f"  {db_name}: BLASTn: {bn}; BLASTx: {bx}; BLASTp: {bp}; DIAMOND: {dm}")
                except Exception:
                    logger.debug(f"Failed to inspect database mapping for {db_name}")
        except Exception:
            logger.debug("Failed to emit database tool availability preflight logs")

        # Determine run modes
        run_dna = True
        run_protein = True

        if options.dna_only:
            run_protein = False
        if options.protein_only:
            run_dna = False

        if not run_dna and not run_protein:
            print("Error: Cannot disable both DNA and protein modes", file=sys.stderr)
            sys.exit(1)



        # Tool sensitivity
        tool_sensitivity_params = {}

        # Map user-friendly sensitivity presets to tool-specific flags. These are conservative defaults
        # and can be overridden per-tool using --tool-param TOOL="args".
        SENSITIVITY_MAP = {
            'very-fast': {'diamond': '--faster', 'bowtie2': '--very-fast'},
            'fast': {'diamond': '--fast', 'bowtie2': '--fast'},
            'sensitive': {'diamond': '--sensitive', 'bowtie2': '--sensitive'},
            'more-sensitive': {'diamond': '--more-sensitive', 'bowtie2': '--sensitive'},
            'very-sensitive': {'diamond': '--very-sensitive', 'bowtie2': '--very-sensitive-local'},
            'ultra-sensitive': {'diamond': '--ultra-sensitive', 'bowtie2': '--very-sensitive-local'},
            'default': {}
        }

        sel = getattr(options, 'sensitivity', 'default') if hasattr(options, 'sensitivity') else 'default'
        if sel in SENSITIVITY_MAP and SENSITIVITY_MAP[sel]:
            for t, flag in SENSITIVITY_MAP[sel].items():
                tool_sensitivity_params[t] = {'sensitivity': flag}
        # Parse extra tool parameters passed via CLI
        tool_extra_params = {}
        try:
            for tp in getattr(options, 'tool_param', []) or []:
                if not tp:
                    continue
                if '=' in tp:
                    key, val = tp.split('=', 1)
                    key = key.strip()
                    val = val.strip()
                    if key:
                        if key in tool_extra_params:
                            tool_extra_params[key] = tool_extra_params[key] + ' ' + val
                        else:
                            tool_extra_params[key] = val
                else:
                    logger.warning(f"Ignored malformed --tool-param value: {tp} (expected KEY=VAL)")
        except Exception:
            logger.debug("Failed to parse --tool-param values")

        ################################ Opening Text File for Logging ###############################
        logger.info("=" * 70)
        logger.info("GeneFíor - The Gene Detection toolkit: " + GENEFIOR_VERSION + " with the AMRfíor Module")
        logger.info("Started at: " + start_time.strftime("%Y-%m-%d %H:%M:%S"))
        logger.info("=" * 70)

        logger.info(f"Threads: {options.threads}")
        logger.info(f"Database(s) chosen:")
        logger.info(f"  ResFinder: {'Yes' if 'resfinder' in databases or 'all' in databases else 'No'}")
        logger.info(f"  CARD: {'Yes' if 'card' in databases or 'all' in databases else 'No'}")
        logger.info(f"  NCBI: {'Yes' if 'ncbi' in databases or 'all' in databases else 'No'}")
        logger.info(f"  User-Provided: {'Yes' if 'user-provided-db' in databases else 'No'}")

        # If user selected Genes-FASTA input, clarify that read-mappers will be skipped
        if getattr(options, 'sequence_type', None) == 'Genes-FASTA':
            logger.info("Sequence type: Genes-FASTA (full-length gene FASTA). Read-mappers (bowtie2, bwa, minimap2) will be skipped by the workflow.")
            genes_type_mapped = ('protein' if getattr(options, 'genes_type', None) == 'aa' else getattr(options, 'genes_type', None))
            if genes_type_mapped == 'protein':
                logger.info("  Genes declared as: AMINO-ACID (aa) — protein searches will be preferred (BLASTp / DIAMOND blastp when available).")
            elif genes_type_mapped == 'dna':
                logger.info("  Genes declared as: DNA — will run BLASTn and translated protein searches (BLASTx / DIAMOND blastx).")
            else:
                logger.info("  Genes type: not specified — workflow will attempt auto-detection to choose BLAST/DIAMOND modes.")

            logger.info(" Tool(s) chosen (note: read-mappers will be skipped for Genes-FASTA):")
            logger.info(f"  BLASTn: {'Yes' if 'blastn' in options.tools else 'No'}")
            logger.info(f"  BLASTx: {'Yes' if 'blastx' in options.tools else 'No'}")
            logger.info(f"  BLASTp: {'Yes' if 'blastp' in options.tools else 'No'}")
            logger.info(f"  DIAMOND: {'Yes' if 'diamond' in options.tools else 'No'}")
            logger.info(f"  Bowtie2: {'Present (will be SKIPPED)' if 'bowtie2' in options.tools else 'No'}")
            logger.info(f"  BWA: {'Present (will be SKIPPED)' if 'bwa' in options.tools else 'No'}")
            logger.info(f"  Minimap2: {'Present (will be SKIPPED)' if 'minimap2' in options.tools else 'No'}")
        else:
            logger.info(f" Tool(s) chosen:")
            logger.info(f"  BLASTn: {'Yes' if 'blastn' in options.tools else 'No'}")
            logger.info(f"  BLASTx: {'Yes' if 'blastx' in options.tools else 'No'}")
            logger.info(f"  BLASTp: {'Yes' if 'blastp' in options.tools else 'No'}")
            logger.info(f"  DIAMOND: {'Yes' if 'diamond' in options.tools else 'No'}")
            logger.info(f"  Bowtie2: {'Yes' if 'bowtie2' in options.tools else 'No'}")
            logger.info(f"  BWA: {'Yes' if 'bwa' in options.tools else 'No'}")
            logger.info(f"  Minimap2: {'Yes' if 'minimap2' in options.tools else 'No'}")

        # logger.info(f"E-value threshold: {evalue}")
        logger.info(f"Min query coverage: {options.query_min_coverage}%")
        logger.info(f"Min query identity: {options.query_min_identity}%")
        logger.info(f"Min detection coverage: {options.detection_min_coverage}%")
        logger.info(f"Min detection identity: {options.detection_min_identity}%")
        # Reflect Genes-FASTA + --genes-type effects in displayed run modes
        eff_run_dna = run_dna
        eff_run_protein = run_protein
        if getattr(options, 'sequence_type', None) == 'Genes-FASTA':
            genes_type_mapped = ('protein' if getattr(options, 'genes_type', None) == 'aa' else getattr(options, 'genes_type', None))
            if genes_type_mapped == 'protein':
                eff_run_dna = False
                eff_run_protein = True
            elif genes_type_mapped == 'dna':
                eff_run_dna = True
                eff_run_protein = True

        logger.info(f"Run DNA mode: {eff_run_dna}")
        logger.info(f"Run Protein mode: {eff_run_protein}")
        params_str = ", ".join(
            f"{tool}: {params}" for tool, params in tool_sensitivity_params.items()
        ) if tool_sensitivity_params else "None"
        logger.info(f"Sensitivity parameters: {options.sensitivity} - {params_str}")
        #logger.info("=" * 70)
        # Compute max_fasta_chunk_bytes early so per-sample loop can use it
        max_fasta_chunk_bytes = get_max_fasta_chunk_bytes(getattr(options, 'max_fasta_chunk_mb', 200.0))
        logger.info(f"Configured max FASTA chunk size: {getattr(options, 'max_fasta_chunk_mb', 200.0)} MiB ({max_fasta_chunk_bytes} bytes)")

        # Enforce rule: BLASTp only makes sense for Genes-FASTA protein queries.
        # If the user requested 'blastp' but did not select Genes-FASTA, warn and remove it
        # from the requested toolset to avoid accidental runs.
        try:
            if 'blastp' in (options.tools or []) and getattr(options, 'sequence_type', None) != 'Genes-FASTA':
                logger.warning(
                    "BLASTp was requested but sequence-type is not Genes-FASTA. BLASTp will be ignored unless -st Genes-FASTA is used.")
                options.tools = [t for t in options.tools if t != 'blastp']
        except Exception:
            pass

        # Input handling: support single input OR directory-based batch inputs
        from copy import deepcopy

        if getattr(options, 'input_dir', None) or getattr(options, 'input_subdirs', None):
            # Discover samples
            samples = []
            if getattr(options, 'input_dir', None):
                samples = discover_samples_from_input_dir(options.input_dir, options.sequence_type, logger)
            else:
                # Prevent discovery from scanning the output directory if it is inside input_subdirs
                samples = discover_samples_from_subdirs(options.input_subdirs, options.sequence_type, logger, exclude_paths=[options.output])

            if not samples:
                logger.error('No samples found in provided input directory')
                sys.exit(1)

            # Report discovered samples count and list to the logger
            try:
                logger.info(f"Discovered {len(samples)} sample(s) to process")
                for idx, (sname, spec) in enumerate(samples, start=1):
                    logger.info(f"  [{idx}/{len(samples)}] {sname}: {spec}")
            except Exception:
                # Ensure discovery logging never interrupts processing
                logger.debug("Failed to log discovered samples list")

            for sample_name, input_spec in samples:
                logger.info(f"Processing sample: {sample_name}")
                sample_opts = deepcopy(options)
                sample_out = os.path.join(options.output, sample_name)
                sample_opts.output = sample_out
                sample_opts.input = input_spec
                if getattr(sample_opts, 'temp_directory', None) is None:
                    sample_opts.temp_directory = sample_out

                os.makedirs(sample_out, exist_ok=True)
                # write run parameters manifest for sample
                try:
                    params = {
                        'tools': options.tools,
                        'query_min_coverage': options.query_min_coverage,
                        'query_min_identity': options.query_min_identity,
                        'detection_min_coverage': options.detection_min_coverage,
                        'detection_min_identity': options.detection_min_identity,
                        'detection_min_base_depth': options.detection_min_base_depth,
                        'detection_min_num_reads': options.detection_min_num_reads,
                        'sensitivity': getattr(options, 'sensitivity', None),
                        'evalue': getattr(options, 'evalue', None),
                        'min_bitscore': getattr(options, 'min_bitscore', None),
                        'tool_sensitivity_params': tool_sensitivity_params,
                        'tool_extra_params': tool_extra_params,
                        'sequence_type': getattr(sample_opts, 'sequence_type', None),
                        'genes_type': ('protein' if getattr(sample_opts, 'genes_type', None) == 'aa' else getattr(sample_opts, 'genes_type', None)),
                        'report_fasta': getattr(options, 'report_fasta', None),
                        'threads': getattr(options, 'threads', None)
                    }
                    with open(os.path.join(sample_out, 'run_parameters.json'), 'w') as pf:
                        json.dump(params, pf, indent=2)
                except Exception:
                    logger.debug("Failed to write run parameters for sample: %s", sample_out)

                try:
                    handle_all_input_files(sample_opts, logger)

                    workflow = Workflow(
                        input_fasta=sample_opts.input_fasta,
                        input_fastq=sample_opts.input_fastq,
                        output_dir=sample_out,
                        databases=databases,
                        threads=options.threads,
                        tool_sensitivity_params=tool_sensitivity_params,
                        extra_tool_params=tool_extra_params,
                        evalue=options.evalue,
                        min_bitscore=options.min_bitscore,
                        genes_type=('protein' if getattr(sample_opts, 'genes_type', None) == 'aa' else getattr(sample_opts, 'genes_type', None)),
                        blastx_task=options.blastx_task,
                        query_min_coverage=options.query_min_coverage,
                        query_min_identity=options.query_min_identity,
                        detection_min_coverage=options.detection_min_coverage,
                        detection_min_identity=options.detection_min_identity,
                        detection_min_base_depth=options.detection_min_base_depth,
                        detection_min_num_reads=options.detection_min_num_reads,
                        run_dna=run_dna,
                        run_protein=run_protein,
                        sequence_type=sample_opts.sequence_type,
                        report_fasta=sample_opts.report_fasta,
                        no_cleanup=sample_opts.no_cleanup,
                        verbose=sample_opts.verbose,
                        logger=logger,
                        temp_directory=sample_opts.temp_directory,
                        chunk_jobs=sample_opts.chunk_jobs,
                        chunk_threads_per_job=sample_opts.chunk_threads_per_job,
                        preserve_chunks=sample_opts.preserve_chunks,
                        max_fasta_chunk_bytes=max_fasta_chunk_bytes
                    )

                    results = workflow.run_workflow(sample_opts)
                    cleanup_all_temp_files(sample_opts, logger)
                finally:
                    # No per-sample handler to detach - global file handler remains in root output
                    pass

            # After all samples processed, combine per-sample detection matrices into combined matrices
            try:
                from .utils import combine_detection_matrices
            except Exception:
                from utils import combine_detection_matrices

            try:
                sample_names = [s for s, _ in samples]
                combine_detection_matrices(options.output, sample_names, logger)
            except Exception:
                logger.debug("Failed to combine detection matrices across samples")

            logger.info("Completed processing all samples")
            return

        else:
            # Single-run case: input handling using the already-created global file handler
            handle_all_input_files(options, logger)

        workflow = Workflow(
            input_fasta=options.input_fasta,
            input_fastq=options.input_fastq,
            output_dir=options.output,
            databases=databases,
            # resfinder_dbs=resfinder_dbs,
            # card_dbs=card_dbs,
            threads=options.threads,
            tool_sensitivity_params=tool_sensitivity_params,
            extra_tool_params=tool_extra_params,
            evalue=options.evalue,
            min_bitscore=options.min_bitscore,
            genes_type=('protein' if getattr(options, 'genes_type', None) == 'aa' else getattr(options, 'genes_type', None)),
            blastx_task=options.blastx_task,
            #max_target_seqs=options.max_target_seqs,
            #evalue=options.evalue,
            query_min_coverage=options.query_min_coverage,
            query_min_identity=options.query_min_identity,
            detection_min_coverage=options.detection_min_coverage,
            detection_min_identity=options.detection_min_identity,
            detection_min_base_depth=options.detection_min_base_depth,
            detection_min_num_reads=options.detection_min_num_reads,
            run_dna=run_dna,
            run_protein=run_protein,
            sequence_type=options.sequence_type,
            report_fasta=options.report_fasta,
            no_cleanup=options.no_cleanup,
            verbose=options.verbose,
            logger=logger,
            temp_directory=options.temp_directory,
            chunk_jobs=options.chunk_jobs,
            chunk_threads_per_job=options.chunk_threads_per_job,
            preserve_chunks=options.preserve_chunks,
            max_fasta_chunk_bytes=max_fasta_chunk_bytes
        )

        ###
        # write run parameters manifest to root output for single-run
        try:
            outdir = options.output
            Path(outdir).mkdir(parents=True, exist_ok=True)
            params = {
                'tools': options.tools,
                'query_min_coverage': options.query_min_coverage,
                'query_min_identity': options.query_min_identity,
                'detection_min_coverage': options.detection_min_coverage,
                'detection_min_identity': options.detection_min_identity,
                'detection_min_base_depth': options.detection_min_base_depth,
                'detection_min_num_reads': options.detection_min_num_reads,
                'sensitivity': getattr(options, 'sensitivity', None),
                'evalue': getattr(options, 'evalue', None),
                'min_bitscore': getattr(options, 'min_bitscore', None),
                'tool_sensitivity_params': tool_sensitivity_params,
                'tool_extra_params': tool_extra_params,
                'sequence_type': getattr(options, 'sequence_type', None),
                'genes_type': ('protein' if getattr(options, 'genes_type', None) == 'aa' else getattr(options, 'genes_type', None)),
                'report_fasta': getattr(options, 'report_fasta', None),
                'threads': getattr(options, 'threads', None)
            }
            with open(os.path.join(outdir, 'run_parameters.json'), 'w') as pf:
                json.dump(params, pf, indent=2)
        except Exception:
            logger.debug("Failed to write run parameters file to output: %s", getattr(options, 'output', ''))

        results = workflow.run_workflow(options)

        failed_tools = []
        all_failed = True

        for db_name, db_results in results.items():
            for tool, val in db_results.items():
                # Normalise possible shapes to (success, genes)
                if isinstance(val, tuple) and len(val) == 2:
                    success, genes = val
                elif isinstance(val, bool):
                    success, genes = val, set()
                elif isinstance(val, set):
                    success, genes = True, val
                else:
                    success, genes = bool(val), set()

                if not success:
                    failed_tools.append((db_name, tool, genes))
                else:
                    all_failed = False

        # Print specific statements for each failed tool
        if failed_tools:
            for db_name, tool, genes in failed_tools:
                gene_count = len(genes) if isinstance(genes, (set, list, tuple, dict)) else 0
                logger.info(f"Tool failure - {tool} (database: {db_name}): detected {gene_count} genes")
                logger.warning(f"  -> {tool} failed for {db_name}")

        if all_failed:
            logger.error("All tools failed - Something went catastrophically wrong. Exiting with error code.")
            sys.exit(1)


    finally:
        # Cleanup temp files
        cleanup_all_temp_files(options, logger)
        # Remove global file handler if created (single-run case)
        try:
            if file_handler:
                remove_file_handler(logger, file_handler)
        except Exception:
            pass
        end_time = datetime.now()
        elapsed = end_time - start_time
        logger.info(f"AMRfíor completed in {elapsed}.")


if __name__ == "__main__":
    main()