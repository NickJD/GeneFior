import argparse
import sys, os
import logging
import json
from datetime import datetime


script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

def _import_local(name: str):
    # Import a local module by name.
    return __import__(name)

# Import constants, databases, workflow, gene_stats, utils in a way that works both as package and script
if __package__:
    # Running as package (python -m GeneFior.GeneFior)
    from .constants import *
    from .databases import (
        RESFINDER_DATABASES, CARD_DATABASES, NCBI_DATABASES,
        gather_databases, gather_multiple_databases, parse_database_paths,
    )
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
        diamond_mode_label,
        sample_temp_directory,
        workflow_has_success,
    )
else:
    # Running as a script (python GeneFior.py) - import local modules by plain name
    try:
        from constants import *
    except Exception:
        # fallback: try package-style (if present in sys.path)
        from GeneFior.constants import *

    try:
        from databases import (
            RESFINDER_DATABASES, CARD_DATABASES, NCBI_DATABASES,
            gather_databases, gather_multiple_databases, parse_database_paths,
        )
    except Exception:
        from GeneFior.databases import (
            RESFINDER_DATABASES, CARD_DATABASES, NCBI_DATABASES,
            gather_databases, gather_multiple_databases, parse_database_paths,
        )

    try:
        from workflow import Workflow
    except Exception:
        from GeneFior.workflow import Workflow

    try:
        from gene_stats import GeneStats
    except Exception:
        from GeneFior.gene_stats import GeneStats

    try:
        from utils import (
            handle_all_input_files,
            cleanup_all_temp_files,
            discover_samples_from_input_dir,
            discover_samples_from_subdirs,
            add_file_handler_for_path,
            remove_file_handler,
            get_max_fasta_chunk_bytes,
            diamond_mode_label,
            sample_temp_directory,
            workflow_has_success,
        )
    except Exception:
        from GeneFior.utils import (
            handle_all_input_files, cleanup_all_temp_files, diamond_mode_label,
            sample_temp_directory, workflow_has_success,
        )


def main():
    parser = argparse.ArgumentParser(
        description='GeneFíor ' + GENEFIOR_VERSION + ' GeneFíor - The Multi-Tool Gene Detection Toolkit.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default tools (runs DNA & protein tools)
  GeneFior -i reads.fasta -st Single-FASTA --db-path <path-to-db> -o results

  # Select specific tools and output detected FASTA sequences
  GeneFior -i reads.fasta -st Single-FASTA --db-path <path-to-db> -o results \
    --tools diamond bowtie2 \
    --report_fasta detected

  # Custom thresholds, paired-fastq input, threads and dna-only mode
  GeneFior -i reads_R1.fastq,reads_R2.fastq -st Paired-FASTQ --db-path <path-to-db> \
    -o results -t 16 --d-min-cov 90 --d-min-id 85 \
    --dna-only
        """
    )

    # Required arguments
    required_group = parser.add_argument_group('Required selection')
    required_group.add_argument('-i', '--input', required=False,
                        help='Input FASTA/FASTAQ file(s) with sequences to analyse - Separate FASTQ R1 and R2 '
                             'with a comma for Paired-FASTQ or single file path for Single-FASTA - .gz files accepted')
    # Optional directory-based inputs: either a flat directory containing multiple FASTA/FASTQ samples
    # or a directory of sample subdirectories (each subdir contains one sample FASTA or paired FASTQ files)
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
    required_group.add_argument('--db-path', dest='user_db_path', type=str, required=False, default=None,
                          help='Path to a user-provided database directory in the format produced by build_databases.sh. '
                               'Supply multiple independent databases as a comma-separated list. '
                               'Not required when --hmmer-only is used with --hmmer-db and/or --hmmer-dna-db.')
    required_group.add_argument('-o', '--output', required=True,
                        help='Output directory for results')

    # Output selection
    output_group = parser.add_argument_group('Output selection')
    output_group.add_argument('--report-fasta',
                            choices=[
                                'all', 'detected', 'detected-all',
                                'evidence', 'evidence-all',
                                'candidate', 'candidate-all',
                                'exact', 'exact-all',
                            ], #, 'hmmer_dna', 'hmmer_protein'],
                            default=None,
                            dest='report_fasta',
                            help='Specify whether to output sequences that "mapped" to genes.'
                                 '"all" should only be used for deep investigation/debugging.'
                                 '"detected" will report the reads that passed detection thresholds for each detected gene.'
                                 '"detected-all" reports all mapped reads for each evidence gene.'
                                 '"evidence" and "evidence-all" are explicit aliases for the detected modes.'
                                 '"candidate" and "candidate-all" restrict output to candidate allele calls.'
                                 '"exact" and "exact-all" restrict output to exact nucleotide alleles.  (default: None)')

    # Tool selection
    tool_group = parser.add_argument_group('Tool selection')
    tool_group.add_argument('--tools', nargs='+',
                            choices=['blastn', 'blastx', 'blastp', 'diamond', 'bowtie2', 'bwa', 'minimap2',
                                     'hmmer_protein', 'hmmer_dna', 'all'],
                            default=['diamond', 'bowtie2', 'bwa', 'minimap2'],
                            help='Specify which tools to run - "all" will run all tools'
                                 ' (default: all except blastx/p as they can be slow!!)')

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
    gene_detection_group.add_argument(
        '--detection-system', '--detection-mode',
        choices=['qualified', 'legacy-relaxed'],
        default='qualified',
        help='Detection interpretation: "qualified" uses evidence, family/mosaic, '
             'and exact-allele resolution (default); "legacy-relaxed" reproduces '
             'the original direct threshold-only detector.')
    gene_detection_group.add_argument('--evidence-corroborating-depth', type=int, default=2,
                              help='Depth required across the gene for a robust read-based allele call (default: 2)')
    gene_detection_group.add_argument('--evidence-exact-identity', type=float, default=100.0,
                              help='Deprecated identity setting retained for compatibility; literal exact calls require 100%% identity and full candidate-depth coverage')
    gene_detection_group.add_argument('--evidence-candidate-depth', type=int, default=3,
                              help='Minimum per-base depth required across 100%% of a nucleotide allele candidate (default: 3)')
    gene_detection_group.add_argument('--evidence-candidate-identity', type=float, default=98.0,
                              help='Minimum mean identity for a nucleotide allele candidate (default: 98.0)')
    gene_detection_group.add_argument('--evidence-max-internal-gap-bp', type=int, default=15,
                              help='Largest unsupported internal gap allowed for an exact call (default: 15 bp)')
    gene_detection_group.add_argument('--evidence-max-internal-gap-fraction', type=float, default=0.02,
                              help='Gene-length fraction allowed as an internal gap (default: 0.02)')
    gene_detection_group.add_argument('--evidence-min-unique-reads', type=int, default=2,
                              help='Minimum competitively unique reads for an exact allele (default: 2)')
    gene_detection_group.add_argument('--evidence-min-unique-fraction', type=float, default=0.10,
                              help='Minimum fraction of passing reads uniquely supporting the allele (default: 0.10)')
    gene_detection_group.add_argument('--evidence-ambiguity-fraction', type=float, default=0.50,
                              help='Ambiguous-best read fraction that triggers an allele warning (default: 0.50)')
    gene_detection_group.add_argument('--evidence-score-tie', type=float, default=1.0,
                              help='Absolute bit-score difference treated as a competitive tie (default: 1.0)')

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
                                   help="Specify extra parameters for a tool as TOOL=\"args\". Can be provided multiple times."
                                        " Example: --tool-param 'diamond=--more-sensitive -e 1e-5'")
    # Tool-wide evalue: belongs with tool-specific parameters
    tool_params_group.add_argument('-e', '--evalue', type=float, default=None,
                                   help='E-value threshold to pass to BLASTn/x and DIAMOND (if not set, the tools\' default is used)')
    tool_params_group.add_argument('--min-bitscore', type=float, default=None,
                                   help='Minimum bitscore threshold to apply to BLAST/DIAMOND hits (if not set, no bitscore filtering is applied)')

    # HMMER-specific arguments
    hmmer_group = parser.add_argument_group('HMMER Parameters')
    hmmer_group.add_argument('--hmmer-only', dest='hmmer_only', action='store_true', default=False,
                             help='Run ONLY HMMER searches (hmmer_protein and/or hmmer_dna). '
                                  'When set, --db-path is not required — provide the HMM database(s) via '
                                  '--hmmer-db and/or --hmmer-dna-db. All other alignment tools are disabled.')
    hmmer_group.add_argument('--hmmer-db', dest='hmmer_db', default=None,
                             help='Path to a protein HMM database file (.hmm) for hmmsearch. '
                                  'Enables HMMER protein search even if "hmmer_protein" is not listed in --tools. '
                                  'When used with --hmmer-only, this is the sole database required. '
                                  'Example: /path/to/biorisk.hmm')
    hmmer_group.add_argument('--hmmer-dna-db', dest='hmmer_dna_db', default=None,
                             help='Path to a DNA HMM database file (.hmm) for nhmmer. '
                                  'Enables HMMER DNA search even if "hmmer_dna" is not listed in --tools. '
                                  'When used with --hmmer-only, this is the sole database required.')
    hmmer_group.add_argument('--hmmer-annotations', dest='hmmer_annotations', default=None,
                             help='Path to a CSV file mapping HMM profile IDs to descriptions '
                                  '(columns: ID, Description, Must flag). Used with --hmmer-db to annotate '
                                  'detected genes. Example: /path/to/biorisk_annotations.csv')
    hmmer_group.add_argument('--hmmer-evalue', dest='hmmer_evalue', type=float, default=None,
                             help='E-value threshold specifically for HMMER (hmmsearch / nhmmer). '
                                  'Overrides the global -e/--evalue for HMMER only, allowing different '
                                  'stringency for HMMER vs BLAST/DIAMOND. If unset, falls back to -e if '
                                  'that is set, otherwise uses HMMER\'s own default (E=10, very permissive). '
                                  'Recommended: 1e-5 for Genes-FASTA, 1e-3 for read inputs. '
                                  'Not used when --hmmer-threshold-mode is tc/ga/nc.')
    hmmer_group.add_argument('--hmmer-threshold-mode', dest='hmmer_threshold_mode',
                             choices=['evalue', 'tc', 'ga', 'nc'],
                             default='evalue',
                             help='Threshold mode for HMMER: "evalue" (default) uses the E-value threshold '
                                  '(-e or --hmmer-evalue). "tc" / "ga" / "nc" use the per-profile trusted / '
                                  'gathering / noise cutoff scores embedded in the HMM database '
                                  '(--cut_tc / --cut_ga / --cut_nc flags passed to hmmsearch/nhmmer). '
                                  'Recommended: "tc" when your HMM database was built with per-profile '
                                  'TC scores (e.g. ibbis biorisk.hmm). When using tc/ga/nc the --hmmer-evalue '
                                  'threshold is ignored and E-value post-filtering is disabled.')
    hmmer_group.add_argument('--no-must-flag', dest='no_must_flag', action='store_true', default=False,
                             help='Disable the must-flag override for HMMER results. '
                                  'By default, genes marked "Must flag" in the annotations CSV are always '
                                  'reported even if they fail the coverage or min-reads thresholds. '
                                  'Setting this flag turns off that behaviour so must-flag genes are subject '
                                  'to the same detection gates as all other genes and no MUST-FLAG ALERT '
                                  'is emitted in the log.')

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
    runtime_group.add_argument(
        '--fastp-trim', '--fastp', '--trim-reads',
        dest='fastp_trim',
        action='store_true',
        default=False,
        help='For Paired-FASTQ input, run fastp automatic adapter and quality trimming before analysis'
    )
    runtime_group.add_argument('--no_cleanup',  action='store_true',)
    runtime_group.add_argument( '--verbose', action='store_true',)

    misc_group = parser.add_argument_group('Miscellaneous Parameters')
    misc_group.add_argument('-v','--version', action='version',
                            version='GeneFíor ' + GENEFIOR_VERSION,
                            help='Show program version and exit')

    options = parser.parse_args()

    # --hmmer-only special case: restrict tools to HMMER only and make --db-path optional
    if getattr(options, 'hmmer_only', False):
        if not getattr(options, 'hmmer_db', None) and not getattr(options, 'hmmer_dna_db', None):
            print("Error: --hmmer-only requires at least one of --hmmer-db or --hmmer-dna-db to be specified", file=sys.stderr)
            sys.exit(1)
        # Override tools to HMMER only
        options.tools = []
        if getattr(options, 'hmmer_db', None):
            options.tools.append('hmmer_protein')
        if getattr(options, 'hmmer_dna_db', None):
            options.tools.append('hmmer_dna')
        # Provide a dummy db-path placeholder so downstream code doesn't fail
        if not getattr(options, 'user_db_path', None):
            options.user_db_path = '__hmmer_only__'
    else:
        # Normal path: --db-path is required
        if not getattr(options, 'user_db_path', None):
            print("Error: --db-path is required unless --hmmer-only is specified", file=sys.stderr)
            sys.exit(1)

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

    options.temp_directory_user_specified = bool(
        getattr(options, 'temp_directory', None)
    )
    # Ensure options.temp_directory is set (default to output dir).
    if getattr(options, 'temp_directory', None) is None:
        options.temp_directory = options.output

    ## Setup logging
    start_time = datetime.now()
    from pathlib import Path
    log_file = Path(options.output) / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = logging.getLogger('GeneFíor')
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

    if options.report_fasta in (None, 'None'):
        options.report_fasta = None


    try:
       # Tool selection
        if options.tools == ['all']:
            # Expand "all" differently depending on sequence type. For Genes-FASTA we only
            # include the tools that make sense for gene FASTA inputs (no read-mappers).
            if getattr(options, 'sequence_type', None) == 'Genes-FASTA':
                genes_type = getattr(options, 'genes_type', None)
                if genes_type == 'aa':
                    # Amino-acid genes: prefer BLASTp and DIAMOND (blastp)
                    options.tools = ['blastp', 'diamond']
                elif genes_type == 'dna':
                    # DNA genes: prefer BLASTn and DIAMOND (blastx)
                    options.tools = ['blastn', 'blastx', 'diamond']
                else:
                    # Unknown genes_type: include nucleotide and translated searches plus DIAMOND
                    options.tools = ['blastn', 'blastx', 'diamond']
            else:
                # Non-Genes inputs: include all DNA/protein/mapping tools
                options.tools = ['blastn', 'blastx', 'diamond', 'bowtie2', 'bwa', 'minimap2']

        # If --hmmer-db is provided, ensure hmmer_protein is in tools (auto-enable)
        if getattr(options, 'hmmer_db', None) and 'hmmer_protein' not in options.tools:
            options.tools = list(options.tools) + ['hmmer_protein']
        if getattr(options, 'hmmer_dna_db', None) and 'hmmer_dna' not in options.tools:
            options.tools = list(options.tools) + ['hmmer_dna']

        # Load one or more independent user-provided databases.
        user_dbs = None
        user_database_sets = {}
        user_database_paths = []
        # Load database from provided path
        if getattr(options, 'hmmer_only', False) and options.user_db_path == '__hmmer_only__':
            # HMMER-only mode: no conventional database directory needed.
            # HMM files will be injected below via --hmmer-db / --hmmer-dna-db.
            user_dbs = {}
        else:
            try:
                user_database_paths = parse_database_paths(options.user_db_path)
                user_database_sets = gather_multiple_databases(
                    options.user_db_path, options.tools
                )
            except ValueError as exc:
                print(f"Error: {exc}", file=sys.stderr)
                sys.exit(1)
            if not user_database_sets:
                print(
                    "Error: Please provide at least one valid database directory "
                    "using --db-path",
                    file=sys.stderr,
                )
                sys.exit(1)
            # Preserve the existing variable for HMM injection and activity logs.
            user_dbs = next(iter(user_database_sets.values()))
        user_database_source_paths = {}
        if user_database_paths:
            user_database_source_paths = dict(
                zip(user_database_sets.keys(), user_database_paths)
            )

        # Inject HMM databases provided via --hmmer-db / --hmmer-dna-db / --hmmer-annotations
        if getattr(options, 'hmmer_db', None):
            if os.path.isfile(options.hmmer_db):
                user_dbs['hmmer_protein'] = options.hmmer_db
            else:
                logger.warning(f"--hmmer-db path not found or not a file: {options.hmmer_db}; ignoring.")
        if getattr(options, 'hmmer_dna_db', None):
            if os.path.isfile(options.hmmer_dna_db):
                user_dbs['hmmer_dna'] = options.hmmer_dna_db
            else:
                logger.warning(f"--hmmer-dna-db path not found or not a file: {options.hmmer_dna_db}; ignoring.")
        if getattr(options, 'hmmer_annotations', None):
            if os.path.isfile(options.hmmer_annotations):
                user_dbs['hmmer_annotations'] = options.hmmer_annotations
            else:
                logger.warning(f"--hmmer-annotations path not found: {options.hmmer_annotations}; ignoring.")
        elif user_dbs.get('hmmer_protein') and 'hmmer_annotations' not in user_dbs:
            # Auto-detect annotations CSV in the same directory as the HMM file
            import glob as _glob
            _hmm_dir = os.path.dirname(user_dbs['hmmer_protein'])
            _csv_candidates = _glob.glob(os.path.join(_hmm_dir, '*.csv'))
            if _csv_candidates:
                user_dbs['hmmer_annotations'] = _csv_candidates[0]
                logger.info(f"  Auto-detected HMMER annotations: {_csv_candidates[0]}")

        if getattr(options, 'hmmer_only', False):
            databases = {'user-provided-db': user_dbs}
        else:
            databases = dict(user_database_sets)
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
                    # Some packaged DB constants expose 'blastx' for protein DBs but not 'blastp'.
                    # If the user explicitly requested 'blastp' we should expose a fallback mapping
                    # from an available 'blastx' entry so BLASTp can be used for protein Genes-FASTA.
                    if isinstance(db_map, dict):
                        if 'blastp' not in db_map and 'blastx' in db_map:
                            db_map['blastp'] = db_map.get('blastx')

                    bn = 'Yes' if db_map.get('blastn') else 'No'
                    bx = 'Yes' if db_map.get('blastx') else 'No'
                    bp = 'Yes' if db_map.get('blastp') else 'No'
                    diamond_label = diamond_mode_label(
                        options.sequence_type, getattr(options, 'genes_type', None)
                    )
                    dm = 'Yes' if db_map.get('diamond') else 'No'
                    logger.info(
                        f"  {db_name}: BLASTn: {bn}; BLASTx: {bx}; "
                        f"BLASTp: {bp}; {diamond_label}: {dm}"
                    )
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

        # Map user-friendly sensitivity presets to tool-specific flags.
        # These can be overridden per-tool using --tool-param TOOL="args".
        SENSITIVITY_MAP = {
            'very-fast': {'diamond': '--fast', 'bowtie2': '--very-fast'},
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
        logger.info("Genefíor - The Gene Detection toolkit: " + GENEFIOR_VERSION)
        logger.info("Started at: " + start_time.strftime("%Y-%m-%d %H:%M:%S"))
        if getattr(options, 'hmmer_only', False):
            logger.info("Mode: HMMER-ONLY — all other alignment tools are disabled.")
        logger.info("=" * 70)
        ###
        # Log input files (handle new FASTA/FASTQ possibilities)
        # if getattr(options, "input_fasta", None):
        #     logger.info(f"Input FASTA: {options.input_fasta}")
        # else:
        #     logger.info("Input FASTA: None")
        #
        # if getattr(options, "input_fastq", None) is None:
        #     logger.info("Input FASTQ: None")
        # else:
        #     if getattr(options, "input_fastq_is_paired", False):
        #         logger.info(f"Input FASTQ (paired): {options.input_fastq[0]}, {options.input_fastq[1]}")
        #     else:
        #         logger.info(f"Input FASTQ (single): {options.input_fastq}")
        ###
        ###

        logger.info(f"Threads: {options.threads}")
        logger.info(f"Database(s) chosen:")
        if getattr(options, 'hmmer_only', False):
            logger.info("  User-Provided: N/A (HMMER-only mode)")
        else:
            for db_name in databases:
                logger.info(
                    f"  {db_name}: "
                    f"{user_database_source_paths.get(db_name, 'unknown path')}"
                )

        # If user selected Genes-FASTA input, clarify that read-mappers will be skipped
        if getattr(options, 'sequence_type', None) == 'Genes-FASTA':
            logger.info("Sequence type: Genes-FASTA (full-length gene FASTA). Read-mappers (bowtie2, bwa, minimap2) will be skipped by the workflow.")
            # Map CLI genes_type 'aa' -> internal 'protein'
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
            diamond_label = diamond_mode_label(
                options.sequence_type, getattr(options, 'genes_type', None)
            )
            logger.info(f"  {diamond_label}: {'Yes' if 'diamond' in options.tools else 'No'}")
            # Indicate mappers are present in the selection but will be skipped
            logger.info(f"  Bowtie2: {'Present (will be SKIPPED)' if 'bowtie2' in options.tools else 'No'}")
            logger.info(f"  BWA: {'Present (will be SKIPPED)' if 'bwa' in options.tools else 'No'}")
            logger.info(f"  Minimap2: {'Present (will be SKIPPED)' if 'minimap2' in options.tools else 'No'}")
            logger.info(f"  HMMER-Protein: {'Yes' if 'hmmer_protein' in options.tools or user_dbs.get('hmmer_protein') else 'No'}")
            logger.info(f"  HMMER-DNA: {'Yes' if 'hmmer_dna' in options.tools or user_dbs.get('hmmer_dna') else 'No'}")
        else:
            logger.info(f" Tool(s) chosen:")
            logger.info(f"  BLASTn: {'Yes' if 'blastn' in options.tools else 'No'}")
            logger.info(f"  BLASTx: {'Yes' if 'blastx' in options.tools else 'No'}")
            logger.info(f"  BLASTp: {'Yes' if 'blastp' in options.tools else 'No'}")
            diamond_label = diamond_mode_label(
                options.sequence_type, getattr(options, 'genes_type', None)
            )
            logger.info(f"  {diamond_label}: {'Yes' if 'diamond' in options.tools else 'No'}")
            logger.info(f"  Bowtie2: {'Yes' if 'bowtie2' in options.tools else 'No'}")
            logger.info(f"  BWA: {'Yes' if 'bwa' in options.tools else 'No'}")
            logger.info(f"  Minimap2: {'Yes' if 'minimap2' in options.tools else 'No'}")
            logger.info(f"  HMMER-Protein: {'Yes' if 'hmmer_protein' in options.tools or user_dbs.get('hmmer_protein') else 'No'}")
            logger.info(f"  HMMER-DNA: {'Yes' if 'hmmer_dna' in options.tools or user_dbs.get('hmmer_dna') else 'No'}")

        # logger.info(f"E-value threshold: {evalue}")
        logger.info(f"Min query coverage: {options.query_min_coverage}%")
        logger.info(f"Min query identity: {options.query_min_identity}%")
        logger.info(f"Min detection coverage: {options.detection_min_coverage}%")
        logger.info(f"Min detection identity: {options.detection_min_identity}%")
        logger.info(f"Detection system: {options.detection_system}")
        # If Genes-FASTA was requested and user provided genes_type, reflect its effect in these logs
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
        # Log HMMER-specific threshold configuration
        _hmmer_tmode = getattr(options, 'hmmer_threshold_mode', 'evalue')
        _hmmer_ev = getattr(options, 'hmmer_evalue', None)
        if 'hmmer_protein' in getattr(options, 'tools', []) or 'hmmer_dna' in getattr(options, 'tools', []):
            logger.info(f"HMMER threshold mode: {_hmmer_tmode}")
            if _hmmer_tmode == 'evalue':
                _ev_display = _hmmer_ev if _hmmer_ev is not None else (
                    getattr(options, 'evalue', None) if getattr(options, 'evalue', None) is not None
                    else "HMMER default (E=10) — consider setting --hmmer-evalue for tighter control"
                )
                logger.info(f"  HMMER E-value: {_ev_display}")
            else:
                logger.info(f"  HMMER using per-profile {_hmmer_tmode.upper()} cutoffs (--cut_{_hmmer_tmode}); E-value post-filtering disabled")
            _mf_status = "ENABLED" if not getattr(options, 'no_must_flag', False) else "DISABLED (--no-must-flag)"
            logger.info(f"  HMMER must-flag override: {_mf_status}")
        #logger.info("=" * 70)
        # Compute max_fasta_chunk_bytes early so batch/sample loop can reuse it
        max_fasta_chunk_bytes = get_max_fasta_chunk_bytes(getattr(options, 'max_fasta_chunk_mb', 200.0))
        logger.info(f"Configured max FASTA chunk size: {getattr(options, 'max_fasta_chunk_mb', 200.0)} MiB ({max_fasta_chunk_bytes} bytes)")

        # Enforce rule: BLASTp only makes sense for Genes-FASTA protein queries.
        # If the user requested 'blastp' but did not select Genes-FASTA, warn and remove it
        # from the requested toolset to avoid accidental runs.
        try:
            if 'blastp' in (options.tools or []) and getattr(options, 'sequence_type', None) != 'Genes-FASTA':
                logger.warning("BLASTp was requested but sequence-type is not Genes-FASTA. BLASTp will be ignored unless -st Genes-FASTA is used.")
                options.tools = [t for t in options.tools if t != 'blastp']
        except Exception:
            pass

        # Input handling: support single input OR directory-based batch inputs
        # If directory-based inputs are provided, discover samples and iterate per-sample.
        from copy import deepcopy

        if getattr(options, 'input_dir', None) or getattr(options, 'input_subdirs', None):
            # Discover samples
            samples = []  # list of (sample_name, input_spec)
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
                logger.debug("Failed to log discovered samples list")

            # Process each sample in turn, placing outputs under output/<sample_name>
            aggregate_results = {}
            failed_samples = []
            for sample_name, input_spec in samples:
                logger.info(f"Processing sample: {sample_name}")
                sample_opts = deepcopy(options)
                sample_out = os.path.join(options.output, sample_name)
                sample_opts.output = sample_out
                # set input spec
                sample_opts.input = input_spec
                sample_opts.temp_directory = sample_temp_directory(
                    options.temp_directory,
                    options.output,
                    sample_out,
                    sample_name,
                    options.temp_directory_user_specified,
                )

                os.makedirs(sample_out, exist_ok=True)
                # Write run parameters manifest for this sample so downstream tools (e.g. GeneFíor-Gene-Stats)
                # can read and apply the same query/detection thresholds used for this run.
                try:
                    params = {
                        'tools': options.tools,
                        'query_min_coverage': options.query_min_coverage,
                        'query_min_identity': options.query_min_identity,
                        'detection_min_coverage': options.detection_min_coverage,
                        'detection_min_identity': options.detection_min_identity,
                        'detection_min_base_depth': options.detection_min_base_depth,
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
                        'sensitivity': getattr(options, 'sensitivity', None),
                        'evalue': getattr(options, 'evalue', None),
                        'min_bitscore': getattr(options, 'min_bitscore', None),
                        'tool_sensitivity_params': tool_sensitivity_params,
                        'tool_extra_params': tool_extra_params,
                        'sequence_type': getattr(sample_opts, 'sequence_type', None),
                        'genes_type': ('protein' if getattr(sample_opts, 'genes_type', None) == 'aa' else getattr(sample_opts, 'genes_type', None)),
                        'report_fasta': getattr(options, 'report_fasta', None),
                        'threads': getattr(options, 'threads', None),
                        'fastp_trim': getattr(options, 'fastp_trim', False),
                        'hmmer_must_flag': not getattr(options, 'no_must_flag', False),
                    }
                    with open(os.path.join(sample_out, 'run_parameters.json'), 'w') as pf:
                        json.dump(params, pf, indent=2)
                except Exception:
                    logger.debug("Failed to write run parameters file for sample: %s", sample_out)

                try:
                    # Perform per-sample input handling and workflow run
                    handle_all_input_files(sample_opts, logger)
                    # Recreate workflow for sample
                    workflow = Workflow(
                        input_fasta=sample_opts.input_fasta,
                        input_fastq=sample_opts.input_fastq,
                        output_dir=sample_out,
                        databases=databases,
                        threads=options.threads,
                        tool_sensitivity_params=tool_sensitivity_params,
                        evalue=options.evalue,
                        min_bitscore=options.min_bitscore,
                        extra_tool_params=tool_extra_params,
                        blastx_task=options.blastx_task,
                        query_min_coverage=options.query_min_coverage,
                        query_min_identity=options.query_min_identity,
                        detection_min_coverage=options.detection_min_coverage,
                        detection_min_identity=options.detection_min_identity,
                        detection_min_base_depth=options.detection_min_base_depth,
                        detection_min_num_reads=options.detection_min_num_reads,
                        detection_system=options.detection_system,
                        run_dna=run_dna,
                        run_protein=run_protein,
                        sequence_type=sample_opts.sequence_type,
                        genes_type=('protein' if getattr(sample_opts, 'genes_type', None) == 'aa' else getattr(sample_opts, 'genes_type', None)),
                        report_fasta=sample_opts.report_fasta,
                        no_cleanup=sample_opts.no_cleanup,
                        verbose=sample_opts.verbose,
                        logger=logger,
                        temp_directory=sample_opts.temp_directory,
                        chunk_jobs=sample_opts.chunk_jobs,
                        chunk_threads_per_job=sample_opts.chunk_threads_per_job,
                        preserve_chunks=sample_opts.preserve_chunks,
                        max_fasta_chunk_bytes=max_fasta_chunk_bytes,
                        hmmer_evalue=getattr(options, 'hmmer_evalue', None),
                        hmmer_threshold_mode=getattr(options, 'hmmer_threshold_mode', 'evalue'),
                        hmmer_must_flag=not getattr(options, 'no_must_flag', False),
                        tools=options.tools,
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
                    )

                    results = workflow.run_workflow(sample_opts)
                    aggregate_results[sample_name] = results
                    if not workflow_has_success(results):
                        failed_samples.append(sample_name)
                        logger.error(
                            f"No selected tool completed successfully for sample "
                            f"'{sample_name}'"
                        )
                finally:
                    cleanup_all_temp_files(sample_opts, logger)

            # After processing all samples, combine per-sample detection matrices into combined matrices
            try:
                from .utils import combine_detection_matrices
            except Exception:
                from utils import combine_detection_matrices

            try:
                sample_names = [s for s, _ in samples]
                combine_detection_matrices(options.output, sample_names, logger)
            except Exception:
                logger.debug("Failed to combine detection matrices across samples")

            # After combining, exit (we do not continue to run a second aggregate workflow)
            logger.info("Completed processing all samples")
            if failed_samples:
                logger.error(
                    "Batch run failed for sample(s): " + ', '.join(failed_samples)
                )
                sys.exit(1)
            return

        else:
            # Single-sample path (existing behavior) - input handling using the already-created global file handler
            handle_all_input_files(options, logger)

        # Run Workflow
        # Convert user-specified chunk size (MiB) into bytes and log it
        max_fasta_chunk_bytes = get_max_fasta_chunk_bytes(getattr(options, 'max_fasta_chunk_mb', 200.0))
        logger.info(f"Configured max FASTA chunk size: {getattr(options, 'max_fasta_chunk_mb', 200.0)} MiB ({max_fasta_chunk_bytes} bytes)")

        workflow = Workflow(
            input_fasta=options.input_fasta,
            input_fastq=options.input_fastq,
            output_dir=options.output,
            databases=databases,
            threads=options.threads,
            tool_sensitivity_params=tool_sensitivity_params,
            evalue=options.evalue,
            min_bitscore=options.min_bitscore,
            extra_tool_params=tool_extra_params,
            blastx_task=options.blastx_task,
            query_min_coverage=options.query_min_coverage,
            query_min_identity=options.query_min_identity,
            detection_min_coverage=options.detection_min_coverage,
            detection_min_identity=options.detection_min_identity,
            detection_min_base_depth=options.detection_min_base_depth,
            detection_min_num_reads=options.detection_min_num_reads,
            detection_system=options.detection_system,
            run_dna=run_dna,
            run_protein=run_protein,
            sequence_type=options.sequence_type,
            genes_type=('protein' if getattr(options, 'genes_type', None) == 'aa' else getattr(options, 'genes_type', None)),
            report_fasta=options.report_fasta,
            no_cleanup=options.no_cleanup,
            verbose=options.verbose,
            logger=logger,
            temp_directory=options.temp_directory,
            chunk_jobs=options.chunk_jobs,
            chunk_threads_per_job=options.chunk_threads_per_job,
            preserve_chunks=options.preserve_chunks,
            max_fasta_chunk_bytes=max_fasta_chunk_bytes,
            hmmer_evalue=getattr(options, 'hmmer_evalue', None),
            hmmer_threshold_mode=getattr(options, 'hmmer_threshold_mode', 'evalue'),
            hmmer_must_flag=not getattr(options, 'no_must_flag', False),
            tools=options.tools,
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
        )

        ###
        # Write run parameters manifest at the root output for single-run mode
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
                'sensitivity': getattr(options, 'sensitivity', None),
                'evalue': getattr(options, 'evalue', None),
                'min_bitscore': getattr(options, 'min_bitscore', None),
                'tool_sensitivity_params': tool_sensitivity_params,
                'tool_extra_params': tool_extra_params,
                    'sequence_type': getattr(options, 'sequence_type', None),
                'genes_type': ('protein' if getattr(options, 'genes_type', None) == 'aa' else getattr(options, 'genes_type', None)),
                'report_fasta': getattr(options, 'report_fasta', None),
                'threads': getattr(options, 'threads', None),
                'fastp_trim': getattr(options, 'fastp_trim', False),
                'hmmer_must_flag': not getattr(options, 'no_must_flag', False),
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
                display_tool = (
                    diamond_mode_label(
                        options.sequence_type,
                        getattr(options, 'genes_type', None),
                    )
                    if tool == 'DIAMOND-AA'
                    else tool
                )
                gene_count = len(genes) if isinstance(genes, (set, list, tuple, dict)) else 0
                logger.info(f"Tool failure - {display_tool} (database: {db_name}): detected {gene_count} genes")
                logger.warning(f"  -> {display_tool} failed for {db_name}")

        if all_failed:
            logger.error("All selected tools failed; no valid detection run was produced.")
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
        logger.info(f"GeneFíor completed in {elapsed}.")


if __name__ == "__main__":
    main()
