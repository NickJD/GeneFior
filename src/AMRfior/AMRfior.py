import argparse
import sys
import os

try:
    from .constants import *
    from .databases import RESFINDER_DATABASES, CARD_DATABASES, NCBI_DATABASES, gather_databases
    from .workflow import AMRWorkflow
    from .gene_stats import GeneStats
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    from databases import RESFINDER_DATABASES, CARD_DATABASES, NCBI_DATABASES, gather_databases
    from workflow import AMRWorkflow
    from gene_stats import GeneStats
    from utils import *

def main():
    parser = argparse.ArgumentParser(
        description='AMRfíor ' + AMRFIOR_VERSION + ' - The Multi-Tool AMR Gene Detection Toolkit.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default tools (runs DNA & protein tools)
  AMRfior -i reads.fasta -st Single-FASTA -o results/

  # Select specific tools and output detected FASTA sequences
  AMRfior -i reads.fasta -st Single-FASTA -o results/ \
    --tools diamond bowtie2 \
    --report_fasta detected

  # Custom thresholds, paired-fastq input, threads and dna-only mode
  AMRfior -i reads_R1.fastq,reads_R2.fastq -st Paired-FASTQ -o results/ \
    -t 16 --d-min-cov 90 --d-min-id 85 \
    --dna-only
        """
    )

    # Required arguments
    required_group = parser.add_argument_group('Required selection')
    required_group.add_argument('-i', '--input', required=True,
                        help='Input FASTA/FASTAQ file(s) with sequences to analyse - Separate FASTQ R1 and R2 '
                             'with a comma for Paired-FASTQ or single file path for Single-FASTA - .gz files accepted')
    required_group.add_argument('-st', '--sequence-type', required=True,
                        choices=['Single-FASTA', 'Paired-FASTQ'],
                        help='Specify the input Sequence Type: Single-FASTA or Paired-FASTQ (R1+R2) - Will '
                             'convert Paired-FASTQ to single combined FASTA for BLAST and DIAMOND analyses (SLOW)')
    required_group.add_argument('-o', '--output', required=True,
                        help='Output directory for results')

    # Output selection
    output_group = parser.add_argument_group('Output selection')
    output_group.add_argument('--report_fasta',
                            choices=['None', 'all', 'detected', 'detected-all'], #, 'hmmer_dna', 'hmmer_protein'],
                            default=[None], #, 'hmmer_dna','hmmer_protein'],
                            help='Specify whether to output sequences that "mapped" to genes.'
                                 '"all" should only be used for deep investigation/debugging.'
                                 '"detected" will report the reads that passed detection thresholds for each detected gene.'
                                 '"detected-all" will report all reads for each detected gene.  (default: None)')

    # Tool selection
    tool_group = parser.add_argument_group('Tool selection')
    tool_group.add_argument('--tools', nargs='+',
                            choices=['blastn', 'blastx', 'diamond', 'bowtie2', 'bwa', 'minimap2', 'all'], #, 'hmmer_dna', 'hmmer_protein'],
                            default=['diamond', 'bowtie2', 'bwa', 'minimap2'], #, 'hmmer_dna','hmmer_protein'],
                            help='Specify which tools to run - "all" will run all tools'
                                 ' (default: all except blastx/n as it is very slow!!)')

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
                                      help='Minimum coverage threshold in percent (default: 40.0)')

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
                            choices=['default', 'conservative', 'sensitive', 'very-sensitive'],
                            help='Preset sensitivity levels - default means each tool uses its own default settings and '
                                 'very-sensitive applies DIAMONDs --ultra-sensitive and Bowtie2s'
                                 ' --very-sensitive-local presets')

    # Tool-specific parameters
    tool_params_group = parser.add_argument_group('Tool-Specific Parameters')
    tool_params_group.add_argument('--minimap2-preset', default='sr',
                                   choices=['sr', 'map-ont', 'map-pb', 'map-hifi'],
                                   help='Minimap2 preset: sr=short reads, map-ont=Oxford Nanopore, '
                                        'map-pb=PacBio, map-hifi=PacBio HiFi (default: sr)')
    # tool_params_group.add_argument('-e', '--evalue', type=float, default=1e-10,
    #                                help='E-value threshold (default: 1e-10)')

    # Runtime parameters
    runtime_group = parser.add_argument_group('Runtime Parameters')
    runtime_group.add_argument('-t', '--threads', type=int, default=4,
                              help='Number of threads to use (default: 4)')
    runtime_group.add_argument('--no_cleanup',  action='store_true',)
    runtime_group.add_argument('-v', '--verbose', action='store_true',)

    options = parser.parse_args()

    if options.report_fasta == ['None']:
        options.report_fasta = None

    # FASTA/FASTQ handling
    if options.sequence_type == 'Paired-FASTQ':
        if any(tool in ('blastn', 'blastx', 'diamond', 'all') for tool in (options.tools or [])):
            FASTQ_to_FASTA(options)
        else:
            options.fastq_input = options.input
            options.fasta_input = None
    else:
        # Check input file exists
        if not os.path.exists(options.input):
            print(f"Error: Input file '{options.input}' not found", file=sys.stderr)
            sys.exit(1)
        options.fasta_input = options.input
        options.fastq_input = None

    if options.tools == ['all']:
        options.tools = ['blastn', 'blastx', 'diamond', 'bowtie2', 'bwa', 'minimap2']  #, 'hmmer_dna','hmmer_protein']

    # Initialise database variables to avoid UnboundLocalError
    resfinder_dbs = None
    card_dbs = None
    ncbi_dbs = None
    user_dbs = None
    # Load database paths from databases
    logger.info("Database selected: " + ', '.join(options.databases))
    if 'resfinder' in options.databases:
        resfinder_dbs = {tool: RESFINDER_DATABASES.get(tool) for tool in options.tools if RESFINDER_DATABASES.get(tool)}
    if 'card' in options.databases:
        card_dbs = {tool: CARD_DATABASES.get(tool) for tool in options.tools if CARD_DATABASES.get(tool)}
    if 'ncbi' in options.databases:
        ncbi_dbs = {tool: NCBI_DATABASES.get(tool) for tool in options.tools if NCBI_DATABASES.get(tool)}
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
        logger.error("Error: At least one database must be specified in databases.py or provided by the user",
                     file=sys.stderr)
        sys.exit(1)

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

    if hasattr(options, 'sensitivity') and options.sensitivity == 'default':
        # Use each tool's default sensitivity settings
        pass
    elif hasattr(options, 'sensitivity') and options.sensitivity == 'very-sensitive':
        # Example: set sensitivity for supported tools
        tool_sensitivity_params['bowtie2'] = {'sensitivity': '--very-sensitive-local'}
        tool_sensitivity_params['diamond'] = {'sensitivity': '--ultra-sensitive'}

    # Run Workflow
    workflow = AMRWorkflow(
        input_fasta=options.fasta_input,
        input_fastq=options.fastq_input,
        output_dir=options.output,
        databases=databases,
        # resfinder_dbs=resfinder_dbs,
        # card_dbs=card_dbs,
        threads=options.threads,
        tool_sensitivity_params=tool_sensitivity_params,
        #max_target_seqs=options.max_target_seqs,
        #evalue=options.evalue,
        detection_min_coverage=options.detection_min_coverage,
        detection_min_identity=options.detection_min_identity,
        detection_min_base_depth=options.detection_min_base_depth,
        detection_min_num_reads=options.detection_min_num_reads,
        query_min_coverage=options.query_min_coverage,
        run_dna=run_dna,
        run_protein=run_protein,
        sequence_type=options.sequence_type,
        report_fasta=options.report_fasta,
        no_cleanup=options.no_cleanup,
        verbose=options.verbose
    )

    ###
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

    # results = workflow.run_workflow(options)
    #
    # # Exit with error code if all tools failed
    # all_failed = True
    # for db_results in results.values():
    #     for success, _ in db_results.values():
    #         if success:
    #             all_failed = False
    #             break
    #     if not all_failed:
    #         break
    #
    # if all_failed:
    #     sys.exit(1)


if __name__ == "__main__":
    main()