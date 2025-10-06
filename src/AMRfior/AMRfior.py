import argparse
import sys
import os

try:
    from .constants import *
    from .databases import RESFINDER_DATABASES, CARD_DATABASES
    from .workflow import AMRWorkflow
    from .gene_stats import GeneStats
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    from databases import RESFINDER_DATABASES, CARD_DATABASES
    from workflow import AMRWorkflow
    from gene_stats import GeneStats


def main():
    parser = argparse.ArgumentParser(
        description='AMR Gene Detection Pipeline - Multi-tool alignment with detection matrices',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default tools
  python amr_pipeline.py -i reads.fasta -o results/

  # Select specific tools
  python amr_pipeline.py -i reads.fasta -o results/ \\
    --tools blastn diamond bowtie2

  # Custom thresholds with long reads
  python amr_pipeline.py -i nanopore.fasta -o results/ \\
    -t 16 --min-cov 90 --min-id 85 \\
    --minimap-preset map-ont --dna-only
        """
    )

    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA file with sequences to analyse')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for results')

    # Tool selection
    tool_group = parser.add_argument_group('Tool selection')
    tool_group.add_argument('--tools', nargs='+',
                            choices=['blastn', 'blastp', 'diamond', 'bowtie2', 'bwa', 'minimap2'], #, 'hmmer_dna', 'hmmer_protein'],
                            default=['blastn', 'blastp', 'diamond', 'bowtie2', 'bwa', 'minimap2'], #, 'hmmer_dna','hmmer_protein'],
                            help='Specify which tools to run (default: all)')


    gene_detection_group = parser.add_argument_group('Gene Detection Parameters')
    gene_detection_group.add_argument('--min-cov', '--min-coverage', type=float, default=80.0,
                              dest='min_coverage',
                              help='Minimum coverage threshold in percent (default: 80.0)')
    gene_detection_group.add_argument('--min-id', '--min-identity', type=float, default=80.0,
                              dest='min_identity',
                              help='Minimum identity threshold in percent (default: 80.0)')
    gene_detection_group.add_argument( '--max_target_seqs', dest='max_target_seqs', type=int, default=100,
                              help='Maximum number of "hits" to return per query sequence (default: 100)')


    # Mode selection
    mode_group = parser.add_argument_group('Mode Selection')
    mode_group.add_argument('--dna-only', action='store_true',
                            help='Run only DNA-based tools')
    mode_group.add_argument('--protein-only', action='store_true',
                            help='Run only protein-based tools')
    mode_group.add_argument('--sensitivity', type=str, default='default',
                            choices=['default', 'conservative', 'sensitive', 'very-sensitive'],
                            help='Preset sensitivity levels - default means each tool uses its own default settings and very-sensitive applies '
                                 'DIAMONDs --ultra-sensitive and Bowtie2s --very-sensitive-local presets')

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

    options = parser.parse_args()

    # Check input file exists
    if not os.path.exists(options.input):
        print(f"Error: Input file '{options.input}' not found", file=sys.stderr)
        sys.exit(1)

    # Load database paths from databases.py
    resfinder_dbs = {tool: RESFINDER_DATABASES.get(tool) for tool in options.tools if RESFINDER_DATABASES.get(tool)}
    card_dbs = {tool: CARD_DATABASES.get(tool) for tool in options.tools if CARD_DATABASES.get(tool)}

    if not resfinder_dbs and not card_dbs:
        print("Error: At least one database must be specified in databases.py", file=sys.stderr)
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
    #
    tool_sensitivity_params = {}

    if hasattr(options, 'sensitivity') and options.sensitivity == 'default':
        # Use each tool's default sensitivity settings
        pass
    elif hasattr(options, 'sensitivity') and options.sensitivity == 'very-sensitive':
        # Example: set sensitivity for supported tools
        tool_sensitivity_params['bowtie2'] = {'sensitivity': '--very-sensitive-local'}
        tool_sensitivity_params['diamond'] = {'sensitivity': '--ultra-sensitive'}







    # Run Workflow
    pipeline = AMRWorkflow(
        input_fasta=options.input,
        output_dir=options.output,
        resfinder_dbs=resfinder_dbs,
        card_dbs=card_dbs,
        threads=options.threads,
        tool_sensitivity_params=tool_sensitivity_params,
        #max_target_seqs=options.max_target_seqs,
        #evalue=options.evalue,
        min_coverage=options.min_coverage,
        min_identity=options.min_identity,
        run_dna=run_dna,
        run_protein=run_protein
    )

    results = pipeline.run_pipeline(options)

    # Exit with error code if all tools failed
    all_failed = True
    for db_results in results.values():
        for success, _ in db_results.values():
            if success:
                all_failed = False
                break
        if not all_failed:
            break

    if all_failed:
        sys.exit(1)


if __name__ == "__main__":
    main()