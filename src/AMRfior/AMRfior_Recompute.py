import argparse
import sys
import os
from pathlib import Path
import logging
from datetime import datetime
import shutil
from collections import defaultdict

try:
    from .gene_stats import GeneStats
    from .constants import AMRFIOR_VERSION
    from .workflow import AMRWorkflow
except (ModuleNotFoundError, ImportError):
    #sys.path.insert(0, str(Path(__file__).parent))
    from gene_stats import GeneStats
    from constants import AMRFIOR_VERSION
    from workflow import AMRWorkflow



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

    # def parse_blast_results(self, output_file: Path, database: str, tool_name: str) -> Set[str]:
    #     # Parse BLAST/DIAMOND tabular output.
    #     detected_genes = set()
    #     gene_lengths = {}
    #     gene_reads = defaultdict(lambda: {'passing': [], 'all': []})
    #
    #     if not output_file.exists():
    #         self.logger.warning(f"File not found: {output_file}")
    #         return detected_genes
    #
    #     self.logger.info(f"Processing {database} - {tool_name}...")
    #
    #     try:
    #         with open(output_file, 'r') as f:
    #             for line in f:
    #                 if line.startswith('#'):
    #                     continue
    #                 fields = line.strip().split('\t')
    #                 if len(fields) < 14:
    #                     continue
    #
    #                 read_name = fields[0]
    #                 gene = fields[1]
    #                 identity = float(fields[2])
    #                 qstart = int(fields[6])
    #                 qend = int(fields[7])
    #                 sstart = int(fields[8])
    #                 send = int(fields[9])
    #                 qlen = int(fields[12])
    #                 slen = int(fields[13])
    #
    #                 # Store gene length
    #                 gene_lengths[gene] = max(gene_lengths.get(gene, 0), slen)
    #
    #                 # Calculate query coverage
    #                 query_coverage = ((abs(qend - qstart) + 1) / qlen) * 100 if qlen else 0
    #
    #                 # Track all reads
    #                 gene_reads[gene]['all'].append(read_name)
    #
    #                 # Apply thresholds
    #                 if identity >= self.detection_min_identity and query_coverage >= self.query_min_coverage:
    #                     if gene not in self.gene_stats[database][tool_name]:
    #                         self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)
    #
    #                     self.gene_stats[database][tool_name][gene].add_hit(
    #                         sstart, send, identity, gene_lengths[gene]
    #                     )
    #                     gene_reads[gene]['passing'].append(read_name)
    #
    #     except Exception as e:
    #         self.logger.error(f"Error parsing {output_file}: {e}")
    #         return detected_genes
    #
    #     # Finalise and detect
    #     for gene in self.gene_stats[database][tool_name]:
    #         stats = self.gene_stats[database][tool_name][gene]
    #         stats.finalise()
    #
    #         if (stats.gene_coverage >= self.detection_min_coverage and
    #                 stats.base_depth >= self.detection_min_base_depth and
    #                 stats.num_sequences >= self.detection_min_num_reads):
    #             detected_genes.add(gene)
    #             self.detections[database][gene][tool_name] = True
    #
    #     self.logger.info(f"  {len(detected_genes)} genes detected")
    #     return detected_genes, gene_reads

    # def parse_bam_results(self, bam_file: Path, database: str, tool_name: str) -> Set[str]:
    #     # Parse BAM file using samtools view.
    #     detected_genes = set()
    #     gene_lengths = {}
    #     gene_reads = defaultdict(lambda: {'passing': [], 'all': []})
    #
    #     if not bam_file.exists():
    #         self.logger.warning(f"File not found: {bam_file}")
    #         return detected_genes
    #
    #     self.logger.info(f"Processing {database} - {tool_name}...")
    #
    #     try:
    #         import subprocess
    #         import re
    #
    #         proc = subprocess.Popen(['samtools', 'view', '-h', str(bam_file)],
    #                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    #
    #         cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')
    #
    #         for line in proc.stdout:
    #             if line.startswith('@SQ'):
    #                 parts = line.strip().split('\t')
    #                 sn = ln = None
    #                 for p in parts:
    #                     if p.startswith('SN:'):
    #                         sn = p.split(':', 1)[1]
    #                     if p.startswith('LN:'):
    #                         ln = int(p.split(':', 1)[1])
    #                 if sn and ln:
    #                     gene_lengths[sn] = ln
    #                 continue
    #             if line.startswith('@'):
    #                 continue
    #
    #             fields = line.rstrip('\n').split('\t')
    #             if len(fields) < 11:
    #                 continue
    #
    #             read_name = fields[0]
    #             flag = int(fields[1])
    #             if flag & 0x4:  # Unmapped
    #                 continue
    #
    #             gene = fields[2]
    #             ref_start = int(fields[3]) - 1
    #             cigar = fields[5]
    #             seq = fields[9]
    #
    #             gene_len = gene_lengths.get(gene, 0)
    #             gene_reads[gene]['all'].append(read_name)
    #
    #             # Get NM tag
    #             nm = 0
    #             for opt in fields[11:]:
    #                 if opt.startswith('NM:i:'):
    #                     nm = int(opt.split(':')[-1])
    #                     break
    #
    #             # Parse CIGAR
    #             ref_pos = ref_start
    #             aligned_positions = set()
    #             alignment_length = 0
    #
    #             for count_str, op in cigar_re.findall(cigar):
    #                 length = int(count_str)
    #                 if op in ('M', '=', 'X'):
    #                     aligned_positions.update(range(ref_pos, ref_pos + length))
    #                     ref_pos += length
    #                     alignment_length += length
    #                 elif op == 'I':  # insertion to reference
    #                     alignment_length += length
    #                 elif op == 'D':  # deletion from reference
    #                     ref_pos += length
    #                     alignment_length += length
    #                 elif op == 'N':
    #                     ref_pos += length
    #                 elif op in ('S', 'H'):
    #                     # soft/hard clip - do not consume reference (H doesn't appear in SEQ)
    #                     pass
    #
    #             identity = ((alignment_length - nm) / alignment_length * 100) if alignment_length > 0 else 0
    #             query_length = len(seq) if seq and seq != '*' else 0
    #             query_coverage = (len(aligned_positions) / query_length * 100) if query_length > 0 else 0
    #
    #             if identity >= self.detection_min_identity and query_coverage >= self.query_min_coverage:
    #                 if gene not in self.gene_stats[database][tool_name]:
    #                     self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)
    #
    #                 self.gene_stats[database][tool_name][gene].add_positions(
    #                     aligned_positions, identity, gene_len
    #                 )
    #                 gene_reads[gene]['passing'].append(read_name)
    #
    #         proc.stdout.close()
    #         proc.wait()
    #
    #     except Exception as e:
    #         self.logger.error(f"Error parsing {bam_file}: {e}")
    #         return detected_genes
    #
    #     # Finalise and detect
    #     for gene in self.gene_stats[database][tool_name]:
    #         stats = self.gene_stats[database][tool_name][gene]
    #         stats.finalise()
    #
    #         if (stats.gene_coverage >= self.detection_min_coverage and
    #                 stats.base_depth >= self.detection_min_base_depth and
    #                 stats.num_sequences >= self.detection_min_num_reads):
    #             detected_genes.add(gene)
    #             self.detections[database][gene][tool_name] = True
    #
    #     self.logger.info(f"  {len(detected_genes)} genes detected")
    #     return detected_genes, gene_reads
    #
    # def write_tool_stats(self, database: str, tool_name: str, gene_reads: dict = None):
    #     # Write detailed statistics for a specific tool to TSV.
    #     stats_file = self.stats_dir / f"{database}_{tool_name}_stats.tsv"
    #
    #     gene_stats = self.gene_stats[database][tool_name]
    #     if not gene_stats:
    #         self.logger.warning(f"No statistics to write for {database} - {tool_name}")
    #         return
    #
    #     with open(stats_file, 'w', newline='') as f:
    #         writer = csv.writer(f, delimiter='\t')
    #
    #         header = ['Gene', 'Gene_Length', 'Num_Sequences_Mapped',
    #                   'Num_Sequences_Passing_Thresholds', 'Gene_Coverage',
    #                   'Base_Coverage', 'Base_Coverage_Hit', 'Avg_Identity', 'Detected']
    #         writer.writerow(header)
    #
    #         genes = sorted(gene_stats.keys())
    #
    #         for gene in genes:
    #             stats = gene_stats[gene]
    #             try:
    #                 detected = self.detections[database][gene][tool_name]
    #             except (KeyError, TypeError):
    #                 detected = False
    #
    #             row = [
    #                 gene,
    #                 stats.gene_length,
    #                 len(gene_reads.get(gene, {}).get('all', [])),
    #                 len(gene_reads.get(gene, {}).get('passing', [])),
    #                 f"{stats.gene_coverage:.2f}",
    #                 f"{stats.base_depth:.2f}",
    #                 f"{stats.base_depth_hit:.2f}",
    #                 f"{stats.avg_identity:.2f}",
    #                 '1' if detected else '0'
    #             ]
    #             writer.writerow(row)
    #
    #     self.logger.info(f"  Stats file: {stats_file}")
    #
    # def generate_detection_matrix(self, database: str):
    #     # Generate TSV matrix of gene detections across tools.
    #     output_file = self.output_dir / f"{database}_detection_matrix.tsv"
    #
    #     all_tools = set()
    #     for gene_detections in self.detections[database].values():
    #         all_tools.update(gene_detections.keys())
    #
    #     if not all_tools:
    #         self.logger.info(f"No detections found for {database} - No matrix generated.")
    #         return
    #
    #     all_tools = sorted(all_tools)
    #
    #     with open(output_file, 'w', newline='') as f:
    #         writer = csv.writer(f, delimiter='\t')
    #
    #         header = ['Gene'] + all_tools + ['Total_Detections']
    #         writer.writerow(header)
    #
    #         if database == 'card':
    #             def get_last_segment(gene_name):
    #                 return gene_name.split('|')[-1] if '|' in gene_name else gene_name
    #
    #             genes = [
    #                 gene for gene in sorted(self.detections[database].keys(), key=get_last_segment)
    #                 if any(self.detections[database][gene][tool] for tool in all_tools)
    #             ]
    #         else:
    #             genes = [
    #                 gene for gene in sorted(self.detections[database].keys())
    #                 if any(self.detections[database][gene][tool] for tool in all_tools)
    #             ]
    #
    #         for gene in genes:
    #             row = [gene]
    #             detections = self.detections[database][gene]
    #
    #             for tool in all_tools:
    #                 row.append('1' if detections[tool] else '0')
    #
    #             total = sum(1 for tool in all_tools if detections[tool])
    #             row.append(str(total))
    #
    #             writer.writerow(row)
    #
    #     self.logger.info(f"Generated detection matrix: {output_file}")
    #     self.logger.info(f"  Total genes detected: {len(genes)}")
    #     self.logger.info(f"  Tools used: {len(all_tools)}")

def run(options, workflow, logger):

    # Storage structures
    gene_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(GeneStats)))
    detections = defaultdict(lambda: defaultdict(lambda: defaultdict(bool)))
    databases_found = set()
    tools_found = defaultdict(set)

    logger.info("=" * 70)
    logger.info(f"AMRfíor-Recompute {AMRFIOR_VERSION}")
    logger.info("=" * 70)
    logger.info(f"Input directory: {options.input}")
    logger.info(f"Output directory: {options.output}")
    logger.info(f"Detection thresholds:")
    logger.info(f"  Query min coverage: {options.query_min_coverage}%")
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
            if tool in ['blastn', 'blastx', 'diamond']:
                #tool_name = {'blast': 'BLASTn', 'blastx': 'BLASTx', 'diamond': 'DIAMOND'}[tool]
                detected, gene_reads = workflow.parse_blast_results(file_path, database, tool)
                workflow.write_tool_stats(database, tool, gene_reads)
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
        description='AMRfíor ' + AMRFIOR_VERSION + ' - AMRfíor-Recompute: Recalculate detection statistics from existing sequence search outputs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Recompute with different thresholds
  AMRfior-recompute -i original_results/ -o recomputed_90_90/ \\
    --d-min-cov 90 --d-min-id 90

  # More stringent depth requirement
  AMRfior-recompute -i original_results/ -o high_depth/ \\
    --d-min-base-depth 5.0 --d-min-reads 10

        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input directory containing AMRfíor results (with raw_outputs/ subdirectory)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for recomputed results')
    parser.add_argument('--tools', nargs='+',
                            choices=['blastn', 'blastx', 'diamond', 'bowtie2', 'bwa', 'minimap2', 'all'], #, 'hmmer_dna', 'hmmer_protein'],
                            default=['all'], #, 'hmmer_dna','hmmer_protein'],
                            help='Specify which tools to recompute - "all" will recompute for all detected tools (default: all)')

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

    # Output selection
    output_group = parser.add_argument_group('Output Parameterts')
    output_group.add_argument('--report-fasta',
                            choices=['None', 'all', 'detected', 'detected-all'], #, 'hmmer_dna', 'hmmer_protein'],
                            default=[None], #, 'hmmer_dna','hmmer_protein'],
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
                            version='AMRfíor ' + AMRFIOR_VERSION,
                            help='Show program version and exit')

    options = parser.parse_args()

    ## Setup logging
    start_time = datetime.now()
    from pathlib import Path
    log_file = Path(options.output) / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = logging.getLogger('AMRfíor-Recompute')
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
        options.tools = ['blastn', 'blastx', 'diamond', 'bowtie2', 'bwa', 'minimap2']  # , 'hmmer_dna','hmmer_protein']
    if options.report_fasta != ['None'] and any(tool in options.tools for tool in ['blastx', 'blastn', 'diamond']):
        if options.query_fasta is None:
            logger.error("Error: --query-fasta must be provided when --report-fasta is used with blast/diamond outputs")
            sys.exit(1)

    # Create output directory
    # options.output.mkdir(parents=True, exist_ok=True)
    # stats_dir = options.output / "recomputed_stats"
    # stats_dir.mkdir(exist_ok=True)

    workflow = AMRWorkflow(
        input_fasta=options.query_fasta,
        input_fastq=None,
        output_dir=options.output,
        databases={},  # to be set based on user input or discovered files
        threads=options.threads if hasattr(options, 'threads') else 4,
        tool_sensitivity_params={},  # to be set if user provides
        detection_min_coverage=options.detection_min_coverage,
        detection_min_identity=options.detection_min_identity,
        detection_min_base_depth=options.detection_min_base_depth,
        detection_min_num_reads=options.detection_min_num_reads,
        query_min_coverage=options.query_min_coverage,
        run_dna=True,
        run_protein=True,
        sequence_type='Single-FASTA',
        report_fasta=options.report_fasta,
        no_cleanup=False,
        verbose=False,
        logger = logger
    )



    run(options, workflow, logger)

    #shutil.rmtree(workflow.raw_dir)

    for handler in logger.handlers:
        handler.flush()
    for handler in logger.handlers:
        handler.close()
    logging.shutdown()

    end_time = datetime.now()
    elapsed = end_time - start_time
    logger.info(f"AMRfíor-Recompute completed in {elapsed}.")



if __name__ == "__main__":
    main()
