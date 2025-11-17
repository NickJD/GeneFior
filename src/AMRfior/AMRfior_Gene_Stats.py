import argparse
import sys
import os
from pathlib import Path
import logging
from datetime import datetime
from typing import Dict, List, Set, Tuple
from collections import defaultdict
import subprocess
import re

try:
    from .constants import *
except (ModuleNotFoundError, ImportError):
    sys.path.insert(0, str(Path(__file__).parent))
    from constants import *


class GenePosition:
    def __init__(self, pos: int):
        self.pos = pos
        self.depth = 0
        self.ref_base = None
        self.bases = defaultdict(int)  # {base: count}
        self.qualities = []

    def add_base(self, base: str, quality: int = None):
        # Add a base observation at this position.
        self.depth += 1
        self.bases[base] += 1
        if quality is not None:
            self.qualities.append(quality)

    def get_consensus(self):
        # Get most common base.
        if not self.bases:
            return None
        return max(self.bases.items(), key=lambda x: x[1])[0]

    def is_variant(self):
        # Check if position has variation.
        if self.ref_base is None or not self.bases:
            return False
        consensus = self.get_consensus()
        return consensus != self.ref_base

    def get_variant_freq(self):
        # Get frequency of variant bases.
        if not self.bases or self.depth == 0:
            return 0.0
        if self.ref_base is None:
            return 0.0
        variant_count = sum(count for base, count in self.bases.items() if base != self.ref_base)
        return (variant_count / self.depth) * 100


class GeneCoverage:
    # Store coverage information for a gene.

    def __init__(self, gene_name: str, gene_length: int, ref_seq: str = None):
        self.gene_name = gene_name
        self.gene_length = gene_length
        self.ref_seq = ref_seq
        self.positions = {i: GenePosition(i) for i in range(gene_length)}
        self.read_count = 0
        self.covered_positions = set()

        # Set reference bases if available
        if ref_seq:
            for i, base in enumerate(ref_seq):
                if i < gene_length:
                    self.positions[i].ref_base = base.upper()

    def add_alignment(self, start: int, end: int, aligned_seq: str = None, ref_start_in_seq: int = 0):
        # Add an alignment to coverage.
        self.read_count += 1
        for pos in range(start, end):
            if 0 <= pos < self.gene_length:
                #self.positions[pos].depth += 1
                self.covered_positions.add(pos)
                # Only increment depth if we're not tracking bases separately
                if aligned_seq is None:
                    self.positions[pos].depth += 1

                # Add base if sequence provided
                if aligned_seq and ref_start_in_seq >= 0:
                    seq_pos = pos - start + ref_start_in_seq
                    if 0 <= seq_pos < len(aligned_seq):
                        base = aligned_seq[seq_pos].upper()
                        self.positions[pos].add_base(base)

    def add_base_at_position(self, pos: int, base: str, quality: int = None):
        # Add a base observation at a specific position.
        if 0 <= pos < self.gene_length:
            self.positions[pos].add_base(base, quality)
            self.covered_positions.add(pos)

    def get_coverage_stats(self):
        # Calculate coverage statistics.
        covered = len(self.covered_positions)
        coverage_pct = (covered / self.gene_length * 100) if self.gene_length > 0 else 0

        depths = [pos.depth for pos in self.positions.values()]
        avg_depth = sum(depths) / len(depths) if depths else 0
        max_depth = max(depths) if depths else 0

        # Find gaps (uncovered regions)
        gaps = []
        in_gap = False
        gap_start = None

        for i in range(self.gene_length):
            if self.positions[i].depth == 0:
                if not in_gap:
                    gap_start = i
                    in_gap = True
            else:
                if in_gap:
                    gaps.append((gap_start, i - 1))
                    in_gap = False

        if in_gap:
            gaps.append((gap_start, self.gene_length - 1))

        # Count how many positions have base-level data
        positions_with_bases = sum(1 for pos in self.positions.values() if len(pos.bases) > 0)
        positions_with_ref = sum(1 for pos in self.positions.values() if pos.ref_base is not None)

        # Find variants
        variants = []
        min_depth = 3 ## Could User Param
        for i, pos in self.positions.items():
            if len(pos.bases) > 0 and pos.depth >= min_depth:
                if pos.ref_base and pos.is_variant():
                        variants.append({
                        'pos': i,
                        'ref': pos.ref_base,
                        'alt': pos.get_consensus(),
                        'freq': pos.get_variant_freq(),
                        'depth': pos.depth,
                        'bases': dict(pos.bases) # Include all base counts
                    })

        return {
            'covered_positions': covered,
            'coverage_percent': coverage_pct,
            'avg_depth': avg_depth,
            'max_depth': max_depth,
            'read_count': self.read_count,
            'gaps': gaps,
            'variants': variants,
            'positions_with_bases': positions_with_bases,
            'positions_with_ref': positions_with_ref
        }


class AMRVisualiser:
    # Generate coverage visualisations for AMR genes.

    def __init__(self, input_dir: str, output_dir: str, genes: List[str],
                 databases: List[str], tools: List[str], ref_fasta: str = None,
                 query_fasta: str = None):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.genes = set(genes) if genes else set()
        self.databases = databases
        self.tools = tools
        self.ref_fasta = ref_fasta
        self.query_fasta = query_fasta

        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.report_dir = self.output_dir / "coverage_reports"
        self.report_dir.mkdir(exist_ok=True)
        #self.plot_dir = self.output_dir / "coverage_plots"
        #self.plot_dir.mkdir(exist_ok=True)

        # Setup logging
        log_file = self.output_dir / f"gene-reports_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

        # Storage
        self.gene_coverages = defaultdict(lambda: defaultdict(dict))  # {database: {tool: {gene: GeneCoverage}}}
        self.gene_sequences = {}  # {gene: sequence} from reference
        self.query_sequences = {}  # Query sequences from FASTA

    def load_reference_sequences(self):
        # Load reference sequences from FASTA.
        if not self.ref_fasta or not Path(self.ref_fasta).exists():
            self.logger.warning("Reference FASTA not provided or not found - variant calling disabled")
            return

        self.logger.info(f"Loading reference sequences from {self.ref_fasta}")

        try:
            current_gene = None
            current_seq = []

            with open(self.ref_fasta, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_gene and current_seq:
                            self.gene_sequences[current_gene] = ''.join(current_seq)
                        current_gene = line[1:].split()[0]  # Get gene ID
                        current_seq = []
                    else:
                        current_seq.append(line)

                if current_gene and current_seq:
                    self.gene_sequences[current_gene] = ''.join(current_seq)

            self.logger.info(f"Loaded {len(self.gene_sequences)} reference sequences")
        except Exception as e:
            self.logger.error(f"Error loading reference FASTA: {e}")

    def load_query_sequences(self):
        # Load query sequences from FASTA.
        if not self.query_fasta or not Path(self.query_fasta).exists():
            self.logger.warning("Query FASTA not provided or not found - BLAST base-level data disabled")
            return

        self.logger.info(f"Loading query sequences from {self.query_fasta}")

        try:
            import gzip
            current_id = None
            current_seq = []
            count = 0

            # Handle gzipped files
            if str(self.query_fasta).endswith('.gz'):
                opener = gzip.open
                mode = 'rt'
            else:
                opener = open
                mode = 'r'

            with opener(self.query_fasta, mode) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id and current_seq:
                            self.query_sequences[current_id] = ''.join(current_seq)
                            count += 1
                        # Get full ID (first word after >)
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line.upper())

                if current_id and current_seq:
                    self.query_sequences[current_id] = ''.join(current_seq)
                    count += 1

            self.logger.info(f"Loaded {count} query sequences")
        except Exception as e:
            self.logger.error(f"Error loading query FASTA: {e}")
            import traceback
            self.logger.error(traceback.format_exc())

    def parse_blast_results(self, blast_file: Path, database: str, tool: str):
        # Parse BLAST format 6 output.
        self.logger.info(f"Parsing BLAST results: {database} - {tool}")

        try:
            alignment_count = 0
            alignments_with_seq = 0

            with open(blast_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 14:
                        continue

                    query_id = fields[0]
                    gene = fields[1]
                    if self.genes and gene not in self.genes:
                        continue

                    # BLAST coordinates
                    qstart = int(fields[6])  # Query start
                    qend = int(fields[7])  # Query end
                    sstart = int(fields[8])  # Subject start
                    send = int(fields[9])  # Subject end
                    slen = int(fields[13])  # Subject length

                    # Initialise coverage if needed
                    if gene not in self.gene_coverages[database][tool]:
                        ref_seq = self.gene_sequences.get(gene)
                        self.gene_coverages[database][tool][gene] = GeneCoverage(gene, slen, ref_seq)
                        if ref_seq:
                            self.logger.debug(
                                f"Initialised coverage for {gene} with reference sequence (length: {len(ref_seq)})")
                        else:
                            self.logger.debug(
                                f"Initialised coverage for {gene} WITHOUT reference sequence (length: {slen})")

                    # Add alignment (convert to 0-based)
                    start = min(sstart, send) - 1
                    end = max(sstart, send)

                    # Try to get query sequence for this alignment
                    query_seq = None
                    if query_id in self.query_sequences:
                        full_query = self.query_sequences[query_id]

                        # Extract aligned portion from query (1-based to 0-based)
                        q_start_idx = min(qstart, qend) - 1
                        q_end_idx = max(qstart, qend)

                        if 0 <= q_start_idx < len(full_query) and q_end_idx <= len(full_query):
                            query_seq = full_query[q_start_idx:q_end_idx]

                            # Handle reverse complement if needed
                            if sstart > send:  # Reverse strand alignment
                                query_seq = self._reverse_complement(query_seq)

                            alignments_with_seq += 1

                            # Add alignment with sequence data
                            # Map query bases to subject positions
                            for i, base in enumerate(query_seq):
                                subject_pos = start + i
                                if 0 <= subject_pos < slen:
                                    self.gene_coverages[database][tool][gene].add_base_at_position(
                                        subject_pos, base
                                    )

                    # Always add the alignment for coverage tracking
                    self.gene_coverages[database][tool][gene].add_alignment(start, end,
                                                                            aligned_seq=None)  # Depth already tracked above
                    alignment_count += 1

                self.logger.info(f"Processed {alignment_count} alignments from {blast_file}")
                if self.query_sequences:
                    self.logger.info(
                        f"  {alignments_with_seq} alignments with sequence data ({alignments_with_seq / alignment_count * 100:.1f}%)")
                else:
                    self.logger.warning(
                        f"  Query sequences not loaded - variant calling not available for {database}/{tool}")


        except Exception as e:
            self.logger.error(f"Error parsing BLAST file {blast_file}: {e}")
            import traceback
            self.logger.error(traceback.format_exc())

    def _reverse_complement(self, seq: str) -> str:
        # Return reverse complement of a DNA sequence.
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq))

    def parse_bam_file(self, bam_file: Path, database: str, tool: str):
        # Parse BAM file for coverage with base-level tracking.
        self.logger.info(f"Parsing BAM file: {database} - {tool}")

        try:
            proc = subprocess.Popen(['samtools', 'view', '-h', str(bam_file)],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            gene_lengths = {}
            cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')
            alignment_count = 0

            for line in proc.stdout:
                # Parse header for gene lengths
                if line.startswith('@SQ'):
                    parts = line.strip().split('\t')
                    gene_name = gene_len = None
                    for p in parts:
                        if p.startswith('SN:'):
                            gene_name = p.split(':', 1)[1]
                        if p.startswith('LN:'):
                            gene_len = int(p.split(':', 1)[1])
                    if gene_name and gene_len:
                        gene_lengths[gene_name] = gene_len
                    continue

                if line.startswith('@'):
                    continue

                fields = line.rstrip('\n').split('\t')
                if len(fields) < 11:
                    continue

                flag = int(fields[1])
                if flag & 0x4:  # Unmapped
                    continue

                gene = fields[2]
                if self.genes and gene not in self.genes:
                    continue

                ref_start = int(fields[3]) - 1  # 0-based
                cigar = fields[5]
                seq = fields[9]
                qual = fields[10] if len(fields) > 10 else None

                # Initialise coverage if needed
                gene_len = gene_lengths.get(gene, 0)
                if gene not in self.gene_coverages[database][tool]:
                    ref_seq = self.gene_sequences.get(gene)
                    self.gene_coverages[database][tool][gene] = GeneCoverage(gene, gene_len, ref_seq)
                    self.logger.debug(
                        f"Initialized coverage for {gene} (length: {gene_len}, ref_seq: {'yes' if ref_seq else 'no'})")

                # Parse CIGAR and track coverage with base information
                ref_pos = ref_start
                seq_pos = 0

                for count_str, op in cigar_re.findall(cigar):
                    length = int(count_str)

                    if op in ('M', '=', 'X'):
                        # Add coverage with sequence information
                        for i in range(length):
                            current_ref_pos = ref_pos + i
                            current_seq_pos = seq_pos + i

                            if 0 <= current_ref_pos < gene_len and current_seq_pos < len(seq):
                                base = seq[current_seq_pos].upper()
                                quality = ord(qual[current_seq_pos]) - 33 if qual and current_seq_pos < len(
                                    qual) else None

                                # Add base observation
                                self.gene_coverages[database][tool][gene].add_base_at_position(
                                    current_ref_pos, base, quality
                                )

                        ref_pos += length
                        seq_pos += length

                    elif op == 'I':
                        seq_pos += length

                    elif op in ('D', 'N'):
                        ref_pos += length

                    elif op == 'S':
                        seq_pos += length

                self.gene_coverages[database][tool][gene].read_count += 1
                alignment_count += 1

            proc.stdout.close()
            proc.wait()

            self.logger.info(f"Processed {alignment_count} alignments from {bam_file}")

        except Exception as e:
            self.logger.error(f"Error parsing BAM file {bam_file}: {e}")
            import traceback
            self.logger.error(traceback.format_exc())

    def generate_text_report(self, gene: str, database: str, tool: str, coverage: GeneCoverage):
        # Generate detailed text report for a gene.
        stats = coverage.get_coverage_stats()


        safe_gene = gene.replace('|', '_').replace('/', '_').replace(':','_').replace('-','_')
        if gene != safe_gene:
            self.logger.info(f"Renaming gene for report file:  {gene}  to  {safe_gene}")

        report_file = self.report_dir / f"{database}_{tool}_{safe_gene}_coverage.txt"

        if tool in ('blastx', 'diamond'): # +3 gene length for AA
            marker = 'aa'
        else:
            marker = 'bp'

        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"AMRfíor Coverage Report\n")
            f.write(f"Gene: {gene}\n")
            f.write(f"Database: {database}\n")
            f.write(f"Tool: {tool}\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 80 + "\n\n")

            # Summary statistics
            f.write("COVERAGE SUMMARY\n")
            f.write("-" * 80 + "\n")
            f.write(f"Gene length:           {coverage.gene_length:,} {marker}\n")
            f.write(f"Covered positions:     {stats['covered_positions']:,} {marker} ({stats['coverage_percent']:.2f}%)\n")
            f.write(f"Total reads mapped:    {stats['read_count']:,}\n")
            f.write(f"Average depth:         {stats['avg_depth']:.2f}×\n")
            f.write(f"Maximum depth:         {stats['max_depth']}×\n")
            f.write(f"Number of gaps:        {len(stats['gaps'])}\n")
            f.write(f"Number of variants:    {len(stats['variants'])}\n")

            # Diagnostic info
            f.write(f"\nVARIANT CALLING STATUS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Positions with base data:     {stats['positions_with_bases']:,}\n")
            f.write(f"Positions with reference:     {stats['positions_with_ref']:,}\n")
            f.write(f"Number of variants detected:  {len(stats['variants'])}\n")

            if stats['positions_with_bases'] == 0:
                f.write(f"Note: {tool} output doesn't contain base-level sequence data\n")
                f.write(f"      Variant calling requires BAM/SAM format (bowtie2, bwa, minimap2)\n")
            elif stats['positions_with_ref'] == 0:
                f.write(f"Note: Reference sequence not provided - variant calling disabled\n")
                f.write(f"      Use --ref-fasta to enable variant detection\n")
            f.write("\n")

            # Gap regions
            if stats['gaps']:
                f.write("UNCOVERED REGIONS (GAPS)\n")
                f.write("-" * 80 + "\n")
                f.write(f"{'Start':<10} {'End':<10} {'Length':<10} {'% of Gene':<12}\n")
                f.write("-" * 80 + "\n")

                for gap_start, gap_end in stats['gaps']:
                    gap_len = gap_end - gap_start + 1
                    gap_pct = (gap_len / coverage.gene_length) * 100
                    f.write(f"{gap_start + 1:<10} {gap_end + 1:<10} {gap_len:<10} {gap_pct:<12.2f}\n")
                f.write("\n")

            # Variants
            if stats['variants']:
                f.write("SEQUENCE VARIANTS\n")
                f.write("-" * 80 + "\n")
                f.write(f"{'Position':<10} {'Ref':<5} {'Alt':<5} {'Depth':<8} {'Var %':<10} {'Base Counts':<20}\n")
                f.write("-" * 80 + "\n")

                for var in sorted(stats['variants'], key=lambda x: x['pos']):
                    base_counts = ','.join([f"{b}:{c}" for b, c in sorted(var['bases'].items())])
                    f.write(f"{var['pos'] + 1:<10} {var['ref']:<5} {var['alt']:<5} "
                            f"{var['depth']:<8} {var['freq']:<10.2f}\n")
                f.write("\n")
            # else:
            #     f.write("SEQUENCE VARIANTS\n")
            #     f.write("-" * 80 + "\n")
            #     if coverage.ref_seq:
            #         f.write("No variants detected (minimum depth: 3×)\n\n")
            #     else:
            #         f.write("Variant calling disabled (no reference sequence provided)\n\n")

            # Coverage visualisation (ASCII art)
            f.write("COVERAGE VISUALISATION\n")
            f.write("-" * 80 + "\n")

            # Bin positions for visualisation
            bin_size = max(1, coverage.gene_length // 100)
            bins = []

            for i in range(0, coverage.gene_length, bin_size):
                bin_end = min(i + bin_size, coverage.gene_length)
                bin_depths = [coverage.positions[j].depth for j in range(i, bin_end)]
                avg_bin_depth = sum(bin_depths) / len(bin_depths) if bin_depths else 0
                bins.append(avg_bin_depth)

            # Normalise for display
            max_display_depth = 50
            if stats['max_depth'] > 0:
                bins_normalised = [min(int((d / stats['max_depth']) * max_display_depth), max_display_depth)
                                   for d in bins]
            else:
                bins_normalised = [0] * len(bins)

            # Draw ASCII histogram
            f.write(f"Position ({marker}):  1 {' ' * 30} {coverage.gene_length // 2} {' ' * 28} {coverage.gene_length}\n")
            f.write(f"Max depth: {stats['max_depth']}×\n")
            f.write("-" * 80 + "\n")

            for row in range(max_display_depth, 0, -5):
                line = f"{row:3}× |"
                for val in bins_normalised:
                    if val >= row:
                        line += "█"
                    else:
                        line += " "
                f.write(line + "\n")

            f.write("     +" + "-" * len(bins) + "\n")

        self.logger.info(f"Generated report: {report_file}")

    def generate_comparison_report(self, gene: str, database: str):
        # Generate comparison report across all tools for a gene.
        safe_gene = gene.replace('|', '_').replace('/', '_').replace(':', '_').replace('-', '_')
        report_file = self.report_dir / f"{database}_{safe_gene}_comparison.txt"

        tools_with_data = [tool for tool in self.tools
                           if gene in self.gene_coverages[database].get(tool, {})]

        if not tools_with_data:
            return

        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"AMRfíor Multi-Tool Comparison\n")
            f.write(f"Gene: {gene}\n")
            f.write(f"Database: {database}\n")
            f.write("=" * 80 + "\n\n")

            # Comparison table
            f.write("COVERAGE COMPARISON\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Tool':<15} {'Coverage %':<12} {'Avg Depth':<12} {'Reads':<10} {'Gaps':<8} {'Variants':<10}\n")
            f.write("-" * 80 + "\n")

            for tool in tools_with_data:
                coverage = self.gene_coverages[database][tool][gene]
                stats = coverage.get_coverage_stats()

                f.write(f"{tool:<15} {stats['coverage_percent']:<12.2f} "
                        f"{stats['avg_depth']:<12.2f} {stats['read_count']:<10} "
                        f"{len(stats['gaps']):<8} {len(stats['variants']):<10}\n")

            f.write("\n")

            # Agreement/disagreement regions
            f.write("TOOL AGREEMENT ANALYSIS\n")
            f.write("-" * 80 + "\n")

            # Find positions covered by all tools vs some tools
            all_positions = set(range(self.gene_coverages[database][tools_with_data[0]][gene].gene_length))
            covered_by_all = all_positions.copy()
            covered_by_any = set()

            for tool in tools_with_data:
                tool_covered = self.gene_coverages[database][tool][gene].covered_positions
                covered_by_all &= tool_covered
                covered_by_any |= tool_covered

            consensus_pct = (len(covered_by_all) / len(all_positions)) * 100 if all_positions else 0
            any_covered_pct = (len(covered_by_any) / len(all_positions)) * 100 if all_positions else 0

            f.write(f"Positions covered by ALL tools:  {len(covered_by_all):,} ({consensus_pct:.2f}%)\n")
            f.write(f"Positions covered by ANY tool:   {len(covered_by_any):,} ({any_covered_pct:.2f}%)\n")
            f.write(f"Tool-specific regions:           {len(covered_by_any - covered_by_all):,}\n")

        self.logger.info(f"Generated comparison report: {report_file}")

    def discover_files(self):
        # Discover alignment files in input directory.
        self.logger.info("Discovering alignment files...")

        raw_dir = self.input_dir / "raw_outputs"
        if not raw_dir.exists():
            self.logger.error(f"Raw outputs directory not found: {raw_dir}")
            return False

        found_files = []

        for database in self.databases:
            for tool in self.tools:
                # Try BLAST format first
                if tool in ['blastn', 'blastx', 'diamond']:
                    pattern = f"{database}_{tool}_results.tsv"
                    files = list(raw_dir.glob(pattern))
                    if files:
                        found_files.append((database, tool, 'blast', files[0]))
                        self.logger.info(f"Found BLAST file: {files[0]}")

                # Try BAM format
                elif tool in ['bowtie2', 'bwa', 'minimap2']:
                    pattern = f"{database}_{tool}_results_sorted.bam"
                    files = list(raw_dir.glob(pattern))
                    if files:
                        found_files.append((database, tool, 'bam', files[0]))
                        self.logger.info(f"Found BAM file: {files[0]}")

        if not found_files:
            self.logger.error("No alignment files found!")
            return False

        self.logger.info(f"Found {len(found_files)} alignment files")
        self.found_files = found_files
        return True

    def run(self):
        # Main execution.
        self.logger.info("=" * 70)
        self.logger.info(f"AMRfíor-Gene-Report {AMRFIOR_VERSION}")
        self.logger.info("=" * 70)
        self.logger.info(f"Input directory: {self.input_dir}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Genes to visualise: {len(self.genes) if self.genes else 'ALL'}")
        self.logger.info(f"Databases: {', '.join(self.databases)}")
        self.logger.info(f"Tools: {', '.join(self.tools)}")
        self.logger.info("=" * 70)


        ############ NOT IMPLEMENTED
        # Load reference sequences
        # if self.ref_fasta:
        #     self.load_reference_sequences()

        # Load query sequences
        # if self.query_fasta:
        #     self.load_query_sequences()
        # else:
        #     self.logger.warning("No query FASTA provided - BLAST variant calling will be limited")

        # Discover files
        if not self.discover_files():
            return False

        # Process files
        for database, tool, file_type, filepath in self.found_files:
            if file_type == 'blast':
                self.parse_blast_results(filepath, database, tool)
            elif file_type == 'bam':
                self.parse_bam_file(filepath, database, tool)

        # Generate reports
        self.logger.info("\nGenerating coverage reports...")

        report_count = 0
        for database in self.databases:
            for tool in self.tools:
                if tool not in self.gene_coverages[database]:
                    continue

                for gene, coverage in self.gene_coverages[database][tool].items():
                    self.generate_text_report(gene, database, tool, coverage)
                    report_count += 1

            # Generate comparison reports
            all_genes = set()
            for tool in self.tools:
                if tool in self.gene_coverages[database]:
                    all_genes.update(self.gene_coverages[database][tool].keys())

            for gene in all_genes:
                self.generate_comparison_report(gene, database)

        # Summary
        self.logger.info("\n" + "=" * 70)
        self.logger.info("GENE REPORT COMPLETE")
        self.logger.info("=" * 70)
        total_reports = len(list(self.report_dir.glob("*.txt")))
        self.logger.info(f"Generated {total_reports} coverage reports")
        self.logger.info(f"Reports saved to: {self.report_dir}")
        self.logger.info("=" * 70)

        return True


def main():
    parser = argparse.ArgumentParser(
        description='AMRfíor ' + AMRFIOR_VERSION + ' - AMRfíor-Gene-Stats: Generate detailed coverage visualisations for AMR genes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Visualise specific genes (FULL NAMES) from all tools
  AMRfior-gene-stats -i results/ -o vis/ \\
    -g "sul1_2_U12338,tet(W)|ARO:3000194" \\
    --databases resfinder card \\
    --tools diamond bowtie2 bwa

  # Visualise from gene (FULL NAMES) list file with reference
  AMRfior-gene-stats -i results/ -o vis/ \\
    -g genes_of_interest.txt \\
    --databases resfinder \\
    --tools blastn diamond 

        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input directory containing AMRfíor results')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for visualisation reports')
    parser.add_argument('-g', '--genes', required=False,
                        help='Comma-separated gene names (FULL NAMES) or path to file with gene names (one per line)')
    parser.add_argument('--databases', nargs='+', required=True,
                        choices=['resfinder', 'card', 'ncbi'],
                        help='Database(s) to interrogate')
    parser.add_argument('--tools', nargs='+', required=True,
                        choices=['blastn', 'blastx', 'diamond', 'bowtie2', 'bwa', 'minimap2', 'all'],
                        help='Tool(s) to interrogate')
    parser.add_argument('--ref-fasta',
                        help='NOT IMPLEMENTED YET - Reference FASTA file for variant calling (optional)')
    parser.add_argument('--query-fasta',
                        help='NOT IMPLEMENTED YET - Query FASTA file (your input reads) for BLAST base-level analysis (optional)')


    options = parser.parse_args()

    # Check input directory
    if not os.path.exists(options.input):
        print(f"Error: Input directory '{options.input}' not found", file=sys.stderr)
        sys.exit(1)

        # Tool selection
    if options.tools == ['all']:
        options.tools = ['blastn', 'blastx', 'diamond', 'bowtie2', 'bwa', 'minimap2']  # , 'hmmer_dna','hmmer_protein']

    # Parse genes
    genes = []
    if options.genes:
        genes_input = Path(options.genes)
        if genes_input.exists() and genes_input.is_file():
            # Read from file
            with open(genes_input, 'r') as f:
                genes = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        else:
            # Parse as comma-separated list
            genes = [g.strip() for g in options.genes.split(',') if g.strip()]

    # Run visualiser
    visualiser = AMRVisualiser(
        input_dir=options.input,
        output_dir=options.output,
        genes=genes,
        databases=options.databases,
        tools=options.tools,
        ref_fasta=options.ref_fasta,
        query_fasta=options.query_fasta
    )

    success = visualiser.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
