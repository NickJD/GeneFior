import argparse
import sys
import os
from pathlib import Path
import logging
from datetime import datetime
from typing import List
from collections import defaultdict
import subprocess
import re
import traceback
import json

try:
    from .constants import *
except (ModuleNotFoundError, ImportError):
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


class GeneVisualiser:
    # Generate coverage visualisations for searched genes.

    def __init__(self, input_dir: str, output_dir: str, genes: List[str],
                 databases: List[str], tools: List[str], ref_fasta: str = None,
                 query_fasta: str = None, plot_per_tool: bool = False):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.genes = set(genes) if genes else set()
        self.databases = databases
        self.tools = tools
        self.ref_fasta = ref_fasta
        self.query_fasta = query_fasta

        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.report_dir = self.output_dir / "gene_reports"
        self.report_dir.mkdir(exist_ok=True)
        #self.plot_dir = self.output_dir / "coverage_plots"
        self.plot_dir = self.output_dir / "gene_plots"
        self.plot_dir.mkdir(exist_ok=True)

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
        # Detected genes cache per database (read from <input_dir>/<database>_detection_matrix.tsv)
        self.detected_genes_by_db = {}
        # Default behavior: only report genes that were marked as detected in the detection matrix
        self.report_only_detected = True
        # Query filtering thresholds - will be loaded from run_parameters.json if present
        self.query_min_coverage = None
        self.query_min_identity = None
        # Try to enable plotting (matplotlib). If not available, continue without plots.
        # By default we DO NOT generate per-tool individual coverage PNGs; only combined comparison plots
        # will be produced. Use the CLI flag --plot-per-tool to enable per-tool PNG generation.
        self.plot_enabled = False
        # Controlled via constructor/CLI: if True generate individual per-tool PNGs in addition to comparison plots
        self.generate_individual_plots = bool(plot_per_tool)
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            self.plt = plt
            self.plot_enabled = True
            logging.getLogger(__name__).info('matplotlib available - plotting enabled (individual per-tool plots disabled by default)')
        except Exception:
            # matplotlib may not be installed in minimal environments
            logging.getLogger(__name__).warning('matplotlib not available - plots will be disabled')

        # Attempt to load run parameters from the input directory (if present)
        try:
            self.load_run_parameters()
        except Exception:
            self.logger.debug('No run parameters file found or failed to load; continuing without applying upstream query filters')

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

                    # If this is a protein alignment (blastx/diamond) the subject
                    # coordinates and length are in amino-acids. Convert to nucleotide
                    # coordinates so coverage plots align with nucleotide-based tools
                    # (e.g. bwa/bowtie2). Conversion: aa position -> nucleotide pos = aa_pos * 3
                    if tool in ('blastx', 'diamond'):
                        try:
                            ref_seq = self.gene_sequences.get(gene)
                            # If a reference nucleotide sequence is available we can
                            # sanity-check the conversion (expected ref len ~= slen*3)
                            if ref_seq is not None:
                                if abs(len(ref_seq) - (slen * 3)) <= 3:
                                    self.logger.debug(f"Converting protein coords -> nucleotides for {gene} (ref present)")
                                else:
                                    self.logger.debug(f"Converting protein coords for {gene} but reference length ({len(ref_seq)}) does not match expected {slen*3}; proceeding anyway")
                            else:
                                self.logger.debug(f"No reference sequence for {gene}; assuming diamond/blastx subject coords are amino-acid positions and converting to nucleotides")

                            sstart = sstart * 3
                            send = send * 3
                            slen = slen * 3
                        except Exception:
                            pass

                    # Optional filtering: pident (identity) and query length
                    identity = None
                    qlen = None
                    try:
                        if len(fields) > 2:
                            identity = float(fields[2])
                    except Exception:
                        identity = None
                    try:
                        if len(fields) > 12:
                            qlen = int(fields[12])
                    except Exception:
                        qlen = None

                    # Compute query coverage if possible
                    query_cov = 0.0
                    try:
                        if qlen and qlen > 0:
                            query_cov = (abs(qend - qstart) + 1) / qlen * 100
                    except Exception:
                        query_cov = 0.0

                    # Apply upstream query filters if they are defined in the run parameters
                    if self.query_min_identity is not None and identity is not None:
                        if identity < float(self.query_min_identity):
                            continue
                    if self.query_min_coverage is not None:
                        if query_cov < float(self.query_min_coverage):
                            continue

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

                # Pre-scan CIGAR to compute aligned positions and alignment length for filtering
                try:
                    cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')
                    ref_pos_tmp = ref_start
                    aligned_positions = set()
                    alignment_length = 0
                    for count_str, op in cigar_re.findall(cigar):
                        length = int(count_str)
                        if op in ('M', '=', 'X'):
                            aligned_positions.update(range(ref_pos_tmp, ref_pos_tmp + length))
                            ref_pos_tmp += length
                            alignment_length += length
                        elif op == 'I':
                            alignment_length += length
                        elif op in ('D', 'N'):
                            ref_pos_tmp += length
                            alignment_length += length
                        elif op == 'S':
                            # soft clip - consumes query but not reference
                            alignment_length += length
                except Exception:
                    aligned_positions = set()
                    alignment_length = 0

                # Parse optional tags (NM) to estimate mismatches
                nm = None
                for opt in fields[11:]:
                    if opt.startswith('NM:i:'):
                        try:
                            nm = int(opt.split(':')[-1])
                        except Exception:
                            nm = None
                        break

                # Compute identity and query coverage if possible
                identity = None
                query_cov = 0.0
                try:
                    if alignment_length > 0 and nm is not None:
                        identity = (alignment_length - nm) / alignment_length * 100
                    query_length = len(seq) if seq and seq != '*' else 0
                    if query_length > 0 and aligned_positions:
                        query_cov = len(aligned_positions) / query_length * 100
                except Exception:
                    identity = None
                    query_cov = 0.0

                # Apply upstream query filters if defined
                if self.query_min_identity is not None and identity is not None:
                    if identity < float(self.query_min_identity):
                        continue
                if self.query_min_coverage is not None:
                    if query_cov < float(self.query_min_coverage):
                        continue

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

        # Coverage/report units are in base pairs (bp). For protein tools
        # (blastx/diamond) we converted coordinates to nucleotide positions
        # when parsing so plots are directly comparable across tools.
        marker = 'bp'

        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"GeneFíor Coverage Report\n")
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
        # Generate PNG coverage plot if plotting enabled and individual plotting requested
        try:
            if getattr(self, 'plot_enabled', False) and getattr(self, 'generate_individual_plots', False):
                self.generate_coverage_plot(gene, database, tool, coverage)
        except Exception as e:
            self.logger.debug(f"Failed to generate coverage plot for {gene} ({database}/{tool}): {e}\n{traceback.format_exc()}")

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
            f.write(f"GeneFíor Multi-Tool Comparison\n")
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
        # Generate comparison plot across tools
        try:
            if getattr(self, 'plot_enabled', False):
                self.generate_comparison_plot(gene, database)
        except Exception as e:
            self.logger.debug(f"Failed to generate comparison plot for {gene} ({database}): {e}\n{traceback.format_exc()}")

    def generate_coverage_plot(self, gene: str, database: str, tool: str, coverage: GeneCoverage):
        """Create a PNG coverage plot for a single gene/tool showing depth across the gene and variant markers."""
        try:
            stats = coverage.get_coverage_stats()
            length = coverage.gene_length
            if length <= 0:
                return

            positions = list(range(1, length + 1))
            depths = [coverage.positions[i].depth for i in range(length)]

            fig, ax = self.plt.subplots(figsize=(12, 4))
            ax.plot(positions, depths, color='tab:blue', linewidth=0.8)
            ax.fill_between(positions, depths, color='tab:blue', alpha=0.1)
            ax.set_xlabel('Position')
            ax.set_ylabel('Depth')
            ax.set_title(f"{database} | {tool} | {gene}")
            ax.grid(alpha=0.3)

            # Mark variant positions
            for var in stats.get('variants', []):
                pos = var.get('pos', None)
                if pos is not None and 0 <= pos < length:
                    ax.axvline(x=pos + 1, color='red', alpha=0.6, linestyle='--')

            # Mark uncovered regions (gaps) as shaded areas
            for gap_start, gap_end in stats.get('gaps', []):
                ax.axvspan(gap_start + 1, gap_end + 1, color='grey', alpha=0.15)

            outpath = self.plot_dir / f"{database}_{tool}_{gene.replace('|','_').replace('/','_').replace(':','_')}_coverage.png"
            fig.tight_layout()
            fig.savefig(str(outpath), dpi=150)
            self.plt.close(fig)
            self.logger.info(f"Coverage plot written: {outpath}")
        except Exception as e:
            self.logger.debug(f"Error generating coverage plot for {gene} ({database}/{tool}): {e}")

    def generate_comparison_plot(self, gene: str, database: str):
        """Create a comparison PNG plotting coverage from all tools for a gene.

        Plots per-tool coverage lines on the same axes to visualise agreement.
        """
        try:
            tools_with_data = [tool for tool in self.tools if gene in self.gene_coverages[database].get(tool, {})]
            if not tools_with_data:
                return

            fig, ax = self.plt.subplots(figsize=(12, 4))
            max_depth = 0
            length = None
            for tool in tools_with_data:
                cov = self.gene_coverages[database][tool][gene]
                length = cov.gene_length if length is None else length
                depths = [cov.positions[i].depth for i in range(cov.gene_length)]
                positions = list(range(1, cov.gene_length + 1))
                ax.plot(positions, depths, label=tool, linewidth=1.0)
                max_depth = max(max_depth, max(depths) if depths else 0)

            ax.set_xlabel('Position')
            ax.set_ylabel('Depth')
            ax.set_title(f"{database} | {gene} | tool comparison")
            ax.legend(loc='upper right', fontsize='small')
            ax.grid(alpha=0.3)

            outpath = self.plot_dir / f"{database}_{gene.replace('|','_').replace('/','_').replace(':','_')}_comparison.png"
            fig.tight_layout()
            fig.savefig(str(outpath), dpi=150)
            self.plt.close(fig)
            self.logger.info(f"Comparison plot written: {outpath}")
        except Exception as e:
            self.logger.debug(f"Error generating comparison plot for {gene} ({database}): {e}")

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

    def load_detected_genes(self, database: str):
        """Load detected genes from a per-sample detection matrix in the input directory.

        Looks for <input_dir>/<database>_detection_matrix.tsv and collects genes that have Total_Detections>0
        or any tool column non-zero. Returns a set of gene names. On error or missing file returns None.
        """
        det_path = Path(self.input_dir) / f"{database}_detection_matrix.tsv"
        if not det_path.exists():
            self.logger.warning(f"Detection matrix not found for '{database}': {det_path}")
            return None
        detected = set()
        try:
            import csv
            with open(det_path, 'r', newline='') as fh:
                reader = csv.reader(fh, delimiter='\t')
                header = next(reader, None)
                if not header:
                    return None
                gene_idx = 0
                total_idx = len(header) - 1
                tool_cols = [i for i in range(len(header)) if i != gene_idx and i != total_idx]
                for row in reader:
                    if not row:
                        continue
                    gene = row[gene_idx]
                    # If Total_Detections column exists and >0 treat as detected
                    try:
                        total = int(row[total_idx])
                    except Exception:
                        # fall back to any tool column non-empty/non-zero
                        total = 0
                        for i in tool_cols:
                            if i < len(row) and row[i] and row[i] != '0':
                                total = 1
                                break
                    if total > 0:
                        detected.add(gene)
            return detected
        except Exception as e:
            self.logger.warning(f"Failed to read detection matrix {det_path}: {e}")
            return None

    def load_run_parameters(self):
        """Attempt to load run parameters (JSON) from the input directory or its subdirectories.

        Looks for a file named 'run_parameters.json' (written by GeneFíor / AMRfíor / Recompute)
        and extracts query filtering thresholds to be applied by the visualiser.
        """
        search_root = Path(self.input_dir)
        candidates = list(search_root.rglob('run_parameters.json'))
        if not candidates:
            # nothing found
            return None

        try:
            with open(candidates[0], 'r') as fh:
                params = json.load(fh)
            self.logger.info(f"Loaded run parameters from: {candidates[0]}")
            # Extract relevant thresholds
            self.query_min_coverage = params.get('query_min_coverage', None)
            self.query_min_identity = params.get('query_min_identity', None)
            # Also allow top-level 'q-min-cov' style keys if present
            if self.query_min_coverage is None and 'q-min-cov' in params:
                self.query_min_coverage = params.get('q-min-cov')
            if self.query_min_identity is None and 'q-min-id' in params:
                self.query_min_identity = params.get('q-min-id')

            return params
        except Exception as e:
            self.logger.warning(f"Failed to parse run parameters file: {e}")
            return None

    def run(self):
        # Main execution.
        self.logger.info("=" * 70)
        self.logger.info(f"GeneFíor-Gene-Report {GENEFIOR_VERSION}")
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
            # Load detected genes for this database (if we're reporting only detected genes)
            detected_set = None
            if self.report_only_detected:
                detected_set = self.load_detected_genes(database)
                if detected_set is None:
                    # If detection matrix missing or unreadable, fallback to reporting all genes but warn
                    self.logger.warning(f"Falling back to reporting all genes for database '{database}' because detection matrix could not be read")
                    detected_set = None
            for tool in self.tools:
                if tool not in self.gene_coverages[database]:
                    continue

                for gene, coverage in self.gene_coverages[database][tool].items():
                    # If requested, only report genes that appear in the detection matrix
                    if self.report_only_detected and detected_set is not None and gene not in detected_set:
                        continue
                    self.generate_text_report(gene, database, tool, coverage)
                    report_count += 1

            # Generate comparison reports
            all_genes = set()
            for tool in self.tools:
                if tool in self.gene_coverages[database]:
                    all_genes.update(self.gene_coverages[database][tool].keys())

            for gene in all_genes:
                if self.report_only_detected and detected_set is not None and gene not in detected_set:
                    continue
                self.generate_comparison_report(gene, database)

        # Summary
        self.logger.info("\n" + "=" * 70)
        self.logger.info("GENE REPORT COMPLETE")
        self.logger.info("=" * 70)
        total_reports = len(list(self.report_dir.glob("*.txt")))
        self.logger.info(f"Generated {total_reports} coverage reports")
        self.logger.info(f"Reports saved to: {self.report_dir}")
        # Report plots summary
        try:
            if getattr(self, 'plot_enabled', False):
                total_plots = len(list(self.plot_dir.glob("*.png")))
                self.logger.info(f"Generated {total_plots} PNG plots")
                self.logger.info(f"Plots saved to: {self.plot_dir}")
        except Exception:
            pass
        self.logger.info("=" * 70)

        return True


def main():
    parser = argparse.ArgumentParser(
        description='GeneFíor ' + GENEFIOR_VERSION + ' - GeneFíor-Gene-Stats: Generate detailed coverage visualisations for searched genes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Visualise specific genes (FULL NAMES) from all tools
  GeneFior-gene-stats -i results/ -o vis/ \\
    -g "sul1_2_U12338,tet(W)|ARO:3000194" \\
    --databases resfinder card \\
    --tools diamond bowtie2 bwa

  # Visualise from gene (FULL NAMES) list file with reference
  GeneFior-gene-stats -i results/ -o vis/ \\
    -g genes_of_interest.txt \\
    --databases resfinder \\
    --tools blastn diamond 

        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input directory containing GeneFíor results')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for visualisation reports')
    parser.add_argument('-g', '--genes', required=False,
                        help='Comma-separated gene names (FULL NAMES) or path to file with gene names (one per line)')
    parser.add_argument('--all-genes', action='store_true', default=False,
                        help='Include all genes found in raw outputs (default: only genes listed as detected in detection_matrix.tsv)')
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
    parser.add_argument('--plot-per-tool', action='store_true', default=False,
                        help='Generate individual per-tool coverage PNGs in addition to combined comparison plots (default: off)')


    options = parser.parse_args()

    # Check input directory
    if not os.path.exists(options.input):
        print(f"Error: Input directory '{options.input}' not found", file=sys.stderr)
        sys.exit(1)

        # Tool selection
    if options.tools == ['all']:
        # Expand "all" differently depending on sequence type. For gene-visualiser
        # invocations that operate on gene FASTA inputs we should avoid including read-mappers.
        if getattr(options, 'sequence_type', None) == 'Genes-FASTA':
            genes_type = getattr(options, 'genes_type', None)
            if genes_type == 'aa':
                options.tools = ['blastp', 'diamond']
            elif genes_type == 'dna':
                options.tools = ['blastn', 'diamond']
            else:
                options.tools = ['blastn', 'blastx', 'diamond']
        else:
            options.tools = ['blastn', 'blastx', 'diamond', 'bowtie2', 'bwa', 'minimap2']

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
    visualiser = GeneVisualiser(
        input_dir=options.input,
        output_dir=options.output,
        genes=genes,
        databases=options.databases,
        tools=options.tools,
        ref_fasta=options.ref_fasta,
        query_fasta=options.query_fasta,
        plot_per_tool=options.plot_per_tool
    )
    # Respect user request to include all genes
    if options.all_genes:
        visualiser.report_only_detected = False

    success = visualiser.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
