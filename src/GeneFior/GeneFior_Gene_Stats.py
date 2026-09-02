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
import tempfile
import math

try:
    from .constants import *
    from .evidence import (
        ALLELE_LIKE,
        CANDIDATE_ALLELE_DETECTED,
        EVIDENCE_PRESENT_STATUSES,
        EXACT_ALLELE_DETECTED,
        MIXED_OR_MOSAIC,
        NOT_DETECTED,
        PARTIAL_OR_DIVERGENT,
    )
    from .read_store import ReadMatchStore, _iter_sequences
except (ModuleNotFoundError, ImportError):
    from constants import *
    from evidence import (
        ALLELE_LIKE,
        CANDIDATE_ALLELE_DETECTED,
        EVIDENCE_PRESENT_STATUSES,
        EXACT_ALLELE_DETECTED,
        MIXED_OR_MOSAIC,
        NOT_DETECTED,
        PARTIAL_OR_DIVERGENT,
    )
    from read_store import ReadMatchStore, _iter_sequences


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

    def __init__(self, gene_name: str, gene_length: int, ref_seq: str = None,
                 unit: str = 'bp'):
        self.gene_name = gene_name
        self.gene_length = gene_length
        self.ref_seq = ref_seq
        self.unit = unit
        self.positions = {i: GenePosition(i) for i in range(gene_length)}
        self.read_count = 0
        self.covered_positions = set()

        # Set reference bases if available
        if ref_seq:
            for i, base in enumerate(ref_seq):
                if i < gene_length:
                    self.positions[i].ref_base = base.upper()

    def add_alignment(self, start: int, end: int, aligned_seq: str = None,
                      ref_start_in_seq: int = 0, increment_depth: bool = True,
                      increment_read_count: bool = True):
        # Add an alignment to coverage.
        if increment_read_count:
            self.read_count += 1
        for pos in range(start, end):
            if 0 <= pos < self.gene_length:
                #self.positions[pos].depth += 1
                self.covered_positions.add(pos)
                # Only increment depth if we're not tracking bases separately
                if aligned_seq is None and increment_depth:
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
        coverage_by_depth = {
            threshold: (
                sum(1 for depth in depths if depth >= threshold)
                / self.gene_length * 100
                if self.gene_length > 0 else 0.0
            )
            for threshold in (1, 2, 3, 5, 10)
        }
        if depths and avg_depth > 0:
            variance = sum(
                (depth - avg_depth) ** 2 for depth in depths
            ) / len(depths)
            depth_cv = math.sqrt(variance) / avg_depth
        else:
            depth_cv = 0.0

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
        internal_gaps = [
            (start, end)
            for start, end in gaps
            if start > 0 and end < self.gene_length - 1
        ]

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
            'coverage_2x': coverage_by_depth[2],
            'coverage_3x': coverage_by_depth[3],
            'coverage_5x': coverage_by_depth[5],
            'coverage_10x': coverage_by_depth[10],
            'avg_depth': avg_depth,
            'max_depth': max_depth,
            'depth_cv': depth_cv,
            'read_count': self.read_count,
            'gaps': gaps,
            'internal_gaps': internal_gaps,
            'longest_internal_gap': max(
                (end - start + 1 for start, end in internal_gaps),
                default=0,
            ),
            'variants': variants,
            'positions_with_bases': positions_with_bases,
            'positions_with_ref': positions_with_ref
        }


class GeneVisualiser:
    # Generate coverage visualisations for searched genes.

    def __init__(self, input_dir: str, output_dir: str, genes: List[str],
                 databases: List[str], tools: List[str], ref_fasta: str = None,
                 query_fasta: str = None, plot_per_tool: bool = False,
                 report_read_names: bool = False, report_fasta: bool = False,
                 read_selection: str = "passing",
                 gene_selection: str = "evidence",
                 raw_dir: str = None):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.raw_dir = self._normalise_raw_dir(raw_dir) if raw_dir else None
        self.raw_dir_used = None
        self.genes = set(genes) if genes else set()
        self.gene_selection = self._normalise_gene_selection(gene_selection)
        self.databases = databases
        self.tools = tools
        self.ref_fasta = ref_fasta
        self.query_fasta = query_fasta
        self.report_read_names = bool(report_read_names)
        self.report_fasta = bool(report_fasta)
        self.read_selection = read_selection
        self.read_store = None

        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.report_dir = self.output_dir / "gene_reports"
        self.report_dir.mkdir(exist_ok=True)
        #self.plot_dir = self.output_dir / "coverage_plots"
        self.plot_dir = self.output_dir / "gene_plots"
        self.plot_dir.mkdir(exist_ok=True)
        self.read_output_dir = self.output_dir / "read_outputs"
        if self.report_read_names or self.report_fasta:
            self.read_output_dir.mkdir(exist_ok=True)
            fd, temp_name = tempfile.mkstemp(
                prefix="genefior_gene_stats_reads_",
                suffix=".sqlite",
                dir=str(self.output_dir),
            )
            os.close(fd)
            os.unlink(temp_name)
            self.read_store = ReadMatchStore(Path(temp_name))

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
        self.detected_by_tool = defaultdict(lambda: defaultdict(dict))
        self.evidence_status = defaultdict(lambda: defaultdict(dict))
        self.evidence_warnings = defaultdict(lambda: defaultdict(dict))
        # Default behavior: only report genes selected from the detection/evidence matrix.
        self.report_only_detected = self.gene_selection != 'all'
        self.selected_genes_by_db = {}
        # Query filtering thresholds - will be loaded from run_parameters.json if present
        self.query_min_coverage = None
        self.query_min_identity = None
        self.detection_min_coverage = 80.0
        self.detection_system = 'qualified'
        self.evidence_corroborating_depth = 3
        self.evidence_candidate_depth = 3
        self.evidence_candidate_identity = 98.0
        self.evidence_max_internal_gap_bp = 15
        self.evidence_max_internal_gap_fraction = 0.02
        self.sequence_type = None
        self.genes_type = None
        self.run_parameters = {}
        self.source_raw_outputs = None
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
            if not self.gene_sequences:
                self.logger.warning(
                    "No reference sequences were loaded. Check that the reference "
                    "file is plain FASTA text and starts with '>'."
                )
        except Exception as e:
            self.logger.error(f"Error loading reference FASTA: {e}")

    def load_query_sequences(self):
        # Load query sequences from FASTA or FASTQ.
        query_paths = [
            Path(part.strip())
            for part in str(self.query_fasta or '').split(',')
            if part.strip()
        ]
        if not query_paths or any(not path.is_file() for path in query_paths):
            self.logger.warning("Query FASTA/FASTQ not provided or not found - BLAST base-level data disabled")
            return

        self.logger.info(
            "Loading query sequences from "
            + ', '.join(str(path) for path in query_paths)
        )

        try:
            count = 0
            for index, query_path in enumerate(query_paths):
                mate_suffix = f"/{index + 1}" if len(query_paths) == 2 else ''
                for record_id, sequence in _iter_sequences(query_path):
                    self.query_sequences[record_id] = sequence.upper()
                    if mate_suffix and not record_id.endswith(('/1', '/2')):
                        self.query_sequences[f"{record_id}{mate_suffix}"] = sequence.upper()
                    count += 1

            self.logger.info(f"Loaded {count} query sequences")
        except Exception as e:
            self.logger.error(f"Error loading query FASTA: {e}")
            import traceback
            self.logger.error(traceback.format_exc())

    @staticmethod
    def _protein_subject_to_nt_coords(sstart: int, send: int):
        """Convert inclusive 1-based amino-acid subject coords to nucleotide coords.

        Returns inclusive 1-based nucleotide coordinates while preserving the
        reported subject direction. Amino acid 1 spans bases 1..3, so the left
        edge is (aa - 1) * 3 + 1 and the right edge is aa * 3.
        """
        if sstart <= send:
            return ((sstart - 1) * 3) + 1, send * 3
        return sstart * 3, ((send - 1) * 3) + 1

    @staticmethod
    def _needs_reverse_complement(tool: str, qstart: int, qend: int,
                                  sstart: int, send: int) -> bool:
        """Return True when the extracted query span must be reverse-complemented."""
        if tool in ('blastx', 'diamond'):
            return qstart > qend
        return (qstart > qend) != (sstart > send)

    def _diamond_mode(self, blast_file: Path) -> str:
        if self.sequence_type == 'Genes-FASTA':
            if self.genes_type in ('aa', 'protein'):
                return 'blastp'
            if self.genes_type == 'dna':
                return 'blastx'
        try:
            with open(blast_file, 'r') as handle:
                for _ in range(200):
                    line = handle.readline()
                    if not line:
                        break
                    if line.startswith('#'):
                        lowered = line.lower()
                        if 'blastp' in lowered:
                            return 'blastp'
                        if 'blastx' in lowered:
                            return 'blastx'
                        continue
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 14:
                        continue
                    alignment_length = float(fields[3])
                    query_span = abs(int(fields[7]) - int(fields[6])) + 1
                    if alignment_length > 0:
                        ratio = query_span / alignment_length
                        if 0.75 <= ratio <= 1.25:
                            return 'blastp'
                        if 2.5 <= ratio <= 3.5:
                            return 'blastx'
        except (OSError, TypeError, ValueError):
            pass
        return 'blastx'

    def parse_blast_results(self, blast_file: Path, database: str, tool: str):
        # Parse BLAST format 6 output.
        self.logger.info(f"Parsing BLAST results: {database} - {tool}")
        if not blast_file.is_file():
            raise FileNotFoundError(f"BLAST results file not found: {blast_file}")

        try:
            alignment_count = 0
            alignments_with_seq = 0
            seen_read_gene_positions = defaultdict(set)
            counted_read_genes = set()
            diamond_mode = self._diamond_mode(blast_file) if tool == 'diamond' else None
            translated_search = (
                tool == 'blastx'
                or (tool == 'diamond' and diamond_mode == 'blastx')
            )
            direct_protein_search = (
                tool == 'blastp'
                or (tool == 'diamond' and diamond_mode == 'blastp')
            )

            with open(blast_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 14:
                        continue

                    query_id = fields[0]
                    gene = fields[1]
                    if not self._should_process_gene(database, gene):
                        continue
                    if self.read_store is not None:
                        self.read_store.add(
                            database, tool, gene, query_id,
                            passing=False,
                        )

                    # BLAST coordinates
                    qstart = int(fields[6])  # Query start
                    qend = int(fields[7])  # Query end
                    sstart = int(fields[8])  # Subject start
                    send = int(fields[9])  # Subject end
                    slen = int(fields[13])  # Subject length

                    # If this is a translated/protein-target alignment, the subject
                    # coordinates and length are in amino acids. Convert the subject
                    # interval to nucleotide coordinates so coverage plots align with
                    # nucleotide-based tools (e.g. bwa/bowtie2).
                    if translated_search:
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

                            sstart, send = self._protein_subject_to_nt_coords(sstart, send)
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
                    if self.read_store is not None:
                        self.read_store.add(
                            database, tool, gene, query_id,
                            passing=True,
                        )

                    # Initialise coverage if needed
                    if gene not in self.gene_coverages[database][tool]:
                        ref_seq = self.gene_sequences.get(gene)
                        self.gene_coverages[database][tool][gene] = GeneCoverage(
                            gene, slen, ref_seq,
                            unit='aa' if direct_protein_search else 'bp',
                        )
                        if ref_seq:
                            self.logger.debug(
                                f"Initialised coverage for {gene} with reference sequence (length: {len(ref_seq)})")
                        else:
                            self.logger.debug(
                                f"Initialised coverage for {gene} WITHOUT reference sequence (length: {slen})")

                    # Add alignment (convert to 0-based)
                    start = min(sstart, send) - 1
                    end = max(sstart, send)
                    read_gene = (query_id, gene)
                    observed_positions = seen_read_gene_positions[read_gene]
                    new_positions = [
                        pos for pos in range(start, end)
                        if 0 <= pos < slen and pos not in observed_positions
                    ]
                    if not new_positions:
                        continue
                    observed_positions.update(new_positions)
                    coverage = self.gene_coverages[database][tool][gene]
                    if read_gene not in counted_read_genes:
                        coverage.read_count += 1
                        counted_read_genes.add(read_gene)

                    # Try to get query sequence for this alignment
                    query_seq = None
                    tracked_depth_from_bases = False
                    if query_id in self.query_sequences:
                        full_query = self.query_sequences[query_id]

                        # Extract aligned portion from query (1-based to 0-based)
                        q_start_idx = min(qstart, qend) - 1
                        q_end_idx = max(qstart, qend)

                        if 0 <= q_start_idx < len(full_query) and q_end_idx <= len(full_query):
                            query_seq = full_query[q_start_idx:q_end_idx]

                            # Handle reverse complement if needed
                            orientation_tool = 'blastx' if translated_search else tool
                            if (not direct_protein_search
                                    and self._needs_reverse_complement(
                                        orientation_tool, qstart, qend,
                                        sstart, send)):
                                query_seq = self._reverse_complement(query_seq)

                            alignments_with_seq += 1

                            # Add alignment with sequence data
                            # Map query bases to subject positions
                            for subject_pos in new_positions:
                                query_index = subject_pos - start
                                if 0 <= query_index < len(query_seq):
                                    coverage.add_base_at_position(
                                        subject_pos, query_seq[query_index]
                                    )
                                    tracked_depth_from_bases = True

                    # Multiple HSP rows from one read are one depth observation.
                    run_start = new_positions[0]
                    previous = run_start
                    for subject_pos in new_positions[1:] + [None]:
                        if subject_pos is not None and subject_pos == previous + 1:
                            previous = subject_pos
                            continue
                        coverage.add_alignment(
                            run_start,
                            previous + 1,
                            aligned_seq=None,
                            increment_depth=not tracked_depth_from_bases,
                            increment_read_count=False,
                        )
                        if subject_pos is not None:
                            run_start = subject_pos
                            previous = subject_pos
                    alignment_count += 1

                self.logger.info(f"Processed {alignment_count} alignments from {blast_file}")
                if self.query_sequences and alignment_count:
                    self.logger.info(
                        f"  {alignments_with_seq} alignments with sequence data ({alignments_with_seq / alignment_count * 100:.1f}%)")
                elif self.query_sequences:
                    self.logger.info("  No qualifying alignments with sequence data")
                else:
                    self.logger.warning(
                        f"  Query sequences not loaded - variant calling not available for {database}/{tool}")


        except Exception as e:
            self.logger.error(f"Error parsing BLAST file {blast_file}: {e}")
            raise RuntimeError(f"Failed to parse BLAST file {blast_file}") from e

    def _reverse_complement(self, seq: str) -> str:
        # Return reverse complement of a DNA sequence.
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq))

    def parse_bam_file(self, bam_file: Path, database: str, tool: str):
        # Parse BAM file for coverage with base-level tracking.
        self.logger.info(f"Parsing BAM file: {database} - {tool}")
        if not bam_file.is_file():
            raise FileNotFoundError(f"BAM file not found: {bam_file}")

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
                if flag & 0x100 or flag & 0x800:
                    continue

                read_name = fields[0]
                if flag & 0x40:
                    read_name += '/1'
                elif flag & 0x80:
                    read_name += '/2'
                gene = fields[2]
                if not self._should_process_gene(database, gene):
                    continue

                ref_start = int(fields[3]) - 1  # 0-based
                cigar = fields[5]
                seq = fields[9]
                qual = fields[10] if len(fields) > 10 else None
                if self.read_store is not None:
                    self.read_store.add(
                        database, tool, gene, read_name,
                        passing=False, sequence=(seq if seq and seq != '*' else None),
                    )

                # Pre-scan CIGAR to compute aligned positions and alignment length for filtering
                try:
                    cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')
                    ref_pos_tmp = ref_start
                    aligned_positions = set()
                    alignment_length = 0
                    query_aligned_bases = 0
                    for count_str, op in cigar_re.findall(cigar):
                        length = int(count_str)
                        if op in ('M', '=', 'X'):
                            aligned_positions.update(range(ref_pos_tmp, ref_pos_tmp + length))
                            ref_pos_tmp += length
                            alignment_length += length
                            query_aligned_bases += length
                        elif op == 'I':
                            alignment_length += length
                            query_aligned_bases += length
                        elif op == 'D':
                            ref_pos_tmp += length
                            alignment_length += length
                        elif op == 'N':
                            ref_pos_tmp += length
                except Exception:
                    aligned_positions = set()
                    alignment_length = 0
                    query_aligned_bases = 0

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
                    if query_length > 0:
                        query_cov = query_aligned_bases / query_length * 100
                except Exception:
                    identity = None
                    query_cov = 0.0

                # Apply upstream query filters if defined
                if self.query_min_identity is not None:
                    if identity is None or identity < float(self.query_min_identity):
                        continue
                if self.query_min_coverage is not None:
                    if query_cov < float(self.query_min_coverage):
                        continue
                if self.read_store is not None:
                    self.read_store.add(
                        database, tool, gene, read_name,
                        passing=True, sequence=(seq if seq and seq != '*' else None),
                    )

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
            stderr = proc.stderr.read() if proc.stderr else ''
            if proc.returncode != 0:
                raise RuntimeError(
                    f"samtools view failed with exit code {proc.returncode}: "
                    f"{stderr.strip()}"
                )

            self.logger.info(f"Processed {alignment_count} alignments from {bam_file}")

        except Exception as e:
            self.logger.error(f"Error parsing BAM file {bam_file}: {e}")
            raise RuntimeError(f"Failed to parse BAM file {bam_file}") from e

    def generate_text_report(self, gene: str, database: str, tool: str, coverage: GeneCoverage):
        # Generate detailed text report for a gene.
        stats = coverage.get_coverage_stats()
        qualified_reporting = self.detection_system != 'legacy-relaxed'
        if qualified_reporting:
            evidence_status, evidence_warnings = self.get_evidence_summary(
                gene, database, tool, coverage
            )


        safe_gene = gene.replace('|', '_').replace('/', '_').replace(':','_').replace('-','_')
        if gene != safe_gene:
            self.logger.info(f"Renaming gene for report file:  {gene}  to  {safe_gene}")

        report_file = self.report_dir / f"{database}_{tool}_{safe_gene}_coverage.txt"

        # Coverage/report units are in base pairs (bp). For protein tools
        # (blastx/diamond) we converted coordinates to nucleotide positions
        # when parsing so plots are directly comparable across tools.
        marker = coverage.unit

        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"GeneFíor Coverage Report\n")
            f.write(f"Gene: {gene}\n")
            f.write(f"Database: {database}\n")
            f.write(f"Tool: {tool}\n")
            f.write(f"Detection system: {self.detection_system}\n")
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
            if qualified_reporting:
                f.write(f"Evidence status:       {evidence_status}\n")
                f.write(f"Evidence warnings:     {evidence_warnings or 'None'}\n")
            else:
                f.write(
                    "Detected:              "
                    f"{self._legacy_detected(database, gene, tool)}\n"
                )
            f.write(f"Coverage at 2×:        {stats['coverage_2x']:.2f}%\n")
            f.write(f"Coverage at 3×:        {stats['coverage_3x']:.2f}%\n")
            f.write(f"Coverage at 5×:        {stats['coverage_5x']:.2f}%\n")
            f.write(f"Coverage at 10×:       {stats['coverage_10x']:.2f}%\n")
            f.write(f"Longest internal gap:  {stats['longest_internal_gap']} {marker}\n")
            f.write(f"Depth CV:              {stats['depth_cv']:.3f}\n")

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
            coverage_bins = []

            for i in range(0, coverage.gene_length, bin_size):
                bin_end = min(i + bin_size, coverage.gene_length)
                bin_depths = [coverage.positions[j].depth for j in range(i, bin_end)]
                avg_bin_depth = sum(bin_depths) / len(bin_depths) if bin_depths else 0
                bins.append(avg_bin_depth)
                covered_in_bin = sum(1 for depth in bin_depths if depth > 0)
                if covered_in_bin == 0:
                    coverage_bins.append(" ")
                elif covered_in_bin == len(bin_depths):
                    coverage_bins.append("█")
                else:
                    coverage_bins.append("▒")

            # Draw ASCII histogram
            f.write(
                f"Gene span: 1..{coverage.gene_length:,} {marker} "
                f"(full length {coverage.gene_length:,} {marker}; "
                f"each plotted bin ~= {bin_size:,} {marker})\n"
            )
            f.write(
                f"Position scale ({marker}): start=1 | "
                f"mid={max(1, coverage.gene_length // 2):,} | "
                f"end={coverage.gene_length:,}\n"
            )
            f.write(f"Max depth: {stats['max_depth']}×\n")
            f.write(
                "ASCII histogram: actual mean depth per plotted bin; "
                "Coverage row shows covered/part-covered bins only\n"
            )
            f.write("-" * 80 + "\n")

            depth_levels = self._ascii_depth_levels(stats['max_depth'])
            label_width = max(
                3,
                len(f"{stats['max_depth']:,}"),
                *(len(f"{level:,}") for level in depth_levels),
            )
            for row in depth_levels:
                line = f"{row:{label_width},}× |"
                for val in bins:
                    if val >= row:
                        line += "█"
                    else:
                        line += " "
                f.write(line + "\n")

            f.write(" " * (label_width + 3) + "+" + "-" * len(bins) + "\n")
            f.write("Coverage |" + "".join(coverage_bins) + "|\n")
            f.write("         █ fully covered bin; ▒ partially covered bin\n")

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
            f.write(f"Detection system: {self.detection_system}\n")
            f.write("=" * 80 + "\n\n")

            # Comparison table
            f.write("COVERAGE COMPARISON\n")
            f.write("-" * 120 + "\n")
            qualified_reporting = self.detection_system != 'legacy-relaxed'
            final_column = (
                f"{'Evidence status':<28}"
                if qualified_reporting
                else f"{'Detected':<10}"
            )
            f.write(
                f"{'Tool':<15} {'Coverage %':<12} {'Coverage 2x':<12} "
                f"{'Avg Depth':<12} {'Reads':<10} {'Gaps':<8} "
                f"{final_column}\n"
            )
            f.write("-" * 120 + "\n")

            for tool in tools_with_data:
                coverage = self.gene_coverages[database][tool][gene]
                stats = coverage.get_coverage_stats()
                if qualified_reporting:
                    status, warnings = self.get_evidence_summary(
                        gene, database, tool, coverage
                    )
                else:
                    status = self._legacy_detected(database, gene, tool)
                    warnings = ''

                f.write(
                    f"{tool:<15} {stats['coverage_percent']:<12.2f} "
                    f"{stats['coverage_2x']:<12.2f} "
                    f"{stats['avg_depth']:<12.2f} {stats['read_count']:<10} "
                    f"{len(stats['gaps']):<8} {status:<28}\n"
                )
                if warnings:
                    f.write(f"  Warnings: {warnings}\n")

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

    @staticmethod
    def _format_length(length: int, unit: str) -> str:
        return f"{length:,} {unit}"

    @staticmethod
    def _ascii_depth_levels(max_depth: int, max_rows: int = 10):
        """Return compact, real depth thresholds for the text histogram."""
        if max_depth <= 0:
            return []
        if max_depth <= max_rows:
            return list(range(max_depth, 0, -1))

        step = max(1, math.ceil(max_depth / max_rows))
        levels = []
        level = max_depth
        while level > 0:
            levels.append(level)
            level -= step
        return levels

    @staticmethod
    def _gene_axis_ticks(length: int):
        if length <= 1:
            return [1]
        tick_candidates = [
            1,
            max(1, round(length * 0.25)),
            max(1, round(length * 0.50)),
            max(1, round(length * 0.75)),
            length,
        ]
        return sorted({tick for tick in tick_candidates if 1 <= tick <= length})

    def _label_gene_axis(self, ax, length: int, unit: str, *,
                         xlabel: bool = False):
        """Make the plotted gene span explicit on a Matplotlib axis."""
        if length <= 0:
            return

        ax.set_xlim(0.5, length + 0.5)
        ax.axvline(1, color='black', linewidth=0.8, linestyle='--', alpha=0.55)
        if length > 1:
            ax.axvline(length, color='black', linewidth=0.8,
                       linestyle='--', alpha=0.55)
        ticks = self._gene_axis_ticks(length)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{tick:,}" for tick in ticks])
        if xlabel:
            ax.set_xlabel(
                f"Gene position ({unit}; full length = "
                f"{self._format_length(length, unit)})"
            )

    def generate_coverage_plot(self, gene: str, database: str, tool: str, coverage: GeneCoverage):
        """Create a PNG coverage plot for a single gene/tool showing depth across the gene and variant markers."""
        try:
            stats = coverage.get_coverage_stats()
            length = coverage.gene_length
            if length <= 0:
                return

            positions = list(range(1, length + 1))
            depths = [coverage.positions[i].depth for i in range(length)]
            low_depth_cap = 10
            clipped_depths = [min(depth, low_depth_cap) for depth in depths]
            qualified_reporting = self.detection_system != 'legacy-relaxed'
            if qualified_reporting:
                evidence_status, evidence_warnings = self.get_evidence_summary(
                    gene, database, tool, coverage
                )
                title_status = f"Evidence: {evidence_status}"
            else:
                evidence_warnings = ''
                title_status = (
                    "Detected: "
                    f"{self._legacy_detected(database, gene, tool)}"
                )

            fig, (ax, low_ax) = self.plt.subplots(
                2, 1, figsize=(12, 6), sharex=True, constrained_layout=True,
                gridspec_kw={'height_ratios': [2.4, 1]},
            )
            ax.plot(positions, depths, color='tab:blue', linewidth=0.8)
            ax.fill_between(positions, depths, color='tab:blue', alpha=0.1)
            ax.set_ylabel('Depth')
            ax.set_title(
                f"{database} | {tool} | {gene}\n"
                f"{title_status} | Gene length: "
                f"{self._format_length(length, coverage.unit)}"
            )
            ax.grid(alpha=0.3)

            low_ax.axhline(1, color='black', linewidth=0.7, linestyle=':', alpha=0.7)
            low_ax.plot(positions, clipped_depths, color='tab:blue', linewidth=1.0)
            low_ax.fill_between(positions, clipped_depths, color='tab:blue', alpha=0.16)
            low_ax.set_ylim(-0.5, low_depth_cap + 0.5)
            low_ax.set_ylabel(f'Depth\n(0-{low_depth_cap}×)')
            low_ax.set_title(
                (
                    f"Low-depth view: values above {low_depth_cap}× are clipped"
                    + (
                        f" | Warnings: {evidence_warnings}"
                        if qualified_reporting and evidence_warnings else ""
                    )
                ),
                fontsize='small',
            )
            low_ax.grid(alpha=0.3)
            self._label_gene_axis(ax, length, coverage.unit)
            self._label_gene_axis(low_ax, length, coverage.unit, xlabel=True)

            # Mark variant positions
            for var in stats.get('variants', []):
                pos = var.get('pos', None)
                if pos is not None and 0 <= pos < length:
                    ax.axvline(x=pos + 1, color='red', alpha=0.6, linestyle='--')
                    low_ax.axvline(x=pos + 1, color='red', alpha=0.6, linestyle='--')

            # Mark uncovered regions (gaps) as shaded areas
            for gap_start, gap_end in stats.get('gaps', []):
                ax.axvspan(gap_start + 1, gap_end + 1, color='grey', alpha=0.15)
                low_ax.axvspan(gap_start + 1, gap_end + 1, color='grey', alpha=0.25)

            outpath = self.plot_dir / f"{database}_{tool}_{gene.replace('|','_').replace('/','_').replace(':','_')}_coverage.png"
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

            low_depth_cap = 10
            coverages = [
                (tool, self.gene_coverages[database][tool][gene])
                for tool in tools_with_data
            ]
            plot_length = max(cov.gene_length for _tool, cov in coverages)
            units = {cov.unit for _tool, cov in coverages}
            axis_unit = next(iter(units)) if len(units) == 1 else "mixed units"
            length_labels = {
                self._format_length(cov.gene_length, cov.unit)
                for _tool, cov in coverages
            }
            if len(length_labels) == 1:
                length_summary = f"Gene length: {next(iter(length_labels))}"
            else:
                per_tool_lengths = ", ".join(
                    f"{tool}={self._format_length(cov.gene_length, cov.unit)}"
                    for tool, cov in coverages
                )
                length_summary = f"Gene lengths differ: {per_tool_lengths}"

            fig, (ax, low_ax) = self.plt.subplots(
                2, 1, figsize=(12, 6), sharex=True, constrained_layout=True,
                gridspec_kw={'height_ratios': [2.4, 1]},
            )
            low_ax.axhline(1, color='black', linewidth=0.7, linestyle=':', alpha=0.7)
            for tool, cov in coverages:
                depths = [cov.positions[i].depth for i in range(cov.gene_length)]
                positions = list(range(1, cov.gene_length + 1))
                if self.detection_system == 'legacy-relaxed':
                    status = self._legacy_detected(database, gene, tool)
                else:
                    status, _warnings = self.get_evidence_summary(
                        gene, database, tool, cov
                    )
                line, = ax.plot(
                    positions,
                    depths,
                    label=f"{tool} [{status}]",
                    linewidth=1.0,
                )
                low_ax.plot(
                    positions,
                    [min(depth, low_depth_cap) for depth in depths],
                    color=line.get_color(),
                    linewidth=1.1,
                )

            ax.set_ylabel('Depth')
            ax.set_title(f"{database} | {gene} | tool comparison\n{length_summary}")
            ax.legend(loc='upper right', fontsize='small')
            ax.grid(alpha=0.3)

            low_ax.set_ylim(-0.5, low_depth_cap + 0.5)
            low_ax.set_ylabel(f'Depth\n(0-{low_depth_cap}×)')
            low_ax.set_title(
                f"Low-depth view: values above {low_depth_cap}× are clipped",
                fontsize='small',
            )
            low_ax.grid(alpha=0.3)
            self._label_gene_axis(ax, plot_length, axis_unit)
            self._label_gene_axis(low_ax, plot_length, axis_unit, xlabel=True)

            outpath = self.plot_dir / f"{database}_{gene.replace('|','_').replace('/','_').replace(':','_')}_comparison.png"
            fig.savefig(str(outpath), dpi=150)
            self.plt.close(fig)
            self.logger.info(f"Comparison plot written: {outpath}")
        except Exception as e:
            self.logger.debug(f"Error generating comparison plot for {gene} ({database}): {e}")

    @staticmethod
    def _normalise_raw_dir(path: str) -> Path:
        raw_path = Path(path).expanduser()
        if raw_path.name == "raw_outputs":
            return raw_path
        nested = raw_path / "raw_outputs"
        return nested if nested.is_dir() else raw_path

    def _candidate_raw_dirs(self):
        candidates = []
        if self.raw_dir is not None:
            candidates.append(self.raw_dir)

        candidates.append(self.input_dir / "raw_outputs")

        if self.source_raw_outputs:
            candidates.append(self._normalise_raw_dir(self.source_raw_outputs))

        source_input = self.run_parameters.get("source_input")
        if source_input:
            candidates.append(self._normalise_raw_dir(source_input))

        unique = []
        seen = set()
        for candidate in candidates:
            try:
                key = str(candidate.resolve())
            except Exception:
                key = str(candidate)
            if key not in seen:
                unique.append(candidate)
                seen.add(key)
        return unique

    def _discover_files_in_raw_dir(self, raw_dir: Path):
        found_files = []
        for database in self.databases:
            for tool in self.tools:
                if tool in ['blastn', 'blastx', 'blastp', 'diamond']:
                    pattern = f"{database}_{tool}_results.tsv"
                    files = list(raw_dir.glob(pattern))
                    if files:
                        found_files.append((database, tool, 'blast', files[0]))
                        self.logger.info(f"Found BLAST file: {files[0]}")
                elif tool in ['bowtie2', 'bwa', 'minimap2']:
                    pattern = f"{database}_{tool}_results_sorted.bam"
                    files = list(raw_dir.glob(pattern))
                    if files:
                        found_files.append((database, tool, 'bam', files[0]))
                        self.logger.info(f"Found BAM file: {files[0]}")
        return found_files

    def discover_files(self):
        # Discover alignment files in input directory.
        self.logger.info("Discovering alignment files...")

        searched = []
        for raw_dir in self._candidate_raw_dirs():
            searched.append(str(raw_dir))
            if not raw_dir.is_dir():
                self.logger.warning(f"Raw outputs directory not found: {raw_dir}")
                continue
            found_files = self._discover_files_in_raw_dir(raw_dir)
            if found_files:
                self.raw_dir_used = raw_dir
                self.logger.info(f"Using raw outputs from: {raw_dir}")
                self.logger.info(f"Found {len(found_files)} alignment files")
                self.found_files = found_files
                return True
            self.logger.warning(f"No alignment files found in: {raw_dir}")

        self.logger.error(
            "No alignment files found. Gene-Stats needs raw BLAST/DIAMOND "
            "TSV or sorted BAM outputs. Searched: " + "; ".join(searched)
        )
        self.logger.error(
            "If this input is a recompute output, pass --raw-dir pointing to "
            "the original run's raw_outputs directory."
        )
        return False

    def load_detected_genes(self, database: str):
        """Load genes selected by the configured evidence tier.

        New runs use <database>_evidence_matrix.tsv so ambiguous family or
        mosaic calls that pass the user thresholds remain visible in
        Gene-Stats. Candidate/exact modes require this qualified matrix.
        Older evidence-mode runs fall back to the strict binary detection matrix.
        """
        evidence_path = Path(self.input_dir) / f"{database}_evidence_matrix.tsv"
        if self.detection_system != 'legacy-relaxed' and evidence_path.exists():
            try:
                import csv
                selected_genes = set()
                with open(evidence_path, 'r', newline='') as handle:
                    reader = csv.DictReader(handle, delimiter='\t')
                    fixed_columns = {
                        'Gene', 'Detection_System', 'Best_Evidence_Status',
                        'Evidence_Tools', 'Evidence_Detections',
                        'Candidate_Allele_Detections',
                        'Exact_Allele_Detections',
                        'Exact_Protein_Detections',
                        'Profile_Detections', 'Strict_Detections',
                        'Always_Flagged', 'Evidence_Warnings',
                    }
                    tool_columns = [
                        column for column in (reader.fieldnames or [])
                        if column not in fixed_columns
                    ]
                    for row in reader:
                        gene = row.get('Gene', '')
                        best = row.get('Best_Evidence_Status', NOT_DETECTED)
                        statuses = [
                            row.get(tool, NOT_DETECTED) or NOT_DETECTED
                            for tool in tool_columns
                        ]
                        if gene and self._row_matches_gene_selection(
                                row, best, statuses):
                            selected_genes.add(gene)
                        warnings = row.get('Evidence_Warnings', '')
                        for tool in tool_columns:
                            normalised = self._normalise_tool_name(tool)
                            status = row.get(tool, NOT_DETECTED) or NOT_DETECTED
                            self.evidence_status[database][gene][normalised] = status
                            self.evidence_warnings[database][gene][normalised] = warnings
                            self.detected_by_tool[database][gene][normalised] = (
                                status in EVIDENCE_PRESENT_STATUSES
                            )
                self.logger.info(
                    f"Loaded {len(selected_genes)} genes for selection "
                    f"'{self.gene_selection}' from {evidence_path}"
                )
                return selected_genes
            except Exception as e:
                self.logger.warning(
                    f"Failed to read evidence matrix {evidence_path}: {e}"
                )

        if self.gene_selection in (
                'candidate', 'exact', 'candidate-or-exact'):
            self.logger.warning(
                f"Cannot select '{self.gene_selection}' genes for database "
                f"'{database}' because qualified evidence matrix was not found: "
                f"{evidence_path}"
            )
            return set()

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
                    for index in tool_cols:
                        if index >= len(row):
                            continue
                        tool = self._normalise_tool_name(header[index])
                        self.detected_by_tool[database][gene][tool] = (
                            row[index].strip() == '1'
                        )
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

    def _legacy_detected(self, database: str, gene: str, tool: str) -> str:
        """Return the legacy binary detection value for one gene/tool."""
        normalised_tool = self._normalise_tool_name(tool)
        return (
            '1'
            if self.detected_by_tool[database][gene].get(
                normalised_tool, False
            )
            else '0'
        )

    @staticmethod
    def _normalise_gene_selection(value: str) -> str:
        name = str(value or 'evidence').strip().lower().replace('_', '-')
        aliases = {
            'detected': 'evidence',
            'candidate-exact': 'candidate-or-exact',
            'candidate-and-exact': 'candidate-or-exact',
            'allele': 'candidate-or-exact',
            'alleles': 'candidate-or-exact',
        }
        name = aliases.get(name, name)
        allowed = {
            'evidence', 'candidate', 'exact',
            'candidate-or-exact', 'all',
        }
        if name not in allowed:
            raise ValueError(
                f"Unsupported gene selection '{value}'. Choose from: "
                f"{', '.join(sorted(allowed))}"
            )
        return name

    @staticmethod
    def _int_field(row: dict, key: str) -> int:
        try:
            return int(row.get(key, 0) or 0)
        except (TypeError, ValueError):
            return 0

    def _row_matches_gene_selection(
            self, row: dict, best_status: str, statuses: List[str]) -> bool:
        if self.gene_selection == 'all':
            return True

        evidence_present = (
            self._int_field(row, 'Evidence_Detections') > 0
            or best_status in EVIDENCE_PRESENT_STATUSES
            or any(status in EVIDENCE_PRESENT_STATUSES for status in statuses)
        )
        always_flagged = self._int_field(row, 'Always_Flagged') > 0
        candidate_present = (
            self._int_field(row, 'Candidate_Allele_Detections') > 0
            or best_status == CANDIDATE_ALLELE_DETECTED
            or any(status == CANDIDATE_ALLELE_DETECTED for status in statuses)
        )
        exact_present = (
            self._int_field(row, 'Exact_Allele_Detections') > 0
            or best_status == EXACT_ALLELE_DETECTED
            or any(status == EXACT_ALLELE_DETECTED for status in statuses)
        )

        if self.gene_selection == 'evidence':
            return evidence_present or always_flagged
        if self.gene_selection == 'candidate':
            return candidate_present
        if self.gene_selection == 'exact':
            return exact_present
        if self.gene_selection == 'candidate-or-exact':
            return candidate_present or exact_present
        return False

    def _should_process_gene(self, database: str, gene: str) -> bool:
        if self.genes and gene not in self.genes:
            return False
        if self.gene_selection == 'all':
            return True
        selected = self.selected_genes_by_db.get(database)
        if selected is not None and gene not in selected:
            return False
        return True

    @staticmethod
    def _normalise_tool_name(tool: str) -> str:
        name = str(tool).lower()
        if 'diamond' in name:
            return 'diamond'
        if 'blastx' in name:
            return 'blastx'
        if 'blastn' in name:
            return 'blastn'
        if 'blastp' in name:
            return 'blastp'
        if 'bowtie2' in name:
            return 'bowtie2'
        if name == 'bwa' or 'bwa-' in name:
            return 'bwa'
        if 'minimap2' in name:
            return 'minimap2'
        if 'hmmer' in name:
            return name
        return name

    def load_run_parameters(self):
        """Attempt to load run parameters (JSON) from the input directory or its subdirectories.

        Looks for a file named 'run_parameters.json' (written by GeneFíor / AMRfíor / Recompute)
        and extracts query filtering thresholds to be applied by the visualiser.
        """
        search_root = Path(self.input_dir)
        direct_manifest = search_root / 'run_parameters.json'
        candidates = (
            [direct_manifest]
            if direct_manifest.is_file()
            else sorted(search_root.rglob('run_parameters.json'))
        )
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
            self.detection_system = params.get(
                'detection_system', self.detection_system
            )
            self.sequence_type = params.get('sequence_type')
            self.genes_type = params.get('genes_type')
            self.run_parameters = params
            self.source_raw_outputs = params.get('source_raw_outputs')
            self.detection_min_coverage = float(
                params.get('detection_min_coverage', self.detection_min_coverage)
            )
            self.evidence_corroborating_depth = int(
                params.get(
                    'evidence_corroborating_depth',
                    self.evidence_corroborating_depth,
                )
            )
            self.evidence_candidate_depth = int(
                params.get(
                    'evidence_candidate_depth',
                    self.evidence_candidate_depth,
                )
            )
            self.evidence_candidate_identity = float(
                params.get(
                    'evidence_candidate_identity',
                    self.evidence_candidate_identity,
                )
            )
            self.evidence_max_internal_gap_bp = int(
                params.get(
                    'evidence_max_internal_gap_bp',
                    self.evidence_max_internal_gap_bp,
                )
            )
            self.evidence_max_internal_gap_fraction = float(
                params.get(
                    'evidence_max_internal_gap_fraction',
                    self.evidence_max_internal_gap_fraction,
                )
            )
            # Also allow top-level 'q-min-cov' style keys if present
            if self.query_min_coverage is None and 'q-min-cov' in params:
                self.query_min_coverage = params.get('q-min-cov')
            if self.query_min_identity is None and 'q-min-id' in params:
                self.query_min_identity = params.get('q-min-id')

            return params
        except Exception as e:
            self.logger.warning(f"Failed to parse run parameters file: {e}")
            return None

    def get_evidence_summary(self, gene: str, database: str, tool: str,
                             coverage: GeneCoverage):
        normalised_tool = self._normalise_tool_name(tool)
        status = self.evidence_status[database][gene].get(normalised_tool)
        warnings = self.evidence_warnings[database][gene].get(
            normalised_tool,
            '',
        )
        if status:
            return status, warnings

        stats = coverage.get_coverage_stats()
        corroborating_key = f"coverage_{self.evidence_corroborating_depth}x"
        corroborating = stats.get(corroborating_key)
        if corroborating is None:
            corroborating = (
                sum(
                    1 for position in coverage.positions.values()
                    if position.depth >= self.evidence_corroborating_depth
                )
                / coverage.gene_length * 100
                if coverage.gene_length > 0 else 0.0
            )
        max_gap = max(
            self.evidence_max_internal_gap_bp,
            int(round(
                coverage.gene_length
                * self.evidence_max_internal_gap_fraction
            )),
        )
        local_warnings = []
        if stats['coverage_percent'] < self.detection_min_coverage:
            status = PARTIAL_OR_DIVERGENT
            local_warnings.append('PARTIAL_COVERAGE')
        elif (corroborating < self.detection_min_coverage
              or stats['longest_internal_gap'] > max_gap):
            status = MIXED_OR_MOSAIC
            local_warnings.extend([
                'LOW_CORROBORATED_COVERAGE',
                'DISCONTINUOUS_COVERAGE',
            ])
        else:
            status = ALLELE_LIKE
            local_warnings.append('COMPETITIVE_ALLELE_SUPPORT_NOT_AVAILABLE')
            if normalised_tool in ('blastx', 'blastp', 'diamond'):
                local_warnings.append('PROTEIN_ONLY_ALLELE_AMBIGUITY')
        return status, ';'.join(sorted(set(local_warnings)))

    def export_read_outputs(self):
        if self.read_store is None:
            return
        self.read_store.flush()
        if self.report_fasta and self.query_fasta:
            matched_records = 0
            updated_rows = 0
            query_paths = [
                Path(part.strip())
                for part in str(self.query_fasta).split(',')
                if part.strip()
            ]
            for index, query_path in enumerate(query_paths):
                mate_suffix = f"/{index + 1}" if len(query_paths) == 2 else None
                matched, updated = self.read_store.populate_sequences(
                    query_path, suffix=mate_suffix
                )
                matched_records += matched
                updated_rows += updated
            self.logger.info(
                f"Matched {matched_records} query records to {updated_rows} "
                "stored read-gene matches."
            )
        if self.report_read_names:
            exported = self.read_store.export_names(
                self.read_output_dir,
                self.read_selection,
            )
            self.logger.info(f"Exported {exported} selected read names.")
        if self.report_fasta:
            exported, missing = self.read_store.export_fasta(
                self.read_output_dir,
                self.read_selection,
            )
            self.logger.info(
                f"Exported {exported} selected FASTA records; "
                f"{missing} records had no sequence available."
            )

    def run(self):
        # Main execution.
        self.logger.info("=" * 70)
        self.logger.info(f"GeneFíor-Gene-Report {GENEFIOR_VERSION}")
        self.logger.info("=" * 70)
        self.logger.info(f"Input directory: {self.input_dir}")
        self.logger.info(f"Output directory: {self.output_dir}")
        genes_label = (
            f"{len(self.genes)} explicit gene(s)"
            if self.genes
            else f"selection='{self.gene_selection}'"
        )
        self.logger.info(f"Genes to visualise: {genes_label}")
        self.logger.info(f"Databases: {', '.join(self.databases)}")
        self.logger.info(f"Tools: {', '.join(self.tools)}")
        self.logger.info("=" * 70)

        # Discover files
        if not self.discover_files():
            if self.read_store is not None:
                self.read_store.close(delete=True)
                self.read_store = None
            return False

        # Load reference sequences
        if self.ref_fasta:
            self.load_reference_sequences()

        # Load query sequences only after raw outputs are confirmed present.
        if self.query_fasta:
            self.load_query_sequences()
        else:
            self.logger.warning("No query FASTA provided - BLAST variant calling will be limited")

        self.selected_genes_by_db = {
            database: self.load_detected_genes(database)
            for database in self.databases
        }

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
            # Use the preloaded evidence metadata as the default report filter
            # unless all genes were explicitly requested.
            selection_set = self.selected_genes_by_db.get(database)
            detected_set = None
            if self.report_only_detected:
                detected_set = selection_set
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

        try:
            self.export_read_outputs()
        finally:
            if self.read_store is not None:
                self.read_store.close(delete=True)
                self.read_store = None

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

  # Visualise every candidate or exact nucleotide allele call
  GeneFior-gene-stats -i results/ -o vis_candidate_exact/ \\
    --gene-selection candidate-or-exact \\
    --databases resfinder card \\
    --tools blastn bwa minimap2

        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input directory containing GeneFíor results')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for visualisation reports')
    parser.add_argument('-g', '--genes', required=False,
                        help='Comma-separated gene names (FULL NAMES) or path to file with gene names (one per line)')
    parser.add_argument('--all-genes', action='store_true', default=False,
                        help='Include all genes found in raw outputs (default: genes with meaningful evidence '
                             'in evidence_matrix.tsv, falling back to strict detections for older runs)')
    parser.add_argument('--gene-selection',
                        choices=[
                            'evidence', 'candidate', 'exact',
                            'candidate-or-exact', 'all',
                        ],
                        default=None,
                        help='Select genes automatically from qualified evidence outputs: '
                             '"evidence" reports threshold-passing evidence genes (default), '
                             '"candidate" reports candidate allele genes, "exact" reports exact allele genes, '
                             '"candidate-or-exact" reports either tier, and "all" reports all raw-output genes. '
                             '--all-genes is an alias for --gene-selection all.')
    parser.add_argument('--databases', nargs='+', required=True,
                        help='Database name(s) to interrogate, or "all" to discover them from raw outputs')
    parser.add_argument('--tools', nargs='+', required=True,
                        choices=['blastn', 'blastx', 'blastp', 'diamond', 'bowtie2', 'bwa', 'minimap2', 'all'],
                        help='Tool(s) to interrogate')
    parser.add_argument('--ref-fasta',
                        help='Plain-text reference FASTA file for variant calling (optional)')
    parser.add_argument('--query-fasta',
                        help='Original query FASTA/FASTQ file. Required for FASTA read export from '
                             'BLAST/DIAMOND results and optional for base-level BLAST analysis.')
    parser.add_argument('--raw-dir',
                        help='Directory containing raw GeneFior alignment outputs. '
                             'May point directly to raw_outputs/ or to a result directory containing raw_outputs/. '
                             'Useful when --input is a recompute output that did not copy large raw files.')
    parser.add_argument('--plot-per-tool', action='store_true', default=False,
                        help='Generate individual per-tool coverage PNGs in addition to combined comparison plots (default: off)')
    parser.add_argument('--report-read-names', action='store_true', default=False,
                        help='Write deduplicated read-name files for the selected genes, databases, and tools.')
    parser.add_argument('--report-fasta', action='store_true', default=False,
                        help='Write deduplicated read FASTA files for the selected genes, databases, and tools. '
                             'BLAST/DIAMOND exports require --query-fasta; BAM-based tools use sequences in the BAM.')
    parser.add_argument('--read-selection', choices=['all', 'passing'], default='passing',
                        help='Which read-gene matches to export: all raw mappings or only mappings that pass '
                             'the original query identity/coverage thresholds (default: passing).')


    options = parser.parse_args()

    # Check input directory
    if not os.path.exists(options.input):
        print(f"Error: Input directory '{options.input}' not found", file=sys.stderr)
        sys.exit(1)

        # Tool selection
    if options.tools == ['all']:
        options.tools = [
            'blastn', 'blastx', 'blastp', 'diamond',
            'bowtie2', 'bwa', 'minimap2',
        ]

    if options.databases == ['all']:
        raw_dir = (
            GeneVisualiser._normalise_raw_dir(options.raw_dir)
            if options.raw_dir
            else Path(options.input) / 'raw_outputs'
        )
        suffixes = (
            '_blastn_results.tsv', '_blastx_results.tsv',
            '_blastp_results.tsv', '_diamond_results.tsv',
            '_bowtie2_results_sorted.bam', '_bwa_results_sorted.bam',
            '_minimap2_results_sorted.bam',
        )
        discovered_databases = set()
        for path in raw_dir.iterdir() if raw_dir.is_dir() else []:
            for suffix in suffixes:
                if path.name.endswith(suffix):
                    discovered_databases.add(path.name[:-len(suffix)])
                    break
        options.databases = sorted(discovered_databases)
        if not options.databases:
            print("Error: no databases could be discovered from raw outputs", file=sys.stderr)
            sys.exit(1)

    if options.report_fasta and any(
            tool in options.tools for tool in ['blastn', 'blastx', 'blastp', 'diamond']):
        if not options.query_fasta:
            print(
                "Error: --query-fasta is required for --report-fasta when "
                "BLAST/DIAMOND tools are selected",
                file=sys.stderr,
            )
            sys.exit(1)
        query_paths = [
            Path(part.strip())
            for part in str(options.query_fasta).split(',')
            if part.strip()
        ]
        if not query_paths or any(not path.is_file() for path in query_paths):
            print(
                f"Error: query FASTA/FASTQ not found: {options.query_fasta}",
                file=sys.stderr,
            )
            sys.exit(1)
    explicit_gene_selection = options.gene_selection is not None
    if options.all_genes:
        if options.gene_selection not in (None, 'all'):
            print(
                "Error: --all-genes cannot be combined with a different "
                "--gene-selection value",
                file=sys.stderr,
            )
            sys.exit(1)
        options.gene_selection = 'all'
    if options.gene_selection is None:
        options.gene_selection = 'evidence'

    if ((options.report_read_names or options.report_fasta)
            and not options.genes
            and not options.all_genes
            and not explicit_gene_selection):
        print(
            "Error: read export requires --genes, --gene-selection, "
            "or explicitly use --all-genes",
            file=sys.stderr,
        )
        sys.exit(1)

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
        plot_per_tool=options.plot_per_tool,
        report_read_names=options.report_read_names,
        report_fasta=options.report_fasta,
        read_selection=options.read_selection,
        gene_selection=options.gene_selection,
        raw_dir=options.raw_dir,
    )

    success = visualiser.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
