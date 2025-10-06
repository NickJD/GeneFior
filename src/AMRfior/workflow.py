import subprocess, sys
import csv
from collections import defaultdict
from pathlib import Path
import logging
from datetime import datetime
from typing import Any, Dict, List, Set, Tuple

try:
    from .gene_stats import GeneStats
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from gene_stats import GeneStats
    from constants import *

class AMRWorkflow:
    """Orchestrates multiple alignment tools for AMR gene detection."""

    def __init__(self, input_fasta: str, output_dir: str,
                 resfinder_dbs: Dict[str, str], card_dbs: Dict[str, str],
                 threads: int = 4, #max_target_seqs: int = 100,
                 tool_sensitivity_params: Dict[str, Dict[str, Any]] = None,
                 #evalue: float = 1e-10,
                 min_coverage: float = 80.0, min_identity: float = 80.0,
                 run_dna: bool = True, run_protein: bool = True):
        self.input_fasta = Path(input_fasta)
        self.output_dir = Path(output_dir)
        self.resfinder_dbs = resfinder_dbs
        self.card_dbs = card_dbs
        self.threads = threads
      #  self.max_target_seqs = max_target_seqs
        self.tool_sensitivity_params = tool_sensitivity_params
       # self.evalue = evalue
        self.min_coverage = min_coverage
        self.min_identity = min_identity
        self.run_dna = run_dna
        self.run_protein = run_protein

        # Create output directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.raw_dir = self.output_dir / "raw_outputs"
        self.raw_dir.mkdir(exist_ok=True)
        self.stats_dir = self.output_dir / "tool_stats"
        self.stats_dir.mkdir(exist_ok=True)

        # Setup logging
        log_file = self.output_dir / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

        # Store detection results: {database: {gene: {tool: bool}}}
        self.detections = {
            'resfinder': defaultdict(lambda: defaultdict(bool)),
            'card': defaultdict(lambda: defaultdict(bool))
        }

        # Store detailed statistics: {database: {tool: {gene: GeneStats}}}
        self.gene_stats = {
            'resfinder': defaultdict(lambda: defaultdict(GeneStats)),
            'card': defaultdict(lambda: defaultdict(GeneStats))
        }

    def run_command(self, cmd: List[str], tool_name: str) -> bool:
        """Run a tool and log the results."""
        self.logger.info(f"Running {tool_name}...")
        self.logger.info(f"Parameters for {tool_name}: {' '.join(cmd)}")
        self.logger.debug(f"Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            self.logger.info(f"{tool_name} completed successfully")
            if result.stdout:
                self.logger.debug(f"{tool_name} stdout: {result.stdout}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"{tool_name} failed with return code {e.returncode}")
            self.logger.error(f"Error message: {e.stderr}")
            return False
        except FileNotFoundError:
            self.logger.error(f"{tool_name} executable not found. Is it in your PATH?")
            return False

    def parse_blast_results(self, output_file: Path, database: str, tool_name: str) -> Set[str]:
        """Parse BLAST/DIAMOND tabular output and extract genes meeting thresholds.
        Detection logic:
        - Only sequences with identity >= min_identity are considered
        - Track which positions on the subject/gene are covered by alignments
        - Gene is detected if combined coverage of gene >= min_coverage
        """
        detected_genes = set()
        gene_lengths = {}  # Store gene lengths

        if not output_file.exists():
            return detected_genes

        try:
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) < 13:
                        continue

                    gene = fields[1]  # sseqid
                    identity = float(fields[2])  # pident
                    sstart = int(fields[8])  # subject start
                    send = int(fields[9])  # subject end
                    slen = int(fields[12])  # subject length (added to output format)

                    # Store gene length
                    if gene in gene_lengths:
                        gene_lengths[gene] = max(gene_lengths[gene], slen)
                    else:
                        gene_lengths[gene] = slen

                    # Only process sequences meeting identity threshold
                    if identity >= self.min_identity:
                        # Initialise stats if first hit for this gene
                        if gene not in self.gene_stats[database][tool_name]:
                            self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)

                        # Add hit to statistics
                        self.gene_stats[database][tool_name][gene].add_hit(
                            sstart, send, identity, gene_lengths[gene]
                        )

            # finalise statistics and determine detection based on gene coverage
            for gene in self.gene_stats[database][tool_name]:
                stats = self.gene_stats[database][tool_name][gene]
                stats.finalise()

                # Gene is detected if gene coverage meets threshold
                if stats.gene_coverage >= self.min_coverage:
                    detected_genes.add(gene)
                    self.detections[database][gene][tool_name] = True

        except Exception as e:
            self.logger.error(f"Error parsing {output_file}: {e}")

        return detected_genes

    def parse_sam_results(self, sam_file: Path, database: str, tool_name: str) -> Set[str]:
        """Parse SAM file from Bowtie2 and extract genes meeting thresholds."""
        detected_genes = set()
        gene_lengths = {}  # Store gene lengths from @SQ headers

        if not sam_file.exists():
            return detected_genes

        try:
            with open(sam_file, 'r') as f:
                for line in f:
                    if line.startswith('@'):
                        # Extract reference lengths from @SQ header lines
                        if line.startswith('@SQ'):
                            parts = line.strip().split('\t')
                            ref_name = None
                            ref_len = None
                            for part in parts:
                                if part.startswith('SN:'):
                                    ref_name = part[3:]
                                elif part.startswith('LN:'):
                                    ref_len = int(part[3:])
                            if ref_name and ref_len:
                                gene_lengths[ref_name] = ref_len
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 11:
                        continue

                    #query = fields[0]
                    gene = fields[2]
                    if gene == '*':  # Unmapped
                        continue

                    pos = int(fields[3])  # 1-based leftmost mapping position
                    cigar = fields[5]
                   # seq = fields[9]

                    # Parse CIGAR to get alignment positions on reference
                    import re
                    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

                    ref_pos = pos
                    ref_end = pos
                    aligned_bases = 0  # For identity calculation

                    for length_str, op in cigar_ops:
                        length = int(length_str)
                        # Operations that consume reference: M, D, N, =, X
                        if op in 'MDN=X':
                            ref_end = ref_pos + length
                            ref_pos = ref_end
                        # M, =, X consume both query and reference (alignment)
                        if op in 'M=X':
                            aligned_bases += length

                    sstart = pos
                    send = ref_end - 1  # -1 because ref_end is exclusive

                    # Get gene length
                    gene_len = gene_lengths.get(gene, send)  # Use end pos if length unknown

                    # Calculate identity from NM tag (edit distance)
                    nm = 0
                    for tag in fields[11:]:
                        if tag.startswith('NM:i:'):
                            nm = int(tag.split(':')[2])
                            break

                    # Identity = (aligned_bases - mismatches) / aligned_bases
                    identity = ((aligned_bases - nm) / aligned_bases) * 100 if aligned_bases > 0 else 0

                    # Only process sequences meeting identity threshold
                    if identity >= self.min_identity:
                        # Initialise stats if first hit for this gene
                        if gene not in self.gene_stats[database][tool_name]:
                            self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)

                        # Add hit to statistics
                        self.gene_stats[database][tool_name][gene].add_hit(
                            sstart, send, identity, gene_len
                        )

            # finalise statistics and determine detection based on gene coverage
            for gene in self.gene_stats[database][tool_name]:
                stats = self.gene_stats[database][tool_name][gene]
                stats.finalise()

                # Gene is detected if gene coverage meets threshold
                if stats.gene_coverage >= self.min_coverage:
                    detected_genes.add(gene)
                    self.detections[database][gene][tool_name] = True

        except Exception as e:
            self.logger.error(f"Error parsing {sam_file}: {e}")

        return detected_genes



    def parse_hmmer_results(self, tbl_file: Path, database: str, tool_name: str) -> Set[str]:
        """Parse HMMER table output and extract genes meeting thresholds."""
        detected_genes = set()
        if not tbl_file.exists():
            return detected_genes

        try:
            with open(tbl_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split()
                    if len(fields) < 6:
                        continue

                    gene = fields[0]  # target name
                    evalue = float(fields[4])
                    score = float(fields[5]) if len(fields) > 5 else 0.0

                    # Initialise stats if first hit for this gene
                    if gene not in self.gene_stats[database][tool_name]:
                        self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)

                    # For HMMER, use score as proxy for coverage/identity
                    # This is not perfect but HMMER doesn't give direct coverage
                    self.gene_stats[database][tool_name][gene].add_hit(score, score)

                    # HMMER doesn't directly give coverage/identity like BLAST
                    # Use E-value as primary filter
                    if evalue <= self.evalue:
                        detected_genes.add(gene)
                        self.detections[database][gene][tool_name] = True

            # finalise statistics
            for gene in self.gene_stats[database][tool_name]:
                self.gene_stats[database][tool_name][gene].finalise()

        except Exception as e:
            self.logger.error(f"Error parsing {tbl_file}: {e}")

        return detected_genes

    def write_tool_stats(self, database: str, tool_name: str):
        """Write detailed statistics for a specific tool to TSV.

        Output columns:
        - Gene: AMR gene name
        - Gene_Length: Length of the gene in the database (bp)
        - Num_Sequences_Mapped: Number of sequences that mapped to this gene with identity >= min_identity
        - Gene_Coverage: Percentage of the gene covered by all qualifying alignments combined (%)
        - Avg_Identity: Average identity across all qualifying sequences (%)
        - Detected: 1 if gene_coverage >= min_coverage threshold, 0 otherwise
        """
        stats_file = self.stats_dir / f"{database}_{tool_name}_stats.tsv"

        gene_stats = self.gene_stats[database][tool_name]
        if not gene_stats:
            self.logger.warning(f"No statistics to write for {database} - {tool_name}")
            return

        with open(stats_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            # Header
            header = ['Gene', 'Gene_Length', 'Num_Sequences_Mapped', 'Gene_Coverage', 'Avg_Identity', 'Detected']
            writer.writerow(header)

            # Sort genes alphabetically
            genes = sorted(gene_stats.keys())

            for gene in genes:
                stats = gene_stats[gene]
                detected = self.detections[database][gene][tool_name]

                row = [
                    gene,
                    stats.gene_length,
                    stats.num_sequences,
                    f"{stats.gene_coverage:.2f}",
                    f"{stats.avg_identity:.2f}",
                    '1' if detected else '0'
                ]
                writer.writerow(row)

        self.logger.info(f"  Stats file: {stats_file}")

    def run_blast(self, db_path: str, database: str, mode: str) -> Tuple[bool, Set[str]]:
        """Run BLAST in DNA or protein mode."""
        if not db_path:
            return False, set()

        blast_cmd = 'blastn' if mode == 'dna' else 'blastx'
        output_file = self.raw_dir / f"{database}_{blast_cmd}_results.tsv"
        tool_name = f"BLAST-{mode.upper()}"

        if mode == 'dna':
            blast_cmd = 'blastn'
            cmd = [
                blast_cmd,
                '-query', str(self.input_fasta),
                '-db', db_path,
                '-out', str(output_file),
                '-outfmt',
                '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen',
                '-perc_identity', str(self.min_identity),
                #'-evalue', str(self.evalue),
                '-num_threads', str(self.threads)#,
               # '-max_target_seqs', str(self.max_target_seqs)
            ]
        else:
            blast_cmd = 'blastx'
            cmd = [
                blast_cmd,
                '-query', str(self.input_fasta),
                '-db', db_path,
                '-out', str(output_file),
                '-outfmt',
                '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen',
             #   '-evalue', str(self.evalue),
                '-num_threads', str(self.threads)#,
              #  '-max_target_seqs', str(self.max_target_seqs)
            ]

        success = self.run_command(cmd, f"{database} - {tool_name}")
        detected = set()
        if success:
            detected = self.parse_blast_results(output_file, database, tool_name)
            self.write_tool_stats(database, tool_name)
        return success, detected

    def run_diamond(self, db_path: str, database: str) -> Tuple[bool, Set[str]]:
        """Run DIAMOND protein search (blastx for DNA->protein)."""
        if not db_path:
            return False, set()

        output_file = self.raw_dir / f"{database}_diamond_results.tsv"
        tool_name = "DIAMOND"

        params = self.tool_sensitivity_params.get('diamond', None)
        sensitivity = params['sensitivity'] if params and 'sensitivity' in params else None


        cmd = [
            'diamond', 'blastx',
            '-q', str(self.input_fasta),
            '-d', db_path,
            '-o', str(output_file),
            '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen',
            '--header',
            '--id', str(self.min_identity),
            #'-e', str(self.evalue),
            '-p', str(self.threads)#,
            #'-k', '10'
        ]
        if sensitivity and sensitivity != 'default':
            cmd.append(sensitivity)

        success = self.run_command(cmd, f"{database} - {tool_name}")
        detected = set()
        if success:
            detected = self.parse_blast_results(output_file, database, tool_name)
            self.write_tool_stats(database, tool_name)
        return success, detected


    def run_bowtie2(self, db_path: str, database: str) -> Tuple[bool, Set[str]]:
        """Run Bowtie2 alignment (DNA mode) and output sorted BAM."""
        if not db_path:
            return False, set()

        sam_file = self.raw_dir / f"{database}_bowtie2_results.sam"
        bam_file = self.raw_dir / f"{database}_bowtie2_results.sorted.bam"
        summary_file = self.raw_dir / f"{database}_bowtie2_summary.txt"
        tool_name = "Bowtie2"

        params = self.tool_sensitivity_params.get('bowtie2', None)
        sensitivity = params['sensitivity'] if params and 'sensitivity' in params else None

        cmd = [
            'bowtie2',
            '-f',
            '-x', db_path,
            '-U', str(self.input_fasta),
            '-S', str(sam_file),
            '-p', str(self.threads),
            '--no-unal',
            '--met-file', str(summary_file)
        ]
        if sensitivity and sensitivity != 'default':
            cmd.append(sensitivity)

        success = self.run_command(cmd, f"{database} - {tool_name}")
        detected = set()
        if success:
            # Convert SAM to sorted BAM
            sort_cmd = [
                'samtools', 'sort',
                '-@', str(self.threads),
                '-o', str(bam_file),
                str(sam_file)
            ]
            sort_success = self.run_command(sort_cmd, f"{database} - samtools sort")
            if sort_success:
                detected = self.parse_sam_results(bam_file, database, tool_name)
                self.write_tool_stats(database, tool_name)
            success = success and sort_success
        return success, detected

    # def run_bwa(self, db_path: str, database: str) -> Tuple[bool, Set[str]]:
    #     """Run BWA alignment (DNA mode)."""
    #     if not db_path:
    #         return False, set()
    #
    #     sam_file = self.raw_dir / f"{database}_bwa_results.sam"
    #     tool_name = "BWA"
    #
    #     cmd = [
    #         'bwa', 'mem',
    #         '-t', str(self.threads),
    #         db_path,
    #         str(self.input_fasta)
    #     ]
    #
    #     # Run BWA and write output to SAM file
    #     try:
    #         with open(sam_file, 'w') as out_f:
    #             success = self.run_command(cmd + ['-o', str(sam_file)], f"{database} - {tool_name}")
    #     except Exception as e:
    #         self.logger.error(f"Error running BWA: {e}")
    #         return False, set()
    #
    #     detected = set()
    #     if success:
    #         detected = self.parse_sam_results(sam_file, database, tool_name)
    #         self.write_tool_stats(database, tool_name)
    #     return success, detected

    def run_bwa(self, db_path: str, database: str) -> Tuple[bool, Set[str]]:
        """Run BWA alignment (DNA mode) and output sorted BAM."""
        if not db_path:
            return False, set()

        sam_file = self.raw_dir / f"{database}_bwa_results.sam"
        bam_file = self.raw_dir / f"{database}_bwa_results.sorted.bam"
        tool_name = "BWA"

        cmd = [
            'bwa', 'mem',
            '-t', str(self.threads),
            db_path,
            str(self.input_fasta)
        ]

        # Run BWA and write output to SAM file
        try:
            success = self.run_command(cmd + ['-o', str(sam_file)], f"{database} - {tool_name}")
        except Exception as e:
            self.logger.error(f"Error running BWA: {e}")
            return False, set()

        detected = set()
        if success:
            # Convert SAM to sorted BAM
            sort_cmd = [
                'samtools', 'sort',
                '-@', str(self.threads),
                '-o', str(bam_file),
                str(sam_file)
            ]
            sort_success = self.run_command(sort_cmd, f"{database} - samtools sort")
            if sort_success:
                detected = self.parse_sam_results(bam_file, database, tool_name)
                self.write_tool_stats(database, tool_name)
            success = success and sort_success
        return success, detected


    def run_minimap2(self, db_path: str, database: str, preset: str = 'sr') -> Tuple[bool, Set[str]]:
        """Run Minimap2 alignment and output sorted BAM."""
        if not db_path:
            return False, set()

        sam_file = self.raw_dir / f"{database}_minimap2_results.sam"
        bam_file = self.raw_dir / f"{database}_minimap2_results.sorted.bam"
        tool_name = "Minimap2"

        cmd = [
            'minimap2',
            '-x', preset,
            '-t', str(self.threads),
            '-a',
            db_path,
            str(self.input_fasta),
            '-o', str(sam_file)
        ]

        success = self.run_command(cmd, f"{database} - {tool_name}")
        detected = set()
        if success:
            # Convert SAM to sorted BAM
            sort_cmd = [
                'samtools', 'sort',
                '-@', str(self.threads),
                '-o', str(bam_file),
                str(sam_file)
            ]
            sort_success = self.run_command(sort_cmd, f"{database} - samtools sort")
            if sort_success:
                detected = self.parse_sam_results(bam_file, database, tool_name)
                self.write_tool_stats(database, tool_name)
            success = success and sort_success
        return success, detected

    def run_hmmer(self, db_path: str, database: str, mode: str) -> Tuple[bool, Set[str]]:
        """Run HMMER profile search."""
        if not db_path:
            return False, set()

        hmmer_cmd = 'nhmmer' if mode == 'dna' else 'hmmsearch'
        output_file = self.raw_dir / f"{database}_{hmmer_cmd}_results.tbl"
        domtbl_file = self.raw_dir / f"{database}_{hmmer_cmd}_domtbl.txt"
        tool_name = f"HMMER-{mode.upper()}"

        cmd = [
            hmmer_cmd,
            '--tblout', str(output_file),
            '--domtblout', str(domtbl_file),
            '-E', str(self.evalue),
            '--cpu', str(self.threads),
            db_path,
            str(self.input_fasta)
        ]

        success = self.run_command(cmd, f"{database} - {tool_name}")
        detected = set()
        if success:
            detected = self.parse_hmmer_results(output_file, database, tool_name)
            self.write_tool_stats(database, tool_name)
        return success, detected

    def generate_detection_matrix(self, database: str):
        """Generate TSV matrix of gene detections across tools."""
        output_file = self.output_dir / f"{database}_detection_matrix.tsv"

        # Get all tools that were run for this database
        all_tools = set()
        for gene_detections in self.detections[database].values():
            all_tools.update(gene_detections.keys())

        if not all_tools:
            self.logger.warning(f"No detections found for {database}")
            return

        all_tools = sorted(all_tools)

        # Write matrix
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            # Header
            header = ['Gene'] + all_tools + ['Total_Detections']
            writer.writerow(header)

            # Only include genes with at least one detection and sort
            if database == 'card':
                def get_last_segment(gene_name):
                    return gene_name.split('|')[-1] if '|' in gene_name else gene_name

                genes = [gene for gene in sorted(
                    self.detections[database].keys(),
                    key=get_last_segment
                ) if any(self.detections[database][gene][tool] for tool in all_tools)]
            else:
                genes = [gene for gene in sorted(self.detections[database].keys())
                         if any(self.detections[database][gene][tool] for tool in all_tools)]

            for gene in genes:
                row = [gene]
                detections = self.detections[database][gene]

                for tool in all_tools:
                    row.append('1' if detections[tool] else '0')

                # Count total detections
                total = sum(1 for tool in all_tools if detections[tool])
                row.append(str(total))

                writer.writerow(row)

        self.logger.info(f"Generated detection matrix: {output_file}")
        self.logger.info(f"  Total genes detected: {len(genes)}")
        self.logger.info(f"  Tools used: {len(all_tools)}")


    def run_pipeline(self,options):

        """Run all configured tools on both databases."""
        self.logger.info("=" * 70)
        self.logger.info("AMRfíor - The AMR Gene Detection tool: " + AMRFIOR_VERSION)
        self.logger.info("=" * 70)
        self.logger.info(f"Input file: {self.input_fasta}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Threads: {self.threads}")
       # self.logger.info(f"E-value threshold: {self.evalue}")
        self.logger.info(f"Min coverage: {self.min_coverage}%")
        self.logger.info(f"Min identity: {self.min_identity}%")
        self.logger.info(f"Run DNA mode: {self.run_dna}")
        self.logger.info(f"Run Protein mode: {self.run_protein}")
        params_str = ", ".join(
            f"{tool}: {params}" for tool, params in self.tool_sensitivity_params.items()
        ) if self.tool_sensitivity_params else "None"
        self.logger.info(f"Sensitivity parameters: {options.sensitivity} - {params_str}")
        self.logger.info("=" * 70)
        results = {'resfinder': {}, 'card': {}}

        # Process ResFinder database
        if self.resfinder_dbs:
            self.logger.info("\n### Processing ResFinder Database ###")

            if self.run_dna and self.resfinder_dbs.get('blastn'):
                results['resfinder']['BLASTn-DNA'] = self.run_blast(
                    self.resfinder_dbs['blastn'], 'resfinder', 'dna')

            if self.run_protein and self.resfinder_dbs.get('blastp'):
                results['resfinder']['BLASTp-AA'] = self.run_blast(
                    self.resfinder_dbs['blastp'], 'resfinder', 'protein')

            if self.run_protein and self.resfinder_dbs.get('diamond'):
                results['resfinder']['DIAMOND-AA'] = self.run_diamond(
                    self.resfinder_dbs['diamond'], 'resfinder')

            if self.run_dna and self.resfinder_dbs.get('bowtie2'):
                results['resfinder']['Bowtie2-DNA'] = self.run_bowtie2(
                    self.resfinder_dbs['bowtie2'], 'resfinder')

            if self.run_dna and self.resfinder_dbs.get('bwa'):
                results['resfinder']['BWA-DNA'] = self.run_bwa(
                    self.resfinder_dbs['bwa'], 'resfinder')

            if self.resfinder_dbs.get('minimap2'):
                results['resfinder']['Minimap2-DNA'] = self.run_minimap2(
                    self.resfinder_dbs['minimap2'], 'resfinder', options.minimap2_preset)

            # if self.run_dna and self.resfinder_dbs.get('hmmer_dna'):
            #     results['resfinder']['HMMER-DNA'] = self.run_hmmer(
            #         self.resfinder_dbs['hmmer_dna'], 'resfinder', 'dna')
            #
            # if self.run_protein and self.resfinder_dbs.get('hmmer_protein'):
            #     results['resfinder']['HMMER-PROTEIN'] = self.run_hmmer(
            #         self.resfinder_dbs['hmmer_protein'], 'resfinder', 'protein')

            self.generate_detection_matrix('resfinder')

        # Process CARD database
        if self.card_dbs:
            self.logger.info("\n### Processing CARD Database ###")

            if self.run_dna and self.card_dbs.get('blastn'):
                results['card']['BLASTn-DNA'] = self.run_blast(
                    self.card_dbs['blastn'], 'card', 'dna')

            if self.run_protein and self.card_dbs.get('blastp'):
                results['card']['BLASTp-AA'] = self.run_blast(
                    self.card_dbs['blastp'], 'card', 'protein')

            if self.run_protein and self.card_dbs.get('diamond'):
                results['card']['DIAMOND-AA'] = self.run_diamond(
                    self.card_dbs['diamond'], 'card')

            if self.run_dna and self.card_dbs.get('bowtie2'):
                results['card']['Bowtie2-DNA'] = self.run_bowtie2(
                    self.card_dbs['bowtie2'], 'card')

            if self.run_dna and self.card_dbs.get('bwa'):
                results['card']['BWA-DNA'] = self.run_bwa(
                    self.card_dbs['bwa'], 'card')

            if self.card_dbs.get('minimap2'):
                results['card']['Minimap2-DNA'] = self.run_minimap2(
                    self.card_dbs['minimap2'], 'card', options.minimap2_preset)

            # if self.run_dna and self.card_dbs.get('hmmer_dna'):
            #     results['card']['HMMER-DNA'] = self.run_hmmer(
            #         self.card_dbs['hmmer_dna'], 'card', 'dna')
            #
            # if self.run_protein and self.card_dbs.get('hmmer_protein'):
            #     results['card']['HMMER-PROTEIN'] = self.run_hmmer(
            #         self.card_dbs['hmmer_protein'], 'card', 'protein')

            self.generate_detection_matrix('card')

        # Final summary
        self.logger.info("\n" + "=" * 70)
        self.logger.info("PIPELINE SUMMARY")
        self.logger.info("=" * 70)

        for db_name in ['resfinder', 'card']:
            if results[db_name]:
                self.logger.info(f"\n{db_name.upper()}:")
                for tool, (success, genes) in results[db_name].items():
                    status = "✓" if success else "✗"
                    gene_count = len(genes) if success else 0
                    self.logger.info(f"  {status} {tool:.<30} {gene_count} genes detected")

        self.logger.info("=" * 70)
        self.logger.info(f"Detection matrices saved to: {self.output_dir}")
        self.logger.info(f"Tool statistics saved to: {self.stats_dir}")
        self.logger.info(f"Raw outputs saved to: {self.raw_dir}")
        self.logger.info("=" * 70)

        return results