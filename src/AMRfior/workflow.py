import subprocess, sys
import csv
from collections import defaultdict
from pathlib import Path
import logging
from datetime import datetime
from typing import Any, Dict, List, Set, Tuple
import re


try:
    from .gene_stats import GeneStats
    from .constants import *
except (ModuleNotFoundError, ImportError) as error:
    from gene_stats import GeneStats
    from constants import *

class AMRWorkflow:
    """Orchestrates multiple alignment tools for AMR gene detection."""

    def __init__(self, input_fasta: str, input_fastq: str, output_dir: str,
                 databases: Dict[str, Dict[str, str]],
                 threads: int = 4, #max_target_seqs: int = 100,
                 tool_sensitivity_params: Dict[str, Dict[str, Any]] = None,
                 #evalue: float = 1e-10,
                 detection_min_coverage: float = 80.0, detection_min_identity: float = 80.0,
                 detection_min_base_depth: float = 1.0,
                 detection_min_num_reads: int = 1,
                 query_min_coverage: float = 50.0,  # NEW: Query coverage threshold

                 run_dna: bool = True, run_protein: bool = True,
                 sequence_type: str = 'Single-FASTA',
                 report_fasta: str = None,
                 no_cleanup: bool = False,
                 verbose: bool = False):
        ### Handle input FASTA and FASTQ
        if input_fasta is not None:
            self.input_fasta = Path(input_fasta)
        if input_fastq is None:
            self.input_fastq = None
            self.input_fastq_is_paired = False
        elif isinstance(input_fastq, (tuple, list)):
            if len(input_fastq) != 2:
                raise ValueError("`input_fastq` tuple/list must contain exactly two paths for paired-end reads")
            self.input_fastq = (Path(input_fastq[0]), Path(input_fastq[1]))
            self.input_fastq_is_paired = True
        elif isinstance(input_fastq, str) and ',' in input_fastq:
            parts = [p.strip() for p in input_fastq.split(',') if p.strip()]
            if len(parts) != 2:
                raise ValueError(
                    "`input_fastq` comma-separated string must contain exactly two paths for paired-end reads")
            self.input_fastq = (Path(parts[0]), Path(parts[1]))
            self.input_fastq_is_paired = True
        else:
            self.input_fastq = Path(input_fastq)
            self.input_fastq_is_paired = False
        ###
        self.output_dir = Path(output_dir)
        self.databases = databases
        self.threads = threads
      #  self.max_target_seqs = max_target_seqs
        self.tool_sensitivity_params = tool_sensitivity_params
       # self.evalue = evalue
        self.detection_min_coverage = detection_min_coverage
        self.detection_min_identity = detection_min_identity
        self.detection_min_base_depth = detection_min_base_depth
        self.detection_min_num_reads = detection_min_num_reads
        self.query_min_coverage = query_min_coverage

        self.run_dna = run_dna
        self.run_protein = run_protein
        self.sequence_type = sequence_type
        self.report_fasta = report_fasta  # 'All', 'Detected', or None

        # Create output directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.raw_dir = self.output_dir / "raw_outputs"
        self.raw_dir.mkdir(exist_ok=True)
        self.stats_dir = self.output_dir / "tool_stats"
        self.stats_dir.mkdir(exist_ok=True)
        if self.report_fasta[0] != None:
            self.fasta_dir = self.output_dir / "fasta_outputs"
            self.fasta_dir.mkdir(exist_ok=True)

        # misc
        self.no_cleanup = no_cleanup
        self.verbose = verbose


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
            db_name: defaultdict(lambda: defaultdict(bool))
            for db_name in self.databases.keys()
        }

        # Store detailed statistics: {database: {tool: {gene: GeneStats}}}
        self.gene_stats = {
            db_name: defaultdict(lambda: defaultdict(GeneStats))
            for db_name in self.databases.keys()
        }

    def check_gzip(self,fasta_path_str):
        self.logger.info(f"Checking if input FASTA ` {self.input_fasta} ` is gzipped and not broken...")
        # Check gzip integrity before attempting to stream; prefer `gzip -t` then fall back to Python check
        gz_ok = True
        if str(fasta_path_str).endswith(('.gz', '.gzip')):
            try:
                # Use system gzip test if available (faster for large files)
                res = subprocess.run(['gzip', '-t', fasta_path_str], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if res.returncode != 0:
                    gz_ok = False
            except FileNotFoundError:
                # Fallback: try opening a byte to ensure file is a valid gzip
                try:
                    import gzip as _gzip
                    with _gzip.open(fasta_path_str, 'rb') as fh:
                        fh.read(1)
                except (OSError, EOFError):
                    gz_ok = False
        else:
            gz_ok = False

        if not gz_ok:
            self.logger.error(f"Input FASTA ` {fasta_path_str} ` appears to be a broken or invalid gzip file.")
            return False
        else:
            return True

    def run_command_gzip(self, cmd: List[str], tool_name: str, fasta_path_str: str) -> bool:
        # If gzipped, stream decompressed FASTA into BLAST stdin to avoid writing a temp file
        self.logger.info(f"Running {tool_name}...")
        self.logger.info(f"Parameters for {tool_name}: {' '.join(str(arg) for arg in cmd)}")
        #self.logger.info(f"Parameters for {tool_name}: {' '.join(cmd)}")
        #self.logger.debug(f"Command: {' '.join(cmd)}")
        self.logger.debug(f"Command: {' '.join(str(arg) for arg in cmd)}")
        is_gzip_valid = self.check_gzip(fasta_path_str)
        if not is_gzip_valid:
            return False
        elif is_gzip_valid:

            self.logger.info(f"Piping decompressed ` {self.input_fasta} ` to BLAST ({tool_name}) to avoid disk IO")
            try:
                gzip_proc = subprocess.Popen(['gzip', '-dc', fasta_path_str], stdout=subprocess.PIPE)
                # Run BLAST reading from stdin
                result = subprocess.run(
                    cmd,
                    stdin=gzip_proc.stdout,
                    check=True,
                    capture_output=True,
                    text=True
                )
                self.logger.info(f"{tool_name} completed successfully")
                # Close gzip stdout to allow gzip to receive SIGPIPE if BLAST exits early
                gzip_proc.stdout.close()
                gzip_proc.wait()
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
            except Exception as e:
                self.logger.error(f"Error running {tool_name} with piped input: {e}")
                return False

    def run_command(self, cmd: List[str], tool_name: str) -> bool:
        # Run a tool and log the results.
        self.logger.info(f"Running {tool_name}...")
        self.logger.info(f"Parameters for {tool_name}: {' '.join(str(arg) for arg in cmd)}")
        #self.logger.info(f"Parameters for {tool_name}: {' '.join(cmd)}")
        #self.logger.debug(f"Command: {' '.join(cmd)}")
        self.logger.debug(f"Command: {' '.join(str(arg) for arg in cmd)}")

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


    def _write_fasta_outputs(self, database: str, tool_name: str, detected_genes: Set[str],
                             gene_reads: dict, all_reads: dict):
        # Method to write FASTA files for mapped reads.
        def sanitise_gene_name(gene: str) -> str:
            safe_gene = gene.replace('|', '_').replace('/', '_').replace(':','_').replace('-','_')
            if not hasattr(self, 'gene_name_changes'):
                self.gene_name_changes = []
            self.gene_name_changes.append((gene, safe_gene))
            return safe_gene


        if self.report_fasta == 'all':
            if getattr(self, "verbose", True):
                self.logger.info(f"Writing FASTA files for all mapped reads in {database}...")
            for gene, read_types in gene_reads.items():
                read_names = set(read_types['all'])  # Use set to avoid duplicates
                if not read_names:
                    continue

                safe_gene = sanitise_gene_name(gene)
                fasta_path = self.fasta_dir / f"{database}_{tool_name}_{safe_gene}_reads.fasta"

                with open(fasta_path, "w") as fasta_out:
                    count = 0
                    for read_name in read_names:
                        if read_name in all_reads:
                            seq = all_reads[read_name]
                            fasta_out.write(f">{read_name}\n{seq}\n")
                            count += 1
                if getattr(self, "verbose", True):
                    self.logger.info(f"  FASTA file: {fasta_path} ({count} reads)")

        elif self.report_fasta == 'detected':
            if getattr(self, "verbose", True):
                self.logger.info(
                    f"Writing FASTA files for threshold-passing reads mapped to detected genes in {database}...")
            for gene in detected_genes:
                read_names = set(gene_reads[gene].get('passing', []))
                if not read_names:
                    continue

                safe_gene = sanitise_gene_name(gene)
                fasta_path = self.fasta_dir / f"{database}_{tool_name}_{safe_gene}_reads.fasta"

                with open(fasta_path, "w") as fasta_out:
                    count = 0
                    for read_name in read_names:
                        if read_name in all_reads:
                            seq = all_reads[read_name]
                            fasta_out.write(f">{read_name}\n{seq}\n")
                            count += 1
                if getattr(self, "verbose", True):
                    self.logger.info(f"  FASTA file: {fasta_path} ({count} reads)")

        elif self.report_fasta == 'detected-all':
            if getattr(self, "verbose", True):
                self.logger.info(f"Writing FASTA files for all reads mapped to detected genes in {database}...")
            for gene in detected_genes:
                read_names = set(gene_reads[gene].get('all', []))
                if not read_names:
                    continue

                safe_gene = sanitise_gene_name(gene)
                fasta_path = self.fasta_dir / f"{database}_{tool_name}_{safe_gene}_reads.fasta"

                with open(fasta_path, "w") as fasta_out:
                    count = 0
                    for read_name in read_names:
                        if read_name in all_reads:
                            seq = all_reads[read_name]
                            fasta_out.write(f">{read_name}\n{seq}\n")
                            count += 1
                if getattr(self, "verbose", True):
                    self.logger.info(f"  FASTA file: {fasta_path} ({count} reads)")




    def run_blast(self, db_path: str, database: str, mode: str) -> Tuple[bool, Set[str]]:
        """Run BLAST in DNA (blastn) or protein (blastx) mode.
            If input FASTA is gzipped, stream decompressed data to BLAST via stdin
            (uses `-query -`) to avoid creating a temporary uncompressed file."""
        blast_cmd = 'blastn' if mode == 'dna' else 'blastx'
        tool_name = 'BLASTn' if mode == 'dna' else 'BLASTx'
        # Extract the path from the dictionary if `database` is a dict
        if isinstance(database, dict):
            database = database.get(blast_cmd)
        # Ensure `db_path` is a string
        if isinstance(db_path, dict):
            db_path = db_path.get(blast_cmd)

        output_file = self.raw_dir / f"{database}_{blast_cmd}_results.tsv"
        #tool_name = f"BLAST-{mode.upper()}"

        # Common outfmt used for both modes
        outfmt_fields = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

        # Determine if input is gzipped
        fasta_path_str = str(self.input_fasta)
        gz_input = fasta_path_str.endswith(('.gz', '.gzip'))

        if mode == 'dna':
            blast_cmd = 'blastn'
            query_arg = '-' if gz_input else fasta_path_str
            cmd = [
                blast_cmd,
                '-query', query_arg,
                '-db', db_path,
                '-out', str(output_file),
                '-outfmt', outfmt_fields,
                '-num_threads', str(self.threads)
            ]
        elif mode == 'protein':
            blast_cmd = 'blastx'
            query_arg = '-' if gz_input else fasta_path_str
            cmd = [
                blast_cmd,
                '-query', query_arg,
                '-db', db_path,
                '-out', str(output_file),
                '-outfmt', outfmt_fields,
                '-num_threads', str(self.threads)
            ]
        else:
            self.logger.error(f"Invalid BLAST'ing' mode: {mode}")
            return False, set()

        success = False

        if gz_input:
            # Gzipped case - stream decompressed FASTA to BLAST stdin
            success = self.run_command_gzip(cmd, f"{database} - {tool_name}", fasta_path_str)
        else:
            # Non-gz case
            success = self.run_command(cmd, f"{database} - {tool_name}")

        detected = set()
        if success:
            detected, gene_reads = self.parse_blast_results(output_file, database, tool_name)
            self.write_tool_stats(database, tool_name, gene_reads)
        return success, detected

    def run_diamond(self, db_path: str, database: str) -> Tuple[bool, Set[str]]:
        # Run DIAMOND protein search (blastx for DNA->protein).
        # Extract the path from the dictionary if `database` is a dict
        if isinstance(database, dict):
            database = database.get('diamond')
        # Ensure `db_path` is a string
        if isinstance(db_path, dict):
            db_path = db_path.get('diamond')

        if not isinstance(db_path, str):
            raise ValueError(f"Invalid database path: {db_path}")

        output_file = self.raw_dir / f"{database}_diamond_results.tsv"
        tool_name = "DIAMOND"

        fasta_path_str = str(self.input_fasta)
        gz_input = fasta_path_str.endswith(('.gz', '.gzip'))
        if gz_input and not self.check_gzip(fasta_path_str): # If input is gzipped, check integrity first
            return False
        else: # Run DIAMOND normally
            params = self.tool_sensitivity_params.get('diamond', None)
            sensitivity = params['sensitivity'] if params and 'sensitivity' in params else None

            cmd = [
                'diamond', 'blastx',
                '-q', str(self.input_fasta),
                '-d', db_path,
                '-o', str(output_file),
                '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen',
                '--header',
                #'--id', str(self.min_identity),
                #'-e', str(self.evalue),
                '-p', str(self.threads)#,
                #'-k', '10'
            ]
            if sensitivity and sensitivity != 'default':
                cmd.append(sensitivity)

            success = self.run_command(cmd, f"{database} - {tool_name}")
            detected = set()
            if success:
                detected, gene_reads = self.parse_blast_results(output_file, database, tool_name)
                self.write_tool_stats(database, tool_name, gene_reads)
            return success, detected


    def run_bowtie2(self, db_path: str, database: str) -> Tuple[bool, Set[str]]:
        # Run Bowtie2 alignment (DNA mode) and output sorted BAM.
        # Extract the path from the dictionary if `database` is a dict
        if isinstance(database, dict):
            database = database.get('bowtie2')
        # Ensure `db_path` is a string
        if isinstance(db_path, dict):
            db_path = db_path.get('bowtie2')

        sam_file = self.raw_dir / f"{database}_bowtie2_results.sam"
        bam_file = self.raw_dir / f"{database}_bowtie2_results.bam"
        sorted_bam_file = self.raw_dir / f"{database}_bowtie2_results_sorted.bam"
        summary_file = self.raw_dir / f"{database}_bowtie2_summary.txt"
        tool_name = "Bowtie2"

        if self.sequence_type == 'Single-FASTA':
            flags = ['-f', '-U', str(self.input_fasta)]
        elif self.sequence_type == 'Paired-FASTQ':
            flags = ['-1', str(self.input_fastq[0]), '-2', str(self.input_fastq[1])]
        else:
            flags = []

        params = self.tool_sensitivity_params.get('bowtie2', None)
        sensitivity = params['sensitivity'] if params and 'sensitivity' in params else None

        cmd = [
            'bowtie2',
              ] + flags + [
            '-x', db_path,
            '-S', str(sam_file),
            '-p', str(self.threads),
            '--no-unal',
            '--met-file', str(summary_file)
        ]
        if sensitivity and sensitivity != 'default':
            cmd.append(sensitivity)

        success = self.run_command(cmd, f"{database} - {tool_name}")
        if not success:
            return False, set()

        # Convert SAM to BAM
        sam_to_bam_cmd = ['samtools', 'view', '-bS', str(sam_file), '-o', str(bam_file)]
        if not self.run_command(sam_to_bam_cmd, f"{database} - SAM to BAM conversion"):
            return False, set()

        # Sort BAM
        sort_cmd = ['samtools', 'sort', str(bam_file), '-o', str(sorted_bam_file)]
        if not self.run_command(sort_cmd, f"{database} - BAM sorting"):
            return False, set()

        # Index BAM (optional but recommended)
        index_cmd = ['samtools', 'index', str(sorted_bam_file)]
        self.run_command(index_cmd, f"{database} - BAM indexing")

        # Parse results
        detected, gene_reads = self.parse_bam_results(sorted_bam_file, database, tool_name)
        self.write_tool_stats(database, tool_name, gene_reads)


        return success, detected



    def run_bwa(self, db_path: str, database: str) -> Tuple[bool, Set[str]]:
        # Run BWA alignment (DNA mode) and output sorted BAM.
        if isinstance(database, dict):
            database = database.get('bwa')
        # Ensure `db_path` is a string
        if isinstance(db_path, dict):
            db_path = db_path.get('bwa')
        sam_file = self.raw_dir / f"{database}_bwa_results.sam"
        bam_file = self.raw_dir / f"{database}_bwa_results.bam"
        sorted_bam = self.raw_dir / f"{database}_bwa_results_sorted.bam"
        tool_name = "BWA"

        if self.sequence_type == 'Single-FASTA':
            flags = [ str(self.input_fasta)]
        elif self.sequence_type == 'Paired-FASTQ':
            flags = [ str(self.input_fastq[0]), str(self.input_fastq[1])]
        else:
            flags = []

        cmd = [
            'bwa', 'mem',
            '-t', str(self.threads),
            db_path,
            ] + flags + [
        ]

        # Run BWA and write output to SAM file
        try:
            success = self.run_command(cmd + ['-o', str(sam_file)], f"{database} - {tool_name}")
        except Exception as e:
            self.logger.error(f"Error running BWA: {e}")
            return False, set()

        if not success:
            return False, set()

        # Convert SAM to BAM
        sam_to_bam_cmd = ['samtools', 'view', '-bS', str(sam_file), '-o', str(bam_file)]
        if not self.run_command(sam_to_bam_cmd, f"{database} - SAM to BAM conversion"):
            return False, set()

        # Sort BAM
        sort_cmd = ['samtools', 'sort', str(bam_file), '-o', str(sorted_bam)]
        if not self.run_command(sort_cmd, f"{database} - BAM sorting"):
            return False, set()

        # Index BAM (optional but recommended)
        index_cmd = ['samtools', 'index', str(sorted_bam)]
        self.run_command(index_cmd, f"{database} - BAM indexing")

        # Parse results
        detected, gene_reads = self.parse_bam_results(sorted_bam, database, tool_name)
        self.write_tool_stats(database, tool_name, gene_reads)

        return success, detected


    def run_minimap2(self, db_path: str, database: str, preset: str = 'sr') -> Tuple[bool, Set[str]]:
        # Run Minimap2 alignment and output sorted BAM.
        if isinstance(database, dict):
            database = database.get('minimap2')
        # Ensure `db_path` is a string
        if isinstance(db_path, dict):
            db_path = db_path.get('minimap2')

        sam_file = self.raw_dir / f"{database}_minimap2_results.sam"
        bam_file = self.raw_dir / f"{database}_minimap2_results.bam"
        sorted_bam = self.raw_dir / f"{database}_minimap2_results_sorted.bam"
        tool_name = "Minimap2"


        if self.sequence_type == 'Single-FASTA':
            flags = [ str(self.input_fasta)]
        elif self.sequence_type == 'Paired-FASTQ':
            flags = [ str(self.input_fastq[0]), str(self.input_fastq[1])]
        else:
            flags = []

        cmd = [
            'minimap2',
            '-x', preset,
            '-t', str(self.threads),
            '-a',
            db_path,
            ] + flags + [
            '-o', str(sam_file)
        ]

        try:
            success = self.run_command(cmd, f"{database} - {tool_name}")
        except Exception as e:
            self.logger.error(f"Error running Minimap2: {e}")
            return False, set()

        if not success:
            return False, set()

        # Convert SAM to BAM
        sam_to_bam_cmd = ['samtools', 'view', '-bS', str(sam_file), '-o', str(bam_file)]
        if not self.run_command(sam_to_bam_cmd, f"{database} - SAM to BAM conversion"):
            return False, set()

        # Sort BAM
        sort_cmd = ['samtools', 'sort', str(bam_file), '-o', str(sorted_bam)]
        if not self.run_command(sort_cmd, f"{database} - BAM sorting"):
            return False, set()

        # Index BAM (optional but recommended)
        index_cmd = ['samtools', 'index', str(sorted_bam)]
        self.run_command(index_cmd, f"{database} - BAM indexing")

        # Parse results
        detected, gene_reads = self.parse_bam_results(sorted_bam, database, tool_name)
        self.write_tool_stats(database, tool_name, gene_reads)

        return success, detected




    def run_hmmer(self, db_path: str, database: str, mode: str) -> Tuple[bool, Set[str]]:
        # Run HMMER search.
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

    def parse_blast_results(self, output_file: Path, database: str, tool_name: str) -> Set[str]:
        """Parse BLAST/DIAMOND tabular output and extract genes meeting thresholds.
        Detection logic:
        - Only sequences with identity >= detection_min_identity are considered
        - Track which positions on the subject/gene are covered by alignments
        - Gene is detected if combined coverage of gene >= detection_min_coverage
        """
        detected_genes = set()
        gene_lengths = {}  # Store gene lengths
        gene_reads = defaultdict(lambda: {'passing': [], 'all': []})  # Track all reads per gene

        if not output_file.exists():
            return detected_genes

        # Load all reads from input FASTA for later FASTA output (cached) - Should only load once
        if not hasattr(self, 'all_reads'):
            self.all_reads = {}
            try:
                import gzip
                fasta_path = str(self.input_fasta)
                # Open gzipped or plain FASTA transparently
                if fasta_path.endswith(('.gz', '.gzip')):
                    fasta_handle = gzip.open(fasta_path, 'rt')
                else:
                    fasta_handle = open(fasta_path, 'r')
                with fasta_handle as fasta_file:
                    read_name = None
                    seq_lines = []
                    for line in fasta_file:
                        # If gzip returns bytes unexpectedly, decode
                        if isinstance(line, bytes):
                            line = line.decode('utf-8')
                        if line.startswith('>'):
                            if read_name in self.all_reads:
                                self.logger.error(f"Warning: Duplicate read name found in FASTA: {read_name}")
                            if read_name and seq_lines:
                                self.all_reads[read_name] = ''.join(seq_lines)
                            read_name = line[1:].strip()
                            seq_lines = []
                        else:
                            seq_lines.append(line.strip())
                    if read_name and seq_lines:
                        self.all_reads[read_name] = ''.join(seq_lines)
            except Exception as e:
                self.logger.error(f"Error reading FASTA file: {e}")
        all_reads = self.all_reads
        mapped_reads = 0
        passing_reads = 0
        try:
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) < 13:
                        continue

                    read_name = fields[0]  # qseqid
                    gene = fields[1]  # sseqid
                    identity = float(fields[2])  # pident
                    qstart = int(fields[6])  # query start
                    qend = int(fields[7])  # query end
                    sstart = int(fields[8])  # subject start
                    send = int(fields[9])  # subject end
                    qlen = int(fields[12])  # query length (added to output format)
                    slen = int(fields[13])  # subject length (added to output format)

                    # Store gene length
                    if gene in gene_lengths:
                        gene_lengths[gene] = max(gene_lengths[gene], slen)
                    else:
                        gene_lengths[gene] = slen

                    # Only process sequences meeting identity and query coverage thresholds - variable redundancy for clarity
                    query_length = qlen
                    query_coverage = ((abs(qend - qstart) + 1) / query_length) * 100 if query_length else 0

                    # Track all reads mapping to this gene
                    gene_reads[gene]['all'].append(read_name)
                    mapped_reads +=1

                    if identity >= self.detection_min_identity and query_coverage >= self.query_min_coverage:
                        # Initialise stats if first hit for this gene
                        if gene not in self.gene_stats[database][tool_name]:
                            self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)

                        # Add hit to statistics
                        self.gene_stats[database][tool_name][gene].add_hit(
                            sstart, send, identity, gene_lengths[gene]
                        )

                        # Track reads that pass thresholds
                        gene_reads[gene]['passing'].append(read_name)
                        passing_reads +=1

        except Exception as e:
            self.logger.error(f"Error parsing {output_file}: {e}")

        # Finalise statistics and determine detection based on gene coverage
        for gene in self.gene_stats[database][tool_name]:
            stats = self.gene_stats[database][tool_name][gene]
            stats.finalise()

            # Gene is detected if gene meets thresholds
            if (stats.gene_coverage >= self.detection_min_coverage and
                    stats.base_depth_hit >= self.detection_min_base_depth and
                    stats.num_sequences >= self.detection_min_num_reads):
                detected_genes.add(gene)
                self.detections[database][gene][tool_name] = True

        self.logger.info(f"Total reads processed in {database.upper()} using {tool_name}: {len(all_reads)}")
        ## Need a fix for this
        self.logger.info(f"Reads that returned a hit in {database.upper()} using {tool_name}: {mapped_reads}")
        self.logger.info(f"Reads passing thresholds in {database.upper()} using {tool_name}: {passing_reads}")
        self.logger.info(f"Detected {len(detected_genes)} genes in {database.upper()} using {tool_name}")

        # Output FASTA files of reads mapping to genes
        if self.report_fasta and detected_genes:
            self._write_fasta_outputs(database, tool_name, detected_genes, gene_reads, all_reads)

        return detected_genes, gene_reads

    def parse_bam_results(self, bam_file: Path, database: str, tool_name: str) -> Set[str]:
        """Parse BAM file from Bowtie2 and extract genes meeting thresholds.

        Detection logic:
        - Only seqs with identity >= detection_min_identity and query coverage >= query-min-coverage are considered
        - Track which positions on the subject/gene are covered by ALIGNED bases only
        - Gene is detected if combined coverage of gene >= detection_min_coverage

        Identity calculation 'matches' DIAMOND/BLAST: (matches / alignment_length) * 100
        Coverage tracks only M/=/X operations (actual aligned bases on reference)
        """
        detected_genes = set()
        if not bam_file.exists():
            self.logger.error(f"BAM file not found: {bam_file}")
            return detected_genes
        gene_lengths = {}  # Store gene lengths from BAM header
        gene_reads = defaultdict(lambda: {'passing': [], 'all': []})  # Track all reads per gene
        all_reads = {}
        mapped_reads = 0
        passing_reads = 0

        try:
            proc = subprocess.Popen(['samtools', 'view', '-h', str(bam_file)],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')

            for line in proc.stdout:
                if line.startswith('@SQ'):
                    # @SQ\tSN:ref\tLN:1234
                    parts = line.strip().split('\t')
                    sn = None
                    ln = None
                    for p in parts:
                        if p.startswith('SN:'):
                            sn = p.split(':', 1)[1]
                        if p.startswith('LN:'):
                            ln = int(p.split(':', 1)[1])
                    if sn and ln:
                        gene_lengths[sn] = ln
                    continue
                if line.startswith('@'):
                    continue  # other headers

                fields = line.rstrip('\n').split('\t')
                if len(fields) < 11:
                    continue

                read_name = fields[0]
                flag = int(fields[1])
                # skip unmapped
                if flag & 0x4:
                    continue

                gene = fields[2]
                try:
                    ref_start = int(fields[3]) - 1  # convert to 0-based
                except ValueError:
                    ref_start = 0
                cigar = fields[5]
                seq = fields[9]
                # store read sequence if available
                if read_name not in all_reads and seq and seq != '*':
                    all_reads[read_name] = seq

                gene_len = gene_lengths.get(gene, 0)
                gene_reads[gene]['all'].append(read_name)
                mapped_reads +=1

                # try to get NM tag from optional fields
                nm = 0
                for opt in fields[11:]:
                    if opt.startswith('NM:i:'):
                        try:
                            nm = int(opt.split(':')[-1])
                        except ValueError:
                            nm = 0
                        break

                # parse CIGAR and compute aligned positions & alignment length
                ref_pos = ref_start
                aligned_positions = set()
                alignment_length = 0

                for count_str, op in cigar_re.findall(cigar):
                    length = int(count_str)
                    if op in ('M', '=', 'X'):
                        aligned_positions.update(range(ref_pos, ref_pos + length))
                        ref_pos += length
                        alignment_length += length
                    elif op == 'I':  # insertion to reference
                        alignment_length += length
                    elif op == 'D':  # deletion from reference
                        ref_pos += length
                        alignment_length += length
                    elif op == 'N':
                        ref_pos += length
                    elif op in ('S', 'H'):
                        # soft/hard clip - do not consume reference (H doesn't appear in SEQ)
                        pass

                if alignment_length == 0:
                    identity = 0.0
                else:
                    matches = max(0, alignment_length - nm)
                    identity = (matches / alignment_length) * 100.0

                query_length = len(seq) if seq and seq != '*' else 0
                query_coverage = (len(aligned_positions) / query_length) * 100.0 if query_length > 0 else 0.0

                if identity >= self.detection_min_identity and query_coverage >= self.query_min_coverage:
                    if gene not in self.gene_stats[database][tool_name]:
                        self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)
                    self.gene_stats[database][tool_name][gene].add_positions(aligned_positions, identity, gene_len)
                    gene_reads[gene]['passing'].append(read_name)
                    passing_reads +=1

            proc.stdout.close()
            proc.wait()
            # optionally capture and log stderr
            stderr = proc.stderr.read() if proc.stderr else ''
            if stderr:
                self.logger.debug(f"samtools stderr: {stderr}")
            if proc.stderr:
                proc.stderr.close()

        except Exception as e:
            self.logger.error(f"Error reading BAM via samtools: {e}")

        # # Open BAM file
        # bamfile = pysam.AlignmentFile(str(bam_file), "rb")
        #
        # # Extract reference lengths from header
        # for ref_name, ref_length in zip(bamfile.references, bamfile.lengths):
        #     gene_lengths[ref_name] = ref_length
        #
        # # Store all reads for later FASTA output
        # all_reads = {}  # {read_name: sequence}
        #
        # # Process alignments
        # try:
        #     for read in bamfile.fetch():
        #         # Store read sequence for later
        #         if read.query_name not in all_reads and read.query_sequence:
        #             all_reads[read.query_name] = read.query_sequence
        #
        #         # Skip unmapped reads
        #         if read.is_unmapped:
        #             continue
        #
        #         gene = read.reference_name
        #         gene_len = gene_lengths.get(gene, 0)
        #
        #         # Track this read maps to this gene (before filtering)
        #         gene_reads[gene]['all'].append(read.query_name)
        #
        #         # Get NM tag (edit distance: mismatches + indels)
        #         try:
        #             nm = read.get_tag('NM')
        #         except KeyError:
        #             nm = 0
        #             print("11")
        #             print(read.query_name)
        #
        #         # Parse CIGAR to get actual aligned positions on reference
        #         # and calculate alignment length for identity
        #         ref_pos = read.reference_start  # 0-based
        #         aligned_positions = set()  # Positions on reference that have aligned bases
        #         alignment_length = 0  # Total alignment columns (for identity calc)
        #
        #             # BAM CIGAR operations:
        #             # 0 = M (match or mismatch) - consumes both query and reference
        #             # 1 = I (insertion to reference) - consumes query only
        #             # 2 = D (deletion from reference) - consumes reference only
        #             # 3 = N (skipped region) - consumes reference only
        #             # 4 = S (soft clip) - consumes query only
        #             # 5 = H (hard clip) - consumes neither
        #             # 6 = P (padding) - consumes neither
        #             # 7 = = (sequence match) - consumes both
        #             # 8 = X (sequence mismatch) - consumes both
        #
        #         for op, length in read.cigartuples:
        #             if op in [0, 7, 8]:  # M, =, X - actual aligned bases
        #                 # Add these reference positions to coverage
        #                 # for i in range(length):
        #                 #     aligned_positions.add(ref_pos + i)
        #                 aligned_positions.update(range(ref_pos, ref_pos + length))
        #                 ref_pos += length
        #                 alignment_length += length
        #
        #             elif op == 1:  # I - insertion (gap in reference)
        #                 # Doesn't consume reference, but counts in alignment length
        #                 alignment_length += length
        #
        #             elif op == 2:  # D - deletion from reference
        #                 # These positions on reference are NOT covered (no bases aligned)
        #                 # But count in alignment length for identity calculation
        #                 ref_pos += length
        #                 alignment_length += length
        #
        #             elif op == 3:  # N (spliced/skipped region)
        #                 ref_pos += length
        #
        #             elif op in [ 4, 5]:  # S, H - soft/hard clips
        #                 # Don't consume reference, don't count in alignment
        #                 pass
        #
        #         # Calculate identity: matches = alignment_length - edit_distance
        #         # This matches BLAST/DIAMOND pident calculation
        #         if alignment_length == 0:
        #             self.logger.warning(
        #                 f"Read {read.query_name} has zero alignment length on gene {gene}. Skipping identity calculation.")
        #             identity = 0
        #         else:
        #             matches = alignment_length - nm
        #             identity = (matches / alignment_length) * 100
        #
        #         if not read.has_tag('NM'):
        #             self.logger.warning(f"Read {read.query_name} missing NM tag for gene {gene}. Assuming NM=0.")
        #
        #         if not read.cigartuples:
        #             self.logger.warning(f"Read {read.query_name} has empty or unusual CIGAR string for gene {gene}.")
        #
        #         # Calculate query coverage
        #         query_length = read.query_length
        #         query_coverage = (len(aligned_positions) / query_length) * 100 if query_length > 0 else 0
        #
        #         # Only process sequences meeting identity threshold
        #         if identity >= self.detection_min_identity and query_coverage >= self.query_min_coverage:
        #             # Initialise stats if first hit for this gene
        #             if gene not in self.gene_stats[database][tool_name]:
        #                 self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)
        #             # Add aligned positions using add_positions method (this handles everything)
        #             self.gene_stats[database][tool_name][gene].add_positions(
        #                 aligned_positions, identity, gene_len
        #             )
        #             # Track reads that pass thresholds
        #             gene_reads[gene]['passing'].append(read.query_name)
        #         else:
        #             if getattr(self, "verbose", True):
        #                 self.logger.info(f"Not passing {gene} {identity} {alignment_length} {read.query_name}")
        #
        # except Exception as e:
        #     self.logger.error(f"Error reading BAM file: {e}")
        # finally:
        #     bamfile.close()

        # Finalise statistics and determine detection based on gene coverage
        for gene in self.gene_stats[database][tool_name]:
            stats = self.gene_stats[database][tool_name][gene]
            stats.finalise()

            # Gene is detected if gene meets thresholds
            if (stats.gene_coverage >= self.detection_min_coverage and
                    stats.base_depth_hit >= self.detection_min_base_depth and
                    stats.num_sequences >= self.detection_min_num_reads):
                detected_genes.add(gene)
                self.detections[database][gene][tool_name] = True

        #self.logger.info(f"Total reads processed in {database.upper()} using {tool_name}: {len(all_reads)}")
        ## Need a fix for this
        self.logger.info(f"Reads that returned a hit in {database.upper()} using {tool_name}: {mapped_reads}")
        self.logger.info(f"Reads passing thresholds in {database.upper()} using {tool_name}: {passing_reads}")
        self.logger.info(f"Detected {len(detected_genes)} genes in {database.upper()} using {tool_name}")

        # Output FASTA files of reads mapping to genes
        if self.report_fasta and detected_genes:
            self._write_fasta_outputs(database, tool_name, detected_genes, gene_reads, all_reads)


        return detected_genes, gene_reads


    def parse_hmmer_results(self, tbl_file: Path, database: str, tool_name: str) -> Set[str]:
        # Parse HMMER table output and extract genes meeting thresholds.
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

        return detected_genes #, gene_reads

    def write_tool_stats(self, database: str, tool_name: str, gene_reads: dict = None):
        """Write detailed statistics for a specific tool to TSV.

        Output columns:
        - Gene: AMR gene name
        - Gene_Length: Length of the gene in the database (bp)
        - Num_Sequences_Mapped: Number of sequences that mapped to this gene with identity >= detection-min-identity
        - Num_Sequences_Passing_Thresholds: Number of sequences that passed all thresholds (identity, coverage, base_coverage, min_reads etc)
        - Gene_Coverage: Percentage of the gene covered by all qualifying alignments combined (%)
        - Base_Coverage: Average base coverage across the entire gene (%)
        - Base_Coverage_Hit: Average base coverage considering only bases with at least one hit (%)
        - Avg_Identity: Average identity across all qualifying sequences (%)
        - Detected: 1 if gene passes all thresholds, 0 otherwise
        """
        stats_file = self.stats_dir / f"{database}_{tool_name}_stats.tsv"

        gene_stats = self.gene_stats[database][tool_name]
        if not gene_stats:
            self.logger.warning(f"No statistics to write for {database} - {tool_name}")
            return

        with open(stats_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            # Header
            header = ['Gene', 'Gene_Length', 'Num_Sequences_Mapped',
                      'Num_Sequences_Passing_Thresholds', 'Gene_Coverage',
                      'Base_Coverage', 'Base_Coverage_Hit', 'Avg_Identity', 'Detected']
            writer.writerow(header)

            # Sort genes alphabetically
            genes = sorted(gene_stats.keys())

            for gene in genes:
                stats = gene_stats[gene]
                try:
                    detected = self.detections[database][gene][tool_name]
                except (KeyError, TypeError):
                    detected = False
                row = [
                    gene,
                    stats.gene_length,
                    len(gene_reads.get(gene, {}).get('all', [])), # 'all' reads mapping to gene
                    len(gene_reads.get(gene, {}).get('passing', [])), # Just those that 'passed' thresholds
                    f"{stats.gene_coverage:.2f}",
                    f"{stats.base_depth:.2f}",
                    f"{stats.base_depth_hit:.2f}",
                    f"{stats.avg_identity:.2f}",
                    '1' if detected else '0'
                ]
                writer.writerow(row)

        self.logger.info(f"  Stats file: {stats_file}")

    def generate_detection_matrix(self, database: str):
        # Generate TSV matrix of gene detections across tools.
        output_file = self.output_dir / f"{database}_detection_matrix.tsv"

        # Get all tools that were run for this database
        all_tools = set()
        for gene_detections in self.detections[database].values():
            all_tools.update(gene_detections.keys())

        if not all_tools:
            self.logger.info(f"No detections found for {database} - No matrix generated.")
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

                genes = [
                    gene for gene in sorted(self.detections[database].keys(), key=get_last_segment)
                    if any(self.detections[database][gene][tool] for tool in all_tools)
                ]
            else:
                genes = [
                    gene for gene in sorted(self.detections[database].keys())
                    if any(self.detections[database][gene][tool] for tool in all_tools)
                ]

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


    def run_workflow(self,options):

        """Run all configured tools on both databases."""
        self.logger.info("=" * 70)
        self.logger.info("AMRfíor - The AMR Gene Detection toolkit: " + AMRFIOR_VERSION)
        self.logger.info("=" * 70)
        ###
        # Log input files (handle new FASTA/FASTQ possibilities)
        if getattr(self, "input_fasta", None):
            self.logger.info(f"Input FASTA: {self.input_fasta}")
        else:
            self.logger.info("Input FASTA: None")

        if getattr(self, "input_fastq", None) is None:
            self.logger.info("Input FASTQ: None")
        else:
            if getattr(self, "input_fastq_is_paired", False):
                self.logger.info(f"Input FASTQ (paired): {self.input_fastq[0]}, {self.input_fastq[1]}")
            else:
                self.logger.info(f"Input FASTQ (single): {self.input_fastq}")
        ###
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Threads: {self.threads}")
        self.logger.info(f"Database(s) chosen:")
        self.logger.info(f"  ResFinder: {'Yes' if 'resfinder'  or 'all' in self.databases else 'No'}")
        self.logger.info(f"  CARD: {'Yes' if 'card' in self.databases or 'all' in self.databases else 'No'}")
        self.logger.info(f"  NCBI: {'Yes' if 'ncbi' in self.databases  or 'all' in self.databases else 'No'}")
        self.logger.info(f"  User-Provided: {'Yes' if 'user_provided' in self.databases else 'No'}")
        self.logger.info(f" Tool(s) chosen:")
        self.logger.info(f"  BLASTn: {'Yes' if 'blastn' in options.tools else 'No'}")
        self.logger.info(f"  BLASTx: {'Yes' if 'blastx' in options.tools else 'No'}")
        self.logger.info(f"  DIAMOND: {'Yes' if 'diamond' in options.tools else 'No'}")
        self.logger.info(f"  Bowtie2: {'Yes' if 'bowtie2' in options.tools else 'No'}")
        self.logger.info(f"  BWA: {'Yes' if 'bwa' in options.tools else 'No'}")
        self.logger.info(f"  Minimap2: {'Yes' if 'minimap2' in options.tools else 'No'}")

       # self.logger.info(f"E-value threshold: {self.evalue}")
        self.logger.info(f"Min query coverage: {self.query_min_coverage}%")
        self.logger.info(f"Min detection coverage: {self.detection_min_coverage}%")
        self.logger.info(f"Min detection identity: {self.detection_min_identity}%")
        self.logger.info(f"Run DNA mode: {self.run_dna}")
        self.logger.info(f"Run Protein mode: {self.run_protein}")
        params_str = ", ".join(
            f"{tool}: {params}" for tool, params in self.tool_sensitivity_params.items()
        ) if self.tool_sensitivity_params else "None"
        self.logger.info(f"Sensitivity parameters: {options.sensitivity} - {params_str}")
        self.logger.info("=" * 70)

        results = {}

        # Iterate over each database in the provided `databases` dictionary
        for db_name, db_paths in self.databases.items():
            results[db_name] = {}
            self.logger.info(f"\n### Processing {db_name.capitalize()} Database ###")

            if self.run_dna and db_paths.get('blastn'):
                results[db_name]['BLASTn-DNA'] = self.run_blast(
                    db_paths['blastn'], db_name, 'dna')

            if self.run_protein and db_paths.get('blastx'):
                results[db_name]['BLASTx-AA'] = self.run_blast(
                    db_paths['blastx'], db_name, 'protein')

            if self.run_protein and db_paths.get('diamond'):
                results[db_name]['DIAMOND-AA'] = self.run_diamond(
                    db_paths['diamond'], db_name)

            if self.run_dna and db_paths.get('bowtie2'):
                results[db_name]['Bowtie2-DNA'] = self.run_bowtie2(
                    db_paths['bowtie2'], db_name)

            if self.run_dna and db_paths.get('bwa'):
                results[db_name]['BWA-DNA'] = self.run_bwa(
                    db_paths['bwa'], db_name)

            if db_paths.get('minimap2'):
                results[db_name]['Minimap2-DNA'] = self.run_minimap2(
                    db_paths['minimap2'], db_name, options.minimap2_preset)

            # Uncomment the following lines if HMMER support is added
            # if self.run_dna and db_paths.get('hmmer_dna'):
            #     results[db_name]['HMMER-DNA'] = self.run_hmmer(
            #         db_paths['hmmer_dna'], db_name, 'dna')
            #
            # if self.run_protein and db_paths.get('hmmer_protein'):
            #     results[db_name]['HMMER-PROTEIN'] = self.run_hmmer(
            #         db_paths['hmmer_protein'], db_name, 'protein')

            self.generate_detection_matrix(db_name)

            if options.report_fasta != None:
                # Write gene name changes to TSV if any changes were made
                if hasattr(self, 'gene_name_changes') and self.gene_name_changes:
                    changes_file = self.fasta_dir / f"{db_name}_gene_name_changes.tsv"
                    with open(changes_file, "w", newline='') as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerow(['Original_Gene_Name', 'Safe_Gene_Name'])
                        writer.writerows(self.gene_name_changes)
                    self.logger.info(f"  Gene name changes file: {changes_file}")
                    self.gene_name_changes.clear()

        # results = {'resfinder': {}, 'card': {}}
        #
        # # Process ResFinder database
        # if self.resfinder_dbs:
        #     self.logger.info("\n### Processing ResFinder Database ###")
        #
        #     if self.run_dna and self.resfinder_dbs.get('blastn'):
        #         results['resfinder']['BLASTn-DNA'] = self.run_blast(
        #             self.resfinder_dbs['blastn'], 'resfinder', 'dna')
        #
        #     if self.run_protein and self.resfinder_dbs.get('blastx'):
        #         results['resfinder']['BLASTx-AA'] = self.run_blast(
        #             self.resfinder_dbs['blastx'], 'resfinder', 'protein')
        #
        #     if self.run_protein and self.resfinder_dbs.get('diamond'):
        #         results['resfinder']['DIAMOND-AA'] = self.run_diamond(
        #             self.resfinder_dbs['diamond'], 'resfinder')
        #
        #     if self.run_dna and self.resfinder_dbs.get('bowtie2'):
        #         results['resfinder']['Bowtie2-DNA'] = self.run_bowtie2(
        #             self.resfinder_dbs['bowtie2'], 'resfinder')
        #
        #     if self.run_dna and self.resfinder_dbs.get('bwa'):
        #         results['resfinder']['BWA-DNA'] = self.run_bwa(
        #             self.resfinder_dbs['bwa'], 'resfinder')
        #
        #     if self.resfinder_dbs.get('minimap2'):
        #         results['resfinder']['Minimap2-DNA'] = self.run_minimap2(
        #             self.resfinder_dbs['minimap2'], 'resfinder', options.minimap2_preset)
        #
        #     # if self.run_dna and self.resfinder_dbs.get('hmmer_dna'):
        #     #     results['resfinder']['HMMER-DNA'] = self.run_hmmer(
        #     #         self.resfinder_dbs['hmmer_dna'], 'resfinder', 'dna')
        #     #
        #     # if self.run_protein and self.resfinder_dbs.get('hmmer_protein'):
        #     #     results['resfinder']['HMMER-PROTEIN'] = self.run_hmmer(
        #     #         self.resfinder_dbs['hmmer_protein'], 'resfinder', 'protein')
        #
        #     self.generate_detection_matrix('resfinder')
        #
        #     if options.report_fasta != None:
        #         # Write gene name changes to TSV if any changes were made
        #         if hasattr(self, 'gene_name_changes') and self.gene_name_changes:
        #             changes_file = self.fasta_dir / f"resfinder_gene_name_changes.tsv"
        #             with open(changes_file, "w", newline='') as f:
        #                 writer = csv.writer(f, delimiter='\t')
        #                 writer.writerow(['Original_Gene_Name', 'Safe_Gene_Name'])
        #                 writer.writerows(self.gene_name_changes)
        #             self.logger.info(f"  Gene name changes file: {changes_file}")
        #             self.gene_name_changes.clear()
        #
        # # Process CARD database
        # if self.card_dbs:
        #     self.logger.info("\n### Processing CARD Database ###")
        #
        #     if self.run_dna and self.card_dbs.get('blastn'):
        #         results['card']['BLASTn-DNA'] = self.run_blast(
        #             self.card_dbs['blastn'], 'card', 'dna')
        #
        #     if self.run_protein and self.card_dbs.get('blastx'):
        #         results['card']['BLASTx-AA'] = self.run_blast(
        #             self.card_dbs['blastx'], 'card', 'protein')
        #
        #     if self.run_protein and self.card_dbs.get('diamond'):
        #         results['card']['DIAMOND-AA'] = self.run_diamond(
        #             self.card_dbs['diamond'], 'card')
        #
        #     if self.run_dna and self.card_dbs.get('bowtie2'):
        #         results['card']['Bowtie2-DNA'] = self.run_bowtie2(
        #             self.card_dbs['bowtie2'], 'card')
        #
        #     if self.run_dna and self.card_dbs.get('bwa'):
        #         results['card']['BWA-DNA'] = self.run_bwa(
        #             self.card_dbs['bwa'], 'card')
        #
        #     if self.card_dbs.get('minimap2'):
        #         results['card']['Minimap2-DNA'] = self.run_minimap2(
        #             self.card_dbs['minimap2'], 'card', options.minimap2_preset)
        #
        #     # if self.run_dna and self.card_dbs.get('hmmer_dna'):
        #     #     results['card']['HMMER-DNA'] = self.run_hmmer(
        #     #         self.card_dbs['hmmer_dna'], 'card', 'dna')
        #     #
        #     # if self.run_protein and self.card_dbs.get('hmmer_protein'):
        #     #     results['card']['HMMER-PROTEIN'] = self.run_hmmer(
        #     #         self.card_dbs['hmmer_protein'], 'card', 'protein')
        #
        #     self.generate_detection_matrix('card')
        #     if options.report_fasta != None:
        #         # Write gene name changes to TSV if any changes were made
        #         if hasattr(self, 'gene_name_changes') and self.gene_name_changes:
        #             changes_file = self.fasta_dir / f"card_gene_name_changes.tsv"
        #             with open(changes_file, "w", newline='') as f:
        #                 writer = csv.writer(f, delimiter='\t')
        #                 writer.writerow(['Original_Gene_Name', 'Safe_Gene_Name'])
        #                 writer.writerows(self.gene_name_changes)
        #             self.logger.info(f"  Gene name changes file: {changes_file}")
        #             self.gene_name_changes.clear()

        # Final summary
        self.logger.info("\n" + "=" * 70)
        self.logger.info("PIPELINE SUMMARY")
        self.logger.info("=" * 70)

        for db_name in self.databases.keys():
            if results[db_name]:
                self.logger.info(f"\n{db_name.upper()}:")
                for tool, val in results[db_name].items():
                    if isinstance(val, tuple) and len(val) == 2:
                        success, genes = val
                    elif isinstance(val, bool):
                        success, genes = val, set()
                    elif isinstance(val, set):
                        success, genes = True, val
                    else:
                        success, genes = bool(val), set()

                    status = "✓" if success else "✗"
                    gene_count = len(genes) if isinstance(genes, (set, list, tuple, dict)) else 0
                    self.logger.info(f"  {status} {tool:.<30} {gene_count} genes detected")
        # for db_name in self.databases.keys():
        #     if results[db_name]:
        #         self.logger.info(f"\n{db_name.upper()}:")
        #         for tool, (success, genes) in results[db_name].items():
        #             status = "✓" if success else "✗"
        #             gene_count = len(genes) if success else 0
        #             self.logger.info(f"  {status} {tool:.<30} {gene_count} genes detected")

        self.logger.info("=" * 70)
        self.logger.info(f"Detection matrices saved to: {self.output_dir}")
        self.logger.info(f"Tool statistics saved to: {self.stats_dir}")
        self.logger.info(f"Raw outputs saved to: {self.raw_dir}")
        if self.report_fasta != None:
            self.logger.info(f"FASTA reports saved to: {self.fasta_dir}")
        self.logger.info("=" * 70)

        return results