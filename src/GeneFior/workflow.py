import subprocess
import sys
import shutil
import threading
import concurrent.futures
import csv
import os
import gzip
import tempfile
import time
import signal
from collections import defaultdict
from pathlib import Path
import logging
import shlex
from typing import Any, Dict, List, Set, Tuple, Iterator
import re
from typing import Optional


_WHOLE_GENOME_CONTIG_AFFIXES = (
    'contig', 'ctg', 'chromosome', 'chrom', 'chr', 'scaffold',
    'scaf', 'replicon', 'segment', 'plasmid',
)
_WHOLE_GENOME_CONTIG_SUFFIX_RE = re.compile(
    r'(?P<genome>.+?)(?:[_-](?P<contig>(?:'
    + '|'.join(_WHOLE_GENOME_CONTIG_AFFIXES)
    + r')(?:[_-]?[A-Za-z0-9.-]+)?))$',
    re.IGNORECASE,
)


def split_whole_genome_reference(reference: str) -> Tuple[str, str]:
    """Split a whole-genome reference into genome and contig identifiers.

    Explicit separators are preferred. Underscore and hyphen splitting is
    only enabled when the suffix is a recognised contig affix, which avoids
    merging arbitrary references that merely share a naming prefix.
    Unrecognised names remain single-contig genomes.
    """
    value = str(reference).strip()
    for separator in ('::', '|', '#', '/'):
        if separator in value:
            genome, contig = value.split(separator, 1)
            if genome and contig:
                return genome, contig
    match = _WHOLE_GENOME_CONTIG_SUFFIX_RE.match(value)
    if match:
        return match.group('genome'), match.group('contig')
    return value, value


try:
    from .gene_stats import GeneStats
    from .evidence import (
        EvidenceConfig,
        DETECTION_SYSTEM_LEGACY_RELAXED,
        DETECTION_SYSTEM_QUALIFIED,
        EXACT_PROTEIN_DETECTED,
        NOT_DETECTED,
        WHOLE_GENOME_MAPPED,
        WHOLE_GENOME_PARTIAL,
        WHOLE_GENOME_NEAR_COMPLETE,
        WHOLE_GENOME_COMPLETE,
        MIXED_OR_MOSAIC,
        PARTIAL_OR_DIVERGENT,
        best_status,
        classify_gene_evidence,
        classify_gene_legacy,
        normalise_detection_system,
        resolve_family_calls,
        scores_tied,
        tool_modality,
    )
    from .constants import *
except (ModuleNotFoundError, ImportError) as error:
    from gene_stats import GeneStats
    from evidence import (
        EvidenceConfig,
        DETECTION_SYSTEM_LEGACY_RELAXED,
        DETECTION_SYSTEM_QUALIFIED,
        EXACT_PROTEIN_DETECTED,
        NOT_DETECTED,
        WHOLE_GENOME_MAPPED,
        WHOLE_GENOME_PARTIAL,
        WHOLE_GENOME_NEAR_COMPLETE,
        WHOLE_GENOME_COMPLETE,
        MIXED_OR_MOSAIC,
        PARTIAL_OR_DIVERGENT,
        best_status,
        classify_gene_evidence,
        classify_gene_legacy,
        normalise_detection_system,
        resolve_family_calls,
        scores_tied,
        tool_modality,
    )
    from constants import *

class Workflow:
    # Orchestrates alignment tools and result parsing for gene detection.

    def __init__(self, input_fasta: str, input_fastq: str, output_dir: str,
                 databases: Dict[str, Dict[str, str]],
                 threads: int = 4,  #max_target_seqs: int = 100,
                 tool_sensitivity_params: Dict[str, Dict[str, Any]] = None,
                 blastx_task: str = 'blastx-fast',
                 #evalue: float = 1e-10,
                 query_min_coverage: float = 40.0,  #Query coverage threshold
                 query_min_identity: float = 80.0,  # Query identity threshold
                 detection_min_coverage: float = 80.0,
                 detection_min_identity: float = 80.0,
                 detection_min_depth: int = 3,
                 detection_min_num_reads: int = 1,
                  evalue: Optional[float] = None,
                  min_bitscore: Optional[float] = None,

                 run_dna: bool = True, run_protein: bool = True,
                 sequence_type: str = 'Single-FASTA',
                 report_fasta: str = None,
                 no_cleanup: bool = False,
                 verbose: bool = False,
                 logger: Optional[logging.Logger] = None,
                 max_fasta_chunk_bytes: int = 200 * 1024 * 1024,
                 temp_directory: Optional[str] = None,
                 chunk_jobs: Optional[int] = None,
                 chunk_threads_per_job: Optional[int] = None,
                   preserve_chunks: bool = False,
                     extra_tool_params: Dict[str, str] = None,
                     genes_type: Optional[str] = None,
                     hmmer_evalue: Optional[float] = None,
                     hmmer_threshold_mode: str = 'evalue',
                     hmmer_must_flag: bool = True,
                     always_flag_genes: Optional[Set[str]] = None,
                     tools: Optional[List[str]] = None,
                     evidence_corroborating_depth: int = 3,
                     evidence_exact_identity: float = 100.0,
                     evidence_candidate_depth: int = 3,
                     evidence_candidate_identity: float = 98.0,
                     evidence_max_internal_gap_bp: int = 15,
                     evidence_max_internal_gap_fraction: float = 0.02,
                     evidence_min_unique_reads: int = 2,
                     evidence_min_unique_fraction: float = 0.10,
                     evidence_ambiguity_fraction: float = 0.50,
                     evidence_score_tie: float = 1.0,
                     detection_system: str = DETECTION_SYSTEM_QUALIFIED,
                     db_whole_genome: bool = False):

                 
        ### Handle input FASTA and FASTQ
        # Always define the attribute `input_fasta` (may be None) to avoid AttributeError
        # when other methods reference it. Use Path when a value is provided.
        self.input_fasta = Path(input_fasta) if input_fasta is not None else None
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
        self.blastx_task = blastx_task
       # self.evalue = evalue
        self.query_min_coverage = query_min_coverage
        self.query_min_identity = query_min_identity
        self.detection_min_coverage = detection_min_coverage
        self.detection_min_identity = detection_min_identity
        self.detection_min_depth = max(3, int(detection_min_depth or 3))
        self.detection_min_num_reads = detection_min_num_reads
        self.detection_system = normalise_detection_system(detection_system)
        self.evidence_config = EvidenceConfig(
            corroborating_depth=evidence_corroborating_depth,
            exact_identity_min=evidence_exact_identity,
            candidate_depth=evidence_candidate_depth,
            candidate_identity_min=evidence_candidate_identity,
            max_internal_gap_bp=evidence_max_internal_gap_bp,
            max_internal_gap_fraction=evidence_max_internal_gap_fraction,
            min_unique_best_reads=evidence_min_unique_reads,
            min_unique_best_fraction=evidence_min_unique_fraction,
            ambiguity_warning_fraction=evidence_ambiguity_fraction,
            score_tie_absolute=evidence_score_tie,
        )

        # Threshold (bytes) to split large FASTA into chunks for parallel BLAST.
        self.max_fasta_chunk_bytes = max_fasta_chunk_bytes

        # Use the user-provided temp directory for chunk files if given, otherwise fall back to output_dir
        if temp_directory:
            self.temp_directory = Path(temp_directory)
        else:
            self.temp_directory = self.output_dir

        # Parallel chunking controls
        # Number of concurrent chunk BLAST jobs to run. If not provided default to 1 (serial)
        # Use explicit None check to avoid treating 0 or other falsy-but-intentional values incorrectly.
        if chunk_jobs is None:
            self.chunk_jobs = 1
        else:
            try:
                self.chunk_jobs = int(chunk_jobs)
            except Exception:
                self.chunk_jobs = 1
        # Track whether the user explicitly provided a chunk_jobs value
        self.chunk_jobs_user_specified = (chunk_jobs is not None)

        # Threads to allocate per chunk job. If None we compute an even split of self.threads
        if chunk_threads_per_job is None:
            self.chunk_threads_per_job = None
        else:
            try:
                self.chunk_threads_per_job = int(chunk_threads_per_job)
            except Exception:
                self.chunk_threads_per_job = None
        # Preserve chunks flag (do not delete part files) - default False
        self.preserve_chunks = bool(preserve_chunks)
        # E-value threshold for tools that accept it. If None, do not pass an
        # explicit e-value to the tool and let the tool use its own default.
        self.evalue = float(evalue) if evalue is not None else None
        # Minimum bitscore threshold for BLAST/DIAMOND hits. If None, no bitscore
        # filtering is applied and tool defaults / post-filtering rules are used.
        self.min_bitscore = float(min_bitscore) if min_bitscore is not None else None
        # Arbitrary extra parameters to append to specific tools (dict: tool -> param string)
        self.extra_tool_params = extra_tool_params or {}
        # HMMER-specific E-value (overrides global evalue for HMMER only).
        # If None, falls back to self.evalue; if both are None HMMER uses its own default (E=10).
        self.hmmer_evalue = float(hmmer_evalue) if hmmer_evalue is not None else None
        # Threshold mode for HMMER: 'evalue' | 'tc' | 'ga' | 'nc'
        # When 'tc'/'ga'/'nc', per-profile cutoffs are passed to hmmsearch/nhmmer
        # and E-value post-filtering is disabled.
        self.hmmer_threshold_mode = hmmer_threshold_mode or 'evalue'
        # The input type is indicated via self.sequence_type. If it equals
        # 'Genes-FASTA' the pipeline treats inputs as full-length gene FASTA(s)
        # (and will skip read-mappers). The optional `genes_type` parameter,
        # when provided, forces interpretation to either 'dna' or 'protein'
        # and bypasses auto-detection.
        self.genes_type = genes_type
        # Track whether this workflow was invoked with Genes-FASTA input. Use an
        # explicit boolean so we don't mutate `self.sequence_type` (which is
        # used elsewhere to decide flags). This preserves the original signal
        # while allowing downstream logic to check `is_genes_fasta`.
        self.is_genes_fasta = (sequence_type == 'Genes-FASTA')
        # Store any detected/declared genes alphabet type ('dna'|'protein'|None)
        self.detected_genes_type = None
        # Whether the provided database(s) are whole-genome references. When True,
        # downstream logic should prefer read-mapping tools and compute per-base
        # coverage summaries instead of gene-centric detection metrics.
        self.db_whole_genome = bool(db_whole_genome)


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
        if self.report_fasta != None:
            self.fasta_dir = self.output_dir / "fasta_outputs"
            self.fasta_dir.mkdir(exist_ok=True)

        # misc
        self.no_cleanup = no_cleanup
        self.verbose = verbose
        # Ensure we always have a logger object to call; fall back to module logger if none provided
        if logger is None:
            self.logger = logging.getLogger('GeneFior')
        else:
            self.logger = logger

        # If logger has no handlers, configure a simple StreamHandler at INFO so runtime tool messages are visible.
        try:
            if not getattr(self.logger, 'handlers', None):
                handler = logging.StreamHandler()
                handler.setLevel(logging.INFO)
                formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
                # If verbose requested, set DEBUG level so internal debug messages are visible
                if getattr(self, 'verbose', False):
                    self.logger.setLevel(logging.DEBUG)
                    handler.setLevel(logging.DEBUG)
                else:
                    self.logger.setLevel(logging.INFO)
                # Avoid double logging to root handlers
                self.logger.propagate = False
        except Exception:
            # If anything goes wrong configuring logger, fall back silently
            pass


        # Store detection results: {database: {gene: {tool: bool}}}
        self.detections = {
            db_name: defaultdict(lambda: defaultdict(bool))
            for db_name in self.databases.keys()
        }
        # Qualified evidence calls parallel the strict legacy boolean matrix.
        # detections=True means threshold-passing evidence; the status records
        # whether it is an allele, protein, profile, or whole-genome mapping call.
        self.evidence_calls = {
            db_name: defaultdict(dict)
            for db_name in self.databases.keys()
        }
        self.family_calls = {
            db_name: defaultdict(dict)
            for db_name in self.databases.keys()
        }

        # Store detailed statistics: {database: {tool: {gene: GeneStats}}}
        self.gene_stats = {
            db_name: defaultdict(lambda: defaultdict(GeneStats))
            for db_name in self.databases.keys()
        }

        # Cache for HMMER annotations: {db_name: {gene_id: {'description': str, 'must_flag': bool}}}
        # Populated lazily when run_hmmer is called and a 'hmmer_annotations' path is present in db_paths.
        self.hmmer_annotations: Dict[str, Dict[str, dict]] = {}

        # Metadata about each loaded annotation file, keyed by db_name.
        # Currently tracks: has_must_flag (bool) — whether the source CSV contained a
        # 'Must flag' column at all.  Used to suppress irrelevant columns/alerts for
        # databases that have no concept of must-flagging.
        self.hmmer_annotations_meta: Dict[str, dict] = {}

        # When False the must-flag override is completely disabled: must-flag genes that
        # fail coverage/min-reads gates are treated identically to all other genes and
        # the MUST-FLAG ALERT log block is suppressed.  True by default.
        self.hmmer_must_flag = hmmer_must_flag
        self.always_flag_genes = set(always_flag_genes or [])
        # User-selected tools list (e.g. ['bowtie2', 'bwa']).
        # When provided, only the listed tools are run regardless of what databases
        # are available.  When None (not provided by caller) all available tools run
        # (legacy behaviour preserved for backwards compatibility).
        self.tools = list(tools) if tools is not None else None

    def _tool_enabled(self, tool_name: str) -> bool:
        """Return True if `tool_name` is in the user-selected tools list.
        If no tools list was provided (self.tools is None) every tool is allowed."""
        if self.tools is None:
            return True
        return tool_name in self.tools

    def _classify_tool_evidence(self, database: str, tool_name: str,
                                must_flag_overrides: Set[str] = None,
                                eligible_genes: Set[str] = None) -> Set[str]:
        """Classify all genes for one tool under the selected detection system."""
        must_flag_overrides = must_flag_overrides or set()
        self.evidence_calls.setdefault(database, defaultdict(dict))
        self.family_calls.setdefault(database, defaultdict(dict))
        self.detections.setdefault(
            database,
            defaultdict(lambda: defaultdict(bool)),
        )

        stats_by_gene = self.gene_stats.get(database, {}).get(tool_name, {})
        calls = {}
        detection_system = normalise_detection_system(
            getattr(self, "detection_system", DETECTION_SYSTEM_QUALIFIED)
        )
        for gene, stats in stats_by_gene.items():
            always_flagged = gene in getattr(self, 'always_flag_genes', set())
            if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED:
                call = classify_gene_legacy(
                    gene=gene,
                    tool_name=tool_name,
                    stats=stats,
                    detection_min_coverage=self.detection_min_coverage,
                    detection_min_identity=self.detection_min_identity,
                    detection_min_depth=self.detection_min_depth,
                    detection_min_num_reads=self.detection_min_num_reads,
                    reads_mode=not getattr(self, 'is_genes_fasta', False),
                )
            else:
                call = classify_gene_evidence(
                    gene=gene,
                    tool_name=tool_name,
                    stats=stats,
                    config=self.evidence_config,
                    detection_min_coverage=self.detection_min_coverage,
                    detection_min_identity=self.detection_min_identity,
                    detection_min_num_reads=self.detection_min_num_reads,
                    reads_mode=not getattr(self, 'is_genes_fasta', False),
                    must_flag_override=(
                        gene in must_flag_overrides or always_flagged
                    ),
                    whole_genome=getattr(self, 'db_whole_genome', False),
                )
            if always_flagged:
                call.warnings = sorted(
                    set(call.warnings + ["ALWAYS_FLAGGED_GENE"])
                )
            if eligible_genes is not None and gene not in eligible_genes:
                call.status = NOT_DETECTED
                call.exact_detected = False
                call.evidence_present = False
                call.warnings = sorted(set(call.warnings + ["SIGNIFICANCE_FILTER_FAILED"]))
            calls[gene] = call

        family_summaries = (
            {}
            if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED
            else resolve_family_calls(
                calls,
                stats_by_gene,
                self.evidence_config,
            )
        )
        self.family_calls[database][tool_name] = family_summaries

        evidence = set()
        status_counts = defaultdict(int)
        review_examples = []
        for gene, call in calls.items():
            self.evidence_calls[database][gene][tool_name] = call
            self.detections[database][gene][tool_name] = call.evidence_present
            status_counts[call.status] += 1
            if call.evidence_present:
                evidence.add(gene)
                if (call.status == "MIXED_OR_MOSAIC"
                        and len(review_examples) < 10):
                    review_examples.append(f"{gene}={call.status}")
        summary = ", ".join(
            f"{status}={count}"
            for status, count in sorted(status_counts.items())
            if count
        )
        if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED:
            self.logger.info(
                f"Detection classification for {database}/{tool_name} "
                f"[{detection_system}]: DETECTED={len(evidence)}"
            )
        else:
            self.logger.info(
                f"Detection classification for {database}/{tool_name} "
                f"[{detection_system}]: {summary}"
            )
        if review_examples:
            self.logger.warning(
                "Priority evidence reviews (first 10): "
                + ", ".join(review_examples)
            )
        return evidence

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
            self.logger.info(f"Input FASTA ` {fasta_path_str} ` appears to be a valid gzip file.")
            return True

    def _terminate_proc(self, proc, name: str = None, timeout: float = 5.0):
        """Attempt to cleanly terminate a subprocess.Popen `proc` and its children.
        Best-effort: try terminate, wait, then kill. Uses process groups when available.
        """
        if proc is None:
            return
        try:
            # If process already finished, nothing to do
            if proc.poll() is not None:
                return
        except Exception:
            pass

        try:
            pid = getattr(proc, 'pid', None)
            if pid is None:
                try:
                    proc.terminate()
                except Exception:
                    pass
                return

            # Try to terminate the process group if the child was started in a new session
            try:
                os.killpg(pid, signal.SIGTERM)
            except Exception:
                # Fall back to terminating the single proc
                try:
                    proc.terminate()
                except Exception:
                    pass

            # Wait up to timeout for process to exit
            end = time.time() + float(timeout)
            while time.time() < end:
                if proc.poll() is not None:
                    break
                time.sleep(0.1)

            if proc.poll() is None:
                # Force kill
                try:
                    os.killpg(pid, signal.SIGKILL)
                except Exception:
                    try:
                        proc.kill()
                    except Exception:
                        pass
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass

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

        gzip_proc = None
        proc = None
        stderr_thread = None
        from collections import deque
        stdout_buffer = deque(maxlen=200)

        try:
            self.logger.info(f"Piping decompressed ` {self.input_fasta} ` to BLAST ({tool_name})")
            # Start child processes in their own session so we can kill process groups if needed
            gzip_proc = subprocess.Popen(['gzip', '-dc', fasta_path_str], stdout=subprocess.PIPE, start_new_session=True, close_fds=True)
            # Run BLAST reading from stdin using Popen so we can stream stderr/stdout live
            proc = subprocess.Popen(cmd, stdin=gzip_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, text=True, start_new_session=True, close_fds=True)

            # Stream stderr in background
            def _stream_stderr_local(stderr_pipe, tname):
                try:
                    if stderr_pipe is None:
                        return
                    buf = ''
                    while True:
                        chunk = stderr_pipe.read(1024)
                        if chunk == '':
                            break
                        if not chunk:
                            continue
                        buf += chunk
                        parts = re.split(r"[\r\n]", buf)
                        if buf and buf[-1] not in ('\n', '\r'):
                            buf = parts.pop()
                        else:
                            buf = ''
                        for part in parts:
                            s = part.strip()
                            if s:
                                try:
                                    self.logger.info(f"{tname} stderr: {s}")
                                except Exception:
                                    pass
                    if buf:
                        s = buf.strip()
                        if s:
                            try:
                                self.logger.info(f"{tname} stderr: {s}")
                            except Exception:
                                pass
                except Exception:
                    pass

            if proc and proc.stderr:
                stderr_thread = threading.Thread(target=_stream_stderr_local, args=(proc.stderr, tool_name), daemon=True)
                stderr_thread.start()

            # Stream stdout line-by-line to avoid buffering large quantities in memory
            if proc and proc.stdout:
                for out_line in proc.stdout:
                    try:
                        if out_line is None:
                            continue
                        s = out_line.rstrip('\n')
                        # BLAST/DIAMOND info lines often start with '#'
                        if s.startswith('#'):
                            try:
                                self.logger.info(f"{tool_name} stdout: {s}")
                            except Exception:
                                pass
                        else:
                            # Store a bounded amount of non-comment stdout for debug logging later
                            stdout_buffer.append(s)
                    except Exception:
                        continue

            if proc:
                proc.wait()

        except subprocess.CalledProcessError as e:
            self.logger.error(f"{tool_name} failed with return code {e.returncode}")
            self.logger.error(f"Error message: {getattr(e, 'stderr', '')}")
            return False
        except FileNotFoundError:
            self.logger.error(f"{tool_name} executable not found. Is it in your PATH?")
            return False
        except Exception as e:
            self.logger.error(f"Error running {tool_name} with piped input: {e}")
            return False
        finally:
            # Ensure subprocesses and related resources are cleaned up
            try:
                if proc:
                    try:
                        if getattr(proc, 'stdout', None):
                            proc.stdout.close()
                    except Exception:
                        pass
                    try:
                        if getattr(proc, 'stderr', None):
                            proc.stderr.close()
                    except Exception:
                        pass
                    # Best-effort terminate
                    try:
                        self._terminate_proc(proc, tool_name)
                    except Exception:
                        pass
            except Exception:
                pass

            try:
                if gzip_proc:
                    try:
                        if getattr(gzip_proc, 'stdout', None):
                            gzip_proc.stdout.close()
                    except Exception:
                        pass
                    try:
                        self._terminate_proc(gzip_proc, 'gzip')
                    except Exception:
                        pass
            except Exception:
                pass

            if stderr_thread and stderr_thread.is_alive():
                try:
                    stderr_thread.join(timeout=1.0)
                except Exception:
                    pass

        # Completed successfully
        self.logger.info(f"{tool_name} completed successfully")
        if stdout_buffer:
            try:
                stdout_text = '\n'.join(list(stdout_buffer))
                self.logger.debug(f"{tool_name} stdout (tail):\n{stdout_text}")
            except Exception:
                pass
        return True

    def _detect_genes_type(self, max_sequences: int = 50) -> Optional[str]:
        """Heuristic to detect whether a provided Genes-FASTA contains DNA or protein sequences.
        Returns 'dna', 'protein', or None if unknown. Reads up to `max_sequences` sequences from the
        input FASTA (handles gzip)."""
        try:
            fasta_path = str(self.input_fasta)
            is_gz = fasta_path.endswith(('.gz', '.gzip'))
            opener = gzip.open if is_gz else open
            dna_chars = set(list('ACGTNacgtn'))
            protein_chars = set(list('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvw y'))
            # Keep simple counters of characters observed
            dna_score = 0
            protein_score = 0
            seq_count = 0
            with opener(fasta_path, 'rt') as fh:
                seq = ''
                for line in fh:
                    if line.startswith('>'):
                        if seq:
                            seq_count += 1
                            # sample characters
                            s = seq.strip()
                            if not s:
                                seq = ''
                                continue
                            # Count letters
                            letters = [c for c in s if c.isalpha()]
                            if not letters:
                                seq = ''
                                continue
                            # Compute simple ratio of DNA-like letters
                            letters_set = set(letters)
                            if letters_set.issubset(dna_chars):
                                dna_score += 1
                            else:
                                # if many letters outside DNA alphabet, treat as protein
                                protein_score += 1
                            seq = ''
                            if seq_count >= max_sequences:
                                break
                    else:
                        seq += line.strip()
                # handle last sequence
                if seq and seq_count < max_sequences:
                    s = seq.strip()
                    letters = [c for c in s if c.isalpha()]
                    if letters:
                        letters_set = set(letters)
                        if letters_set.issubset(dna_chars):
                            dna_score += 1
                        else:
                            protein_score += 1
            # Decide
            if dna_score > 0 and protein_score == 0:
                return 'dna'
            if protein_score > 0 and dna_score == 0:
                return 'protein'
            # Mixed or unknown
            return None
        except Exception:
            return None

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
            # Log any informational output written to stderr by the tool
            if result.stderr:
                try:
                    stderr_text = result.stderr.strip()
                    if stderr_text:
                        self.logger.info(f"{tool_name} stderr:\n{stderr_text}")
                except Exception:
                    pass
            if result.stdout:
                self.logger.debug(f"{tool_name} stdout: {result.stdout}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"{tool_name} failed with return code {e.returncode}")
            self.logger.error(f"Error message: {e.stderr}")
            # If extra parameters were supplied for this tool, warn the user that
            # the failure may be caused by incompatible/invalid parameters.
            try:
                lowered = tool_name.lower()
                for key, params in (self.extra_tool_params or {}).items():
                    if key and key.lower() in lowered:
                        self.logger.warning("\n" + "*" * 80)
                        self.logger.warning(f"WARNING: Additional parameters provided for tool '{key}': {params}")
                        self.logger.warning("These parameters may be incompatible with the tool and could have caused the failure.")
                        self.logger.warning("Please re-run without the extra params or verify they are valid for this tool.")
                        self.logger.warning("" + "*" * 80 + "\n")
                        break
            except Exception:
                pass
            return False
        except FileNotFoundError:
            self.logger.error(f"{tool_name} executable not found. Is it in your PATH?")
            return False

    def _stream_filter_hits(self, cmd: List[str], output_file: Path, tool_name: str, gz_input: bool = False, fasta_path_str: str = None) -> bool:
        """Run a tabular-output tool and filter hits on-the-fly.

        Writes passing lines to output_file. If gz_input is True and fasta_path_str is
        given, streams decompressed FASTA into the tool stdin.
        """
        gzip_proc = None
        proc = None
        try:
            # Normalise command elements to strings to avoid accidental bools or other types
            cmd = [str(x) for x in cmd]

            stdin_src = None
            # Normalise fasta_path_str to either None or a string; avoid bools or other types
            if isinstance(fasta_path_str, bool):
                fasta_path_str = None
            if fasta_path_str is not None:
                try:
                    fasta_path_str = str(fasta_path_str)
                except Exception:
                    fasta_path_str = None

            if gz_input:
                # If fasta_path_str not provided, attempt to derive from self.input_fasta
                if not fasta_path_str and hasattr(self, 'input_fasta'):
                    try:
                        fasta_path_str = str(self.input_fasta)
                    except Exception:
                        fasta_path_str = None

                # Only attempt gzip piping if we have a valid existing file path
                if not fasta_path_str or not os.path.isfile(fasta_path_str):
                    self.logger.warning(f"gz_input requested but fasta path is not a valid file: {repr(fasta_path_str)}; proceeding without gzip piping")
                    gz_input = False

            if gz_input and fasta_path_str:
                # start gzip in its own session so we can terminate its group if needed
                gzip_proc = subprocess.Popen(['gzip', '-dc', fasta_path_str], stdout=subprocess.PIPE, start_new_session=True, close_fds=True)
                stdin_src = gzip_proc.stdout

            # Ensure parent directory exists
            output_file.parent.mkdir(parents=True, exist_ok=True)


            # Use line-buffered I/O to improve live logging responsiveness
            # Try to force line-buffered IO from the child process so we can stream progress messages.
            # If `stdbuf` is available, prepend it to the command: `stdbuf -oL -eL <cmd>`
            try:
                if shutil.which('stdbuf'):
                    cmd = ['stdbuf', '-oL', '-eL'] + cmd
                    try:
                        self.logger.debug(f"Prepended stdbuf for line-buffering: {' '.join(cmd)}")
                    except Exception:
                        pass
                elif shutil.which('unbuffer'):
                    # `unbuffer` from expect/tcl can also be used to disable buffering
                    cmd = ['unbuffer'] + cmd
                    try:
                        self.logger.debug(f"Prepended unbuffer for line-buffering: {' '.join(cmd)}")
                    except Exception:
                        pass
            except Exception:
                pass

            proc = subprocess.Popen(cmd, stdin=stdin_src, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, text=True, start_new_session=True, close_fds=True)

            # Start a background thread to stream stderr lines to the logger in near-real-time
            stderr_thread = None
            def _stream_stderr(stderr_pipe, tname):
                # Read stderr in small chunks and log segments split on both CR and LF.
                try:
                    if stderr_pipe is None:
                        return
                    buf = ''
                    while True:
                        chunk = stderr_pipe.read(1024)
                        if chunk == '':
                            # EOF reached
                            break
                        if not chunk:
                            # No data available; continue
                            continue
                        buf += chunk
                        # Split on CR or LF while keeping partial tail in buf
                        parts = re.split(r"[\r\n]", buf)
                        # If last char was not newline/CR, last element is partial; keep it
                        if buf[-1] not in ('\n', '\r'):
                            buf = parts.pop()
                        else:
                            buf = ''

                        for part in parts:
                            s = part.strip()
                            if s:
                                try:
                                    self.logger.info(f"{tname} stderr: {s}")
                                except Exception:
                                    pass
                    # Flush any remaining buffer
                    if buf:
                        s = buf.strip()
                        if s:
                            try:
                                self.logger.info(f"{tname} stderr: {s}")
                            except Exception:
                                pass
                except Exception:
                    pass

            if proc.stderr:
                stderr_thread = threading.Thread(target=_stream_stderr, args=(proc.stderr, tool_name), daemon=True)
                stderr_thread.start()
                try:
                    self.logger.debug(f"Started stderr thread for {tool_name}, pid={proc.pid}")
                except Exception:
                    pass

            kept = 0
            total = 0
            with open(output_file, 'w') as outf:
                assert proc.stdout is not None
                for line in proc.stdout:
                    # Preserve comment/info lines from BLAST/DIAMOND stdout (they often start with '#')
                    if line.startswith('#'):
                        try:
                            self.logger.info(f"{tool_name} stdout: {line.strip()}")
                        except Exception:
                            pass
                        continue
                    total += 1
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 13:
                        # If line too short, skip
                        continue
                    try:
                        identity = float(fields[2])
                        qstart = int(fields[6])
                        qend = int(fields[7])
                        qlen = int(fields[12])
                        # bitscore is present in the standard outfmt used by BLAST/DIAMOND
                        bitscore = float(fields[11]) if len(fields) > 11 else None
                    except Exception:
                        continue

                    # Stage-1 query coverage — span-based, valid for all search modes.
                    # qstart/qend are always in the query's own coordinate system (NT for blastn
                    # and blastx; AA for blastp), so qlen and the span are in the same unit.
                    query_coverage = ((abs(qend - qstart) + 1) / qlen) * 100 if qlen else 0.0

                    if identity >= getattr(self, 'query_min_identity', 0.0) and query_coverage >= getattr(self, 'query_min_coverage', 0.0):
                        # If the user requested a minimum bitscore, enforce it here
                        if getattr(self, 'min_bitscore', None) is not None:
                            try:
                                if bitscore is None or bitscore < float(self.min_bitscore):
                                    continue
                            except Exception:
                                # If bitscore parsing fails, be conservative and skip
                                continue
                        outf.write(line)
                        kept += 1

            # Wait for process to finish; stderr thread will exit when pipe closes
            if proc:
                proc.wait()

            if gzip_proc:
                try:
                    gzip_proc.stdout.close()
                    gzip_proc.wait()
                except Exception:
                    pass

            # Ensure stderr thread has finished
            if stderr_thread and stderr_thread.is_alive():
                try:
                    stderr_thread.join(timeout=1.0)
                except Exception:
                    pass

            if proc.returncode != 0:
                self.logger.error(f"{tool_name} failed with return code {proc.returncode}")
                # If extra parameters were supplied for this tool, warn the user that
                # the failure may be caused by incompatible/invalid parameters.
                try:
                    lowered = tool_name.lower()
                    for key, params in (self.extra_tool_params or {}).items():
                        if key and key.lower() in lowered:
                            self.logger.warning("\n" + "*" * 80)
                            self.logger.warning(f"WARNING: Additional parameters provided for tool '{key}': {params}")
                            self.logger.warning("These parameters may be incompatible with the tool and could have caused the failure.")
                            self.logger.warning("Please re-run without the extra params or verify they are valid for this tool.")
                            self.logger.warning("" + "*" * 80 + "\n")
                            break
                except Exception:
                    pass
                return False

            self.logger.info(f"Filtered {tool_name} output {output_file}: kept {kept}/{total} hits (identity>={self.query_min_identity}, qcov>={self.query_min_coverage}%)")
            return True

        except FileNotFoundError:
            self.logger.error(f"{tool_name} executable not found. Is it in your PATH?")
            return False
        except TypeError as e:
            # Common cause: passing a bool where a path was expected
            self.logger.error(f"TypeError running {tool_name} stream/filter: {e} -- cmd types: {[type(x) for x in cmd]}, fasta_path_str type: {type(fasta_path_str)}")
            return False
        except Exception as e:
            self.logger.error(f"Error running {tool_name} stream/filter: {e}")
            return False
        finally:
            # Close and terminate spawned processes
            try:
                if proc:
                    try:
                        if getattr(proc, 'stderr', None):
                            proc.stderr.close()
                    except Exception:
                        pass
                    try:
                        if getattr(proc, 'stdout', None):
                            proc.stdout.close()
                    except Exception:
                        pass
                    try:
                        self._terminate_proc(proc, tool_name)
                    except Exception:
                        pass
            except Exception:
                pass

            try:
                if gzip_proc:
                    try:
                        if getattr(gzip_proc, 'stdout', None):
                            gzip_proc.stdout.close()
                    except Exception:
                        pass
                    try:
                        self._terminate_proc(gzip_proc, 'gzip')
                    except Exception:
                        pass
            except Exception:
                pass

    def _split_fasta_chunks(self, fasta_path: str, target_chunks: Optional[int] = None) -> List[Path]:
        """Split a FASTA file into multiple chunk files under a temporary directory.
        Returns a list of Path objects for the created chunk FASTA files.
        Handles gzipped input transparently. Splits approximately by bytes, breaking between sequences.
        """
        src = str(fasta_path)
        is_gz = src.endswith(('.gz', '.gzip'))

        # Create a chunks directory under the user-specified temp directory to hold pieces
        # Keep chunk files inside the same temp directory used for other temporary files
        ts = int(time.time())
        chunks_dir = Path(str(self.temp_directory)) / f'genefior_chunks_{ts}'
        chunks_dir.mkdir(parents=True, exist_ok=True)
        chunk_paths: List[Path] = []

        opener = gzip.open if is_gz else open
        mode = 'rt' if is_gz else 'r'

        chunk_idx = 0
        bytes_written = 0
        out_fh = None
        current_chunk_path = None

        try:
            # If target_chunks is provided, split sequences evenly by count
            if target_chunks and target_chunks > 1:
                # Count total sequences first
                total_seqs = 0
                with opener(src, mode) as inf_count:
                    for line in inf_count:
                        if line.startswith('>'):
                            total_seqs += 1

                if total_seqs == 0:
                    return []

                per_chunk = max(1, (total_seqs + target_chunks - 1) // target_chunks)
                # Now iterate and write per_chunk sequences
                seqs_in_current = 0
                with opener(src, mode) as inf:
                    for line in inf:
                        if line.startswith('>') and seqs_in_current >= per_chunk:
                            # rotate chunk
                            if out_fh is not None:
                                out_fh.close()
                            out_fh = None
                            chunk_idx += 1
                            seqs_in_current = 0

                        if out_fh is None:
                            current_chunk_path = chunks_dir / f'chunk_{chunk_idx:04d}.fasta'
                            out_fh = open(current_chunk_path, 'w')
                            chunk_paths.append(current_chunk_path)

                        out_fh.write(line)
                        if line.startswith('>'):
                            seqs_in_current += 1

                if out_fh is not None:
                    out_fh.close()

            else:
                with opener(src, mode) as inf:
                    for line in inf:
                        # If header line and current chunk size exceeded threshold, rotate chunk
                        if line.startswith('>') and out_fh is not None and bytes_written >= self.max_fasta_chunk_bytes:
                            out_fh.close()
                            out_fh = None
                            chunk_idx += 1
                            bytes_written = 0

                        if out_fh is None:
                            current_chunk_path = chunks_dir / f'chunk_{chunk_idx:04d}.fasta'
                            out_fh = open(current_chunk_path, 'w')
                            chunk_paths.append(current_chunk_path)

                        out_fh.write(line)
                        bytes_written += len(line.encode('utf-8'))

            if out_fh is not None:
                out_fh.close()

        except Exception as e:
            # Clean up partial chunks on error
            try:
                if out_fh is not None:
                    out_fh.close()
            except Exception:
                pass
            for p in chunk_paths:
                try:
                    os.remove(str(p))
                except Exception:
                    pass
            try:
                os.rmdir(str(chunks_dir))
            except Exception:
                pass
            raise

        return chunk_paths

    def _iter_fasta_chunks(self, fasta_path: str, target_chunks: Optional[int] = None) -> Iterator[Path]:
        """Generator that yields chunk file Paths as they are created.
        This creates chunk files under self.temp_directory/genefior_chunks_<ts> and yields each path
        so callers can submit jobs immediately without creating all chunks up-front.
        """
        src = str(fasta_path)
        is_gz = src.endswith(('.gz', '.gzip'))

        ts = int(time.time())
        chunks_dir = Path(str(self.temp_directory)) / f'genefior_chunks_{ts}'
        chunks_dir.mkdir(parents=True, exist_ok=True)

        opener = gzip.open if is_gz else open
        mode = 'rt' if is_gz else 'r'

        chunk_idx = 0
        bytes_written = 0
        out_fh = None
        current_chunk_path = None

        try:
            if target_chunks and target_chunks > 1:
                # Count total sequences first
                total_seqs = 0
                with opener(src, mode) as inf_count:
                    for line in inf_count:
                        if line.startswith('>'):
                            total_seqs += 1

                if total_seqs == 0:
                    return

                per_chunk = max(1, (total_seqs + target_chunks - 1) // target_chunks)
                seqs_in_current = 0
                with opener(src, mode) as inf:
                    for line in inf:
                        if line.startswith('>') and seqs_in_current >= per_chunk:
                            if out_fh is not None:
                                out_fh.close()
                                yield current_chunk_path
                            out_fh = None
                            chunk_idx += 1
                            seqs_in_current = 0

                        if out_fh is None:
                            current_chunk_path = chunks_dir / f'chunk_{chunk_idx:06d}.fasta'
                            out_fh = open(current_chunk_path, 'w')

                        out_fh.write(line)
                        if line.startswith('>'):
                            seqs_in_current += 1

                    if out_fh is not None:
                        out_fh.close()
                        yield current_chunk_path

            else:
                with opener(src, mode) as inf:
                    for line in inf:
                        if line.startswith('>') and out_fh is not None and bytes_written >= self.max_fasta_chunk_bytes:
                            out_fh.close()
                            yield current_chunk_path
                            out_fh = None
                            chunk_idx += 1
                            bytes_written = 0

                        if out_fh is None:
                            current_chunk_path = chunks_dir / f'chunk_{chunk_idx:06d}.fasta'
                            out_fh = open(current_chunk_path, 'w')
                        out_fh.write(line)
                        bytes_written += len(line.encode('utf-8'))

                    if out_fh is not None:
                        out_fh.close()
                        yield current_chunk_path

        except Exception:
            try:
                if out_fh is not None:
                    out_fh.close()
            except Exception:
                pass
            # cleanup partial chunks on error
            try:
                for p in chunks_dir.iterdir():
                    try:
                        p.unlink()
                    except Exception:
                        pass
                chunks_dir.rmdir()
            except Exception:
                pass
            raise


    def _write_fasta_outputs(self, database: str, tool_name: str, detected_genes: Set[str],
                             gene_reads: dict, all_reads: dict):
        # Method to write FASTA files for mapped reads.
        evidence_genes = {
            gene
            for gene, gene_calls in self.evidence_calls.get(database, {}).items()
            if tool_name in gene_calls and gene_calls[tool_name].evidence_present
        }
        exact_genes = {
            gene
            for gene, gene_calls in self.evidence_calls.get(database, {}).items()
            if tool_name in gene_calls
            and gene_calls[tool_name].exact_allele_detected
        }
        candidate_genes = {
            gene
            for gene, gene_calls in self.evidence_calls.get(database, {}).items()
            if tool_name in gene_calls
            and gene_calls[tool_name].candidate_allele_detected
        }
        def sanitise_gene_name(gene: str) -> str:
            safe_gene = gene.replace('|', '_').replace('/', '_').replace(':','_').replace('-','_')
            if not hasattr(self, 'gene_name_changes'):
                self.gene_name_changes = []
            self.gene_name_changes.append((gene, safe_gene))
            return safe_gene

        def determine_read_type(read_name: str) -> str:
            # Determine if read is R1 or R2 based on suffix
            if read_name.endswith('_R1') or '/1' in read_name:
                return ' [R1]'
            elif read_name.endswith('_R2') or '/2' in read_name:
                return ' [R2]'
            return ''

        def count_read_types(read_names: set) -> tuple:
            # Count R1 and R2 reads in a set
            r1_count = sum(1 for rn in read_names if rn.endswith('_R1') or '/1' in rn)
            r2_count = sum(1 for rn in read_names if rn.endswith('_R2') or '/2' in rn)
            return r1_count, r2_count


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
                            read_type = determine_read_type(read_name)
                            fasta_out.write(f">{read_name}{read_type}\n{seq}\n")
                            count += 1
                if getattr(self, "verbose", True):
                    r1_count, r2_count = count_read_types(read_names)
                    self.logger.info(f"  FASTA file: {fasta_path} ({count} reads: {r1_count} R1, {r2_count} R2)")


        elif self.report_fasta in ('detected', 'evidence', 'candidate', 'exact'):
            selected_genes = (
                exact_genes
                if self.report_fasta == 'exact'
                else candidate_genes
                if self.report_fasta == 'candidate'
                else evidence_genes
            )
            if getattr(self, "verbose", True):
                self.logger.info(
                    f"Writing FASTA files for threshold-passing reads mapped to "
                    f"{self.report_fasta} genes in {database}...")
            for gene in selected_genes:
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
                            read_type = determine_read_type(read_name)
                            fasta_out.write(f">{read_name}{read_type}\n{seq}\n")
                            count += 1
                if getattr(self, "verbose", True):
                    r1_count, r2_count = count_read_types(read_names)
                    self.logger.info(f"  FASTA file: {fasta_path} ({count} reads: {r1_count} R1, {r2_count} R2)")


        elif self.report_fasta in (
                'detected-all', 'evidence-all', 'candidate-all',
                'exact-all'):
            selected_genes = (
                exact_genes
                if self.report_fasta == 'exact-all'
                else candidate_genes
                if self.report_fasta == 'candidate-all'
                else evidence_genes
            )
            if getattr(self, "verbose", True):
                self.logger.info(
                    f"Writing FASTA files for all reads mapped to "
                    f"{self.report_fasta} genes in {database}..."
                )
            for gene in selected_genes:
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
                            read_type = determine_read_type(read_name)
                            fasta_out.write(f">{read_name}{read_type}\n{seq}\n")
                            count += 1

                if getattr(self, "verbose", True):
                    r1_count, r2_count = count_read_types(read_names)
                    self.logger.info(f"  FASTA file: {fasta_path} ({count} reads: {r1_count} R1, {r2_count} R2)")

    def run_blast(self, db_path: str, database: str, mode: str) -> Tuple[bool, Set[str]]:
        """Run BLAST using the requested blast tool.
        mode may be one of: 'dna'|'blastn', 'protein'|'blastx', or 'blastp'.
        If input FASTA is gzipped, stream decompressed data to BLAST via stdin (uses `-query -`)."""
        # Normalise mode into concrete blast program name
        m = str(mode).lower() if mode is not None else 'dna'
        if m in ('dna', 'blastn'):
            blast_cmd = 'blastn'
            tool_name = 'BLASTn'
        elif m in ('protein', 'blastx'):
            blast_cmd = 'blastx'
            tool_name = 'BLASTx'
        elif m == 'blastp':
            blast_cmd = 'blastp'
            tool_name = 'BLASTp'
        else:
            self.logger.error(f"Invalid BLAST mode requested: {mode}")
            return False, set()
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

        # Determine if input is gzipped. Use a safe string when input_fasta is None.
        fasta_path_str = str(self.input_fasta) if getattr(self, 'input_fasta', None) else ''
        gz_input = bool(fasta_path_str and fasta_path_str.endswith(('.gz', '.gzip')))

        query_arg = '-' if gz_input else fasta_path_str
        # Build common command base
        # Build base command. Do NOT pass '-perc_identity' to BLAST programs that don't accept it
        # (blastp in particular raises an "Unknown argument: perc_identity" error).
        cmd = [
            blast_cmd,
            '-query', query_arg,
            '-db', db_path,
            '-out', '-',
            '-outfmt', outfmt_fields,
            '-num_threads', str(self.threads),
            '-qcov_hsp_perc', str(int(self.query_min_coverage)),
        ]
        # Only include an explicit evalue if the user supplied one; otherwise
        # allow BLAST to use its internal default
        if getattr(self, 'evalue', None) is not None:
            cmd.extend(['-evalue', str(self.evalue)])

        # Only add '-perc_identity' for BLAST programs that accept it (blastn).
        # We rely on post-filtering in parse_blast_results for identity thresholds
        # for other programs (blastx/blastp/diamond).
        try:
            if blast_cmd == 'blastn':
                cmd.extend(['-perc_identity', str(self.query_min_identity)])
        except Exception:
            pass

        # Add blastx-specific task if requested
        if blast_cmd == 'blastx':
            try:
                cmd.insert(3, '-task')
                cmd.insert(4, str(self.blastx_task))
            except Exception:
                pass

        # Append any user-specified extra parameters for this blast_cmd
        try:
            extra = (self.extra_tool_params or {}).get(blast_cmd)
            if extra:
                extra_args = shlex.split(str(extra))
                cmd.extend(extra_args)
        except Exception:
            pass

        success = False
        # If input FASTA file is large, split into chunks and run BLAST per-chunk to reduce memory/IO
        try:
            file_size = None
            try:
                # Only attempt to stat the file if we have a real path
                file_size = os.path.getsize(fasta_path_str) if fasta_path_str and os.path.isfile(fasta_path_str) else None
            except Exception:
                file_size = None

            # Determine if input FASTA is a regular file that can be chunked
            input_path_str = fasta_path_str if fasta_path_str else None
            input_is_regular_file = False
            try:
                input_is_regular_file = bool(input_path_str and input_path_str != '-' and os.path.isfile(input_path_str))
            except Exception:
                input_is_regular_file = False

            # If input is '-' (stdin) and user requested chunking, attempt to materialise stdin to a temp file
            stdin_materialised_path = None
            if input_path_str == '-' and self.chunk_jobs_user_specified and self.chunk_jobs > 1:
                # Only attempt to read stdin if data is piped (not a TTY)
                try:
                    if not sys.stdin.isatty():
                        self.logger.info("Materialising stdin to temporary FASTA for chunking as --chunk-jobs was requested...")
                        tf = tempfile.NamedTemporaryFile(delete=False, dir=str(self.temp_directory), prefix='genefior_stdin_', suffix='.fasta')
                        stdin_materialised_path = tf.name
                        # Read all stdin bytes and write
                        data = sys.stdin.buffer.read()
                        tf.write(data)
                        tf.flush()
                        tf.close()
                        input_is_regular_file = True
                        input_path_str = stdin_materialised_path
                        self.logger.info(f"Wrote piped stdin to temporary FASTA: {stdin_materialised_path}")
                    else:
                        self.logger.warning(f"User requested chunking (--chunk-jobs={self.chunk_jobs}) but stdin is a TTY; cannot materialise. Falling back to single-run BLAST.")
                except Exception as _e:
                    self.logger.warning(f"Failed to materialise stdin for chunking: {_e}; falling back to single-run BLAST")

            # Decide whether to chunk: either file is large enough OR user explicitly requested multiple chunk jobs
            force_by_jobs = (self.chunk_jobs_user_specified and self.chunk_jobs > 1 and input_is_regular_file)
            if (self.chunk_jobs_user_specified and self.chunk_jobs > 1) and not input_is_regular_file:
                self.logger.warning(f"User requested chunking (--chunk-jobs={self.chunk_jobs}) but input FASTA is not a regular file or is '-' — cannot chunk stdin. Falling back to single-run BLAST.")
            file_too_large = bool(file_size and file_size > self.max_fasta_chunk_bytes)
            self.logger.debug(f"Chunk decision vars: file_size={file_size}, max_fasta_chunk_bytes={self.max_fasta_chunk_bytes}, file_too_large={file_too_large}, chunk_jobs={self.chunk_jobs}, chunk_jobs_user_specified={getattr(self, 'chunk_jobs_user_specified', False)})")

            if file_too_large or force_by_jobs:
                # Chunking path (on-the-fly): create chunks incrementally and submit jobs as they are yielded
                self.logger.info(f"Chunking FASTA (on-the-fly): file_size={file_size}, threshold={self.max_fasta_chunk_bytes}, chunk_jobs={self.chunk_jobs}, force_by_jobs={force_by_jobs}")

                # Decide whether to request a target_chunks split (only when user forced by jobs and file not too large)
                target_chunks = self.chunk_jobs if (force_by_jobs and not file_too_large) else None

                # Determine concurrency and threads-per-job robustly.
                # Start with safe defaults
                concurrency = 1
                threads_per_job_list = [max(1, int(self.threads))]

                # If user requested a specific threads-per-job, compute a cap
                max_concurrency_by_threads = None
                if self.chunk_threads_per_job:
                    try:
                        t_per = max(1, int(self.chunk_threads_per_job))
                        max_concurrency_by_threads = max(1, int(self.threads) // t_per)
                    except Exception:
                        max_concurrency_by_threads = None

                # Determine desired concurrency: honour explicit --chunk-jobs if provided
                if getattr(self, 'chunk_jobs_user_specified', False) and self.chunk_jobs and self.chunk_jobs > 0:
                    desired = int(self.chunk_jobs)
                else:
                    # Auto derive from available threads: try to give a reasonable number of concurrent jobs
                    if int(self.threads) >= 8:
                        desired = max(1, int(self.threads) // 4)
                    elif int(self.threads) >= 4:
                        desired = 2
                    else:
                        desired = 1

                # Cap desired by thread-based cap if applicable
                if max_concurrency_by_threads is not None:
                    concurrency = min(desired, max_concurrency_by_threads)
                else:
                    concurrency = desired

                # Ensure at least 1
                concurrency = max(1, int(concurrency))

                # Compute threads-per-job distribution
                if self.chunk_threads_per_job:
                    t_per = max(1, int(self.chunk_threads_per_job))
                    # If requested threads-per-job exceed total threads, fall back to single concurrency
                    if t_per * concurrency > int(self.threads):
                        concurrency = max(1, int(self.threads) // t_per)
                        if concurrency < 1:
                            concurrency = 1
                    threads_per_job_list = [max(1, t_per)] * concurrency
                else:
                    # Evenly distribute available threads across concurrency slots
                    base = max(1, int(self.threads) // concurrency)
                    remainder = int(self.threads) - (base * concurrency)
                    threads_per_job_list = [base + (1 if i < remainder else 0) for i in range(concurrency)]

                self.logger.info(f"Starting on-the-fly chunking with concurrency={concurrency}, threads per job distribution={threads_per_job_list}, max_concurrency_by_threads={max_concurrency_by_threads}")

                part_outputs = []
                chunk_paths_created = []
                failed = False

                chunk_gen = self._iter_fasta_chunks(self.input_fasta, target_chunks=target_chunks)
                with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as executor:
                    running = {}
                    i = 0
                    for chunk in chunk_gen:
                        # record chunk
                        chunk_paths_created.append(chunk)

                        # Wait for a slot if we've reached concurrency capacity
                        while len(running) >= concurrency:
                            done, _ = concurrent.futures.wait(running.keys(), return_when=concurrent.futures.FIRST_COMPLETED)
                            for fut in done:
                                idx, chunk_output, chunk_path = running.pop(fut)
                                try:
                                    ok = fut.result()
                                except Exception as e:
                                    self.logger.error(f"Chunk {idx} exception: {e}")
                                    ok = False
                                if not ok:
                                    self.logger.error(f"Chunk {idx} failed: {chunk_path}")
                                    failed = True
                                    break
                                part_outputs.append(chunk_output)
                                # remove chunk file unless preserving
                                if not getattr(self, 'preserve_chunks', False):
                                    try:
                                        os.remove(str(chunk_path))
                                    except Exception:
                                        pass
                            if failed:
                                break
                        if failed:
                            break

                        slot = i % concurrency
                        threads_for_this = threads_per_job_list[slot]
                        chunk_output = output_file.parent / f"{output_file.stem}.part{i:04d}.tsv"
                        # Build per-chunk command for the concrete blast_cmd
                        per_cmd = [
                            blast_cmd,
                            '-query', str(chunk),
                            '-db', db_path,
                            '-out', '-',
                            '-outfmt', outfmt_fields,
                            '-num_threads', str(threads_for_this),
                            '-qcov_hsp_perc', str(int(self.query_min_coverage)),
                        ]
                        if getattr(self, 'evalue', None) is not None:
                            per_cmd.extend(['-evalue', str(self.evalue)])
                        if blast_cmd == 'blastx':
                            per_cmd.insert(3, '-task')
                            per_cmd.insert(4, str(self.blastx_task))
                        try:
                            extra = (self.extra_tool_params or {}).get(blast_cmd)
                            if extra:
                                per_cmd.extend(shlex.split(str(extra)))
                        except Exception:
                            pass

                        self.logger.info(f"Submitting chunk {i+1} (slot {slot}, threads={threads_for_this}): {chunk}")
                        fut = executor.submit(self._stream_filter_hits, per_cmd, chunk_output, f"{database} - {tool_name} (chunk {i})", False, str(chunk))
                        running[fut] = (i, chunk_output, chunk)
                        i += 1

                    # After generator exhausted, wait for remaining running futures to complete
                    if not failed:
                        while running:
                            done, _ = concurrent.futures.wait(running.keys(), return_when=concurrent.futures.FIRST_COMPLETED)
                            for fut in done:
                                idx, chunk_output, chunk_path = running.pop(fut)
                                try:
                                    ok = fut.result()
                                except Exception as e:
                                    self.logger.error(f"Chunk {idx} exception: {e}")
                                    ok = False
                                if not ok:
                                    self.logger.error(f"Chunk {idx} failed: {chunk_path}")
                                    failed = True
                                    break
                                part_outputs.append(chunk_output)
                                if not getattr(self, 'preserve_chunks', False):
                                    try:
                                        os.remove(str(chunk_path))
                                    except Exception:
                                        pass

                    # If failed, cleanup and abort
                    if failed:
                        for p in part_outputs:
                            try:
                                os.remove(str(p))
                            except Exception:
                                pass
                        try:
                            if chunk_paths_created:
                                chunks_dir = Path(chunk_paths_created[0]).parent
                                for p in chunk_paths_created:
                                    try:
                                        os.remove(str(p))
                                    except Exception:
                                        pass
                                try:
                                    os.rmdir(str(chunks_dir))
                                except Exception:
                                    pass
                        except Exception:
                            pass
                        return False, set()

                num_chunks = len(chunk_paths_created)
                self.logger.info(f"Completed on-the-fly chunking: created {num_chunks} chunks, combined parts: {len(part_outputs)}")

                # Concatenate part outputs to final output_file
                try:
                    with open(output_file, 'w') as outf:
                        for p in part_outputs:
                            if os.path.isfile(str(p)):
                                with open(str(p), 'r') as inf:
                                    for ln in inf:
                                        outf.write(ln)
                except Exception as e:
                    self.logger.error(f"Failed to write concatenated BLAST output: {e}")
                    return False, set()

                # Clean up part outputs and chunk files unless user requested preservation
                if not getattr(self, 'preserve_chunks', False):
                    for p in part_outputs:
                        try:
                            os.remove(str(p))
                        except Exception:
                            pass
                    try:
                        if chunk_paths_created:
                            chunks_dir = Path(chunk_paths_created[0]).parent
                            for p in chunk_paths_created:
                                try:
                                    os.remove(str(p))
                                except Exception:
                                    pass
                            try:
                                os.rmdir(str(chunks_dir))
                            except Exception:
                                pass
                    except Exception:
                        pass
                else:
                    if chunk_paths_created:
                        self.logger.info(f"Preserving chunk files in {Path(chunk_paths_created[0]).parent} as requested (--preserve-chunks)")
                    else:
                        self.logger.info(f"Preserve-chunks requested but no chunk files were created")

                success = True
            else:
                # Log the command we will run and whether input will be gzipped-streamed
                try:
                    self.logger.info(f"Parameters for {database} - {tool_name}: {' '.join(str(x) for x in cmd)}")
                    self.logger.info(f"gz_input={gz_input}, fasta={fasta_path_str}")
                except Exception:
                    pass
                # Stream BLAST output to our filter (handles gzipped input when query is '-')
                success = self._stream_filter_hits(cmd, output_file, f"{database} - {tool_name}", gz_input, fasta_path_str)

        except Exception as e:
            self.logger.error(f"Error during BLAST chunking/run: {e}")
            return False, set()

        detected = set()
        if success:
            detected, gene_reads = self.parse_blast_results(output_file, database, tool_name)
            self.write_tool_stats(database, tool_name, gene_reads)
        return success, detected

    def run_diamond(self, db_path: str, database: str, query_mode: str = 'blastx') -> Tuple[bool, Set[str]]:
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
        tool_name = f"DIAMOND-{str(query_mode).upper()}"

        fasta_path_str = str(self.input_fasta)
        gz_input = fasta_path_str.endswith(('.gz', '.gzip'))
        if gz_input and not self.check_gzip(fasta_path_str): # If input is gzipped, check integrity first
            return False, set()
        else: # Run DIAMOND normally
            params = (self.tool_sensitivity_params or {}).get('diamond', None)
            sensitivity = params['sensitivity'] if params and 'sensitivity' in params else None

            cmd = [
                'diamond', str(query_mode),
                '-q', str(self.input_fasta),
                '-d', db_path,
                '-o', str(output_file),
                '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen',
                '--id', str(self.query_min_identity),
                '--query-cover', str(self.query_min_coverage),
                #'-e', str(self.evalue),
                '-p', str(self.threads)#,
                #'-k', '10'
            ]
            if sensitivity and sensitivity != 'default':
                cmd.append(sensitivity)

            # Add evalue for DIAMOND only if the user supplied one. Otherwise
            # leave DIAMOND to use its default e-value.
            try:
                if getattr(self, 'evalue', None) is not None:
                    cmd.extend(['-e', str(self.evalue)])
            except Exception:
                pass
            # Append any user-specified extra parameters for DIAMOND
            try:
                extra = (self.extra_tool_params or {}).get('diamond')
                if extra:
                    extra_args = shlex.split(str(extra))
                    cmd.extend(extra_args)
            except Exception:
                pass

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

        if self.sequence_type == 'Single-FASTA' or getattr(self, 'is_genes_fasta', False):
            flags = ['-f', '-U', str(self.input_fasta)]
        elif self.sequence_type == 'Paired-FASTQ':
            flags = ['-1', str(self.input_fastq[0]), '-2', str(self.input_fastq[1])]
        else:
            flags = []

        params = (self.tool_sensitivity_params or {}).get('bowtie2', None)
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

        # Append any user-specified extra parameters for Bowtie2
        try:
            extra = (self.extra_tool_params or {}).get('bowtie2')
            if extra:
                extra_args = shlex.split(str(extra))
                cmd.extend(extra_args)
        except Exception:
            pass

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
        # Generate per-base coverage for this mapping result (best-effort)
        try:
            self.generate_coverage_from_bam(sorted_bam_file, database, tool_name)
        except Exception:
            self.logger.debug(f"Coverage generation failed or skipped for {database} using {tool_name}")

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

        if self.sequence_type == 'Single-FASTA' or getattr(self, 'is_genes_fasta', False):
            flags = [ str(self.input_fasta)]
        elif self.sequence_type == 'Paired-FASTQ':
            flags = [ str(self.input_fastq[0]), str(self.input_fastq[1])]
        else:
            flags = []

        cmd = [
            'bwa', 'mem',
            '-t', str(self.threads),
            '-o', str(sam_file),
            db_path,
            ] + flags + [
        ]

        # Append any user-specified extra parameters for BWA
        try:
            extra = (self.extra_tool_params or {}).get('bwa')
            if extra:
                extra_args = shlex.split(str(extra))
                # Insert extra args after the 'mem' subcommand and before db/flags
                # For simplicity append at end which is usually acceptable
                cmd.extend(extra_args)
        except Exception:
            pass

        # Run BWA and write output to SAM file
        try:
            success = self.run_command(cmd, f"{database} - {tool_name}")
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
        # Generate per-base coverage for this mapping result (best-effort)
        try:
            self.generate_coverage_from_bam(sorted_bam, database, tool_name)
        except Exception:
            self.logger.debug(f"Coverage generation failed or skipped for {database} using {tool_name}")

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


        if self.sequence_type == 'Single-FASTA' or getattr(self, 'is_genes_fasta', False):
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
            # Append any user-specified extra parameters for Minimap2
            try:
                extra = (self.extra_tool_params or {}).get('minimap2')
                if extra:
                    extra_args = shlex.split(str(extra))
                    cmd.extend(extra_args)
            except Exception:
                pass

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
        # Generate per-base coverage for this mapping result (best-effort)
        try:
            self.generate_coverage_from_bam(sorted_bam, database, tool_name)
        except Exception:
            self.logger.debug(f"Coverage generation failed or skipped for {database} using {tool_name}")

        return success, detected




    # Column name variants accepted for the optional must-flag field.
    _MUST_FLAG_COL_VARIANTS = ('Must flag', 'Must Flag', 'must_flag', 'must flag', 'MUST FLAG', 'MustFlag')

    def _load_hmmer_annotations(self, annotations_path: str) -> Tuple[dict, dict]:
        """Load a HMMER annotations CSV.

        The CSV must have an ``ID`` column (gene/profile name) and optionally a
        ``Description`` column and a ``Must flag`` column (case/space-insensitive).

        The ``Must flag`` column is **optional** — databases that have no concept
        of priority flagging simply omit it and the feature is transparently disabled.

        Returns
        -------
        annotations : dict
            {gene_id: {'description': str, 'must_flag': bool}}
        meta : dict
            Metadata about the loaded file.  Currently:
              has_must_flag (bool) — whether the source CSV contained a must-flag column.
        """
        annotations: dict = {}
        meta: dict = {'has_must_flag': False}
        try:
            with open(annotations_path, 'r', newline='', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                # Detect must-flag column once using the actual CSV fieldnames
                mf_col = None
                if reader.fieldnames:
                    for candidate in self._MUST_FLAG_COL_VARIANTS:
                        if candidate in reader.fieldnames:
                            mf_col = candidate
                            break
                meta['has_must_flag'] = mf_col is not None

                for row in reader:
                    gene_id = row.get('ID', '').strip()
                    if not gene_id:
                        continue
                    description = row.get('Description', '').strip()
                    if mf_col is not None:
                        must_flag = row.get(mf_col, '').strip().upper() in ('TRUE', 'YES', '1')
                    else:
                        must_flag = False
                    annotations[gene_id] = {'description': description, 'must_flag': must_flag}

            n_must = sum(1 for v in annotations.values() if v['must_flag'])
            flag_note = f", {n_must} must-flag" if meta['has_must_flag'] else " (no must-flag column)"
            self.logger.info(f"  Loaded {len(annotations)} HMMER annotations from {annotations_path}{flag_note}")
        except Exception as e:
            self.logger.warning(f"Could not load HMMER annotations from {annotations_path}: {e}")
        return annotations, meta

    def _write_hmmer_annotations_report(self, database: str, tool_name: str):
        """Write a TSV report mapping detected HMMER gene IDs to their annotations.

        Produces:
        - A full annotations TSV (all genes in the annotation file, detected or not).
        - Console/log warnings when any 'Must Flag' genes are detected.

        The 'Must_Flag' column is only included in the TSV — and must-flag alerts are
        only emitted — when the source annotation CSV actually contained a must-flag
        column.  This makes the method universally usable for any HMM database, not
        only the ibbis biorisk database.
        """
        annotations = self.hmmer_annotations.get(database, {})
        if not annotations:
            return

        meta = self.hmmer_annotations_meta.get(database, {})
        has_must_flag = meta.get('has_must_flag', False)

        report_file = self.output_dir / f"{database}_{tool_name}_annotations.tsv"
        detected_in_tool = {
            gene for gene, tools in self.detections[database].items()
            if tools.get(tool_name)
        }

        try:
            with open(report_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                # Only include Must_Flag column when the source CSV had one
                headers = ['Gene_ID', 'Description']
                if has_must_flag:
                    headers.append('Must_Flag')
                headers.append('Detected')
                writer.writerow(headers)

                for gene_id, info in sorted(annotations.items()):
                    detected = 1 if gene_id in detected_in_tool else 0
                    row = [gene_id, info.get('description', '')]
                    if has_must_flag:
                        row.append('TRUE' if info.get('must_flag') else 'FALSE')
                    row.append(detected)
                    writer.writerow(row)
            self.logger.info(f"  HMMER annotations report: {report_file}")
        except Exception as e:
            self.logger.warning(f"Could not write HMMER annotations report: {e}")

        # Compute annotation summary counts
        n_total_detected = len(detected_in_tool)
        n_annotated_detected = sum(1 for g in detected_in_tool if g in annotations)
        n_unannotated_detected = n_total_detected - n_annotated_detected

        # Build base summary line
        summary = (
            f"  Detection summary [{tool_name}]: {n_total_detected} gene(s) detected — "
            f"{n_annotated_detected} have annotation entries, {n_unannotated_detected} have none"
        )

        if has_must_flag:
            n_must_flag_in_db = sum(1 for info in annotations.values() if info.get('must_flag'))
            must_flag_detected = sorted([
                (gid, annotations[gid]['description'])
                for gid in detected_in_tool
                if annotations.get(gid, {}).get('must_flag')
            ])
            n_must_flag_detected = len(must_flag_detected)
            mf_enabled = getattr(self, 'hmmer_must_flag', True)
            summary += f"; {n_must_flag_detected} of {n_must_flag_in_db} must-flag gene(s) detected"
            if not mf_enabled:
                summary += " [must-flag override DISABLED]"
            self.logger.info(summary + ".")

            if must_flag_detected and mf_enabled:
                self.logger.warning("!" * 70)
                self.logger.warning(
                    f"  *** MUST-FLAG ALERT: {n_must_flag_detected} of {n_total_detected} detected gene(s) "
                    f"are flagged as high-priority genes [{database} / {tool_name}] ***"
                )
                for gid, desc in must_flag_detected:
                    self.logger.warning(f"    ⚑  {gid}: {desc}")
                self.logger.warning("!" * 70)
            elif must_flag_detected and not mf_enabled:
                self.logger.info(
                    f"  {n_must_flag_detected} must-flag gene(s) detected "
                    f"(must-flag override is disabled — no threshold bypass applied)."
                )
            else:
                self.logger.info(f"  No must-flag genes detected in {database} [{tool_name}].")
        else:
            self.logger.info(summary + ".")

    # Standard genetic code codon table
    _CODON_TABLE: Dict[str, str] = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    @staticmethod
    def _revcomp(seq: str) -> str:
        """Return the reverse complement of a DNA sequence."""
        comp = str.maketrans('ACGTacgt', 'TGCAtgca')
        return seq.translate(comp)[::-1]

    def _translate_frame(self, dna: str, offset: int, reverse: bool) -> str:
        """Translate one reading frame. Returns amino-acid string (no stop markers)."""
        if reverse:
            dna = self._revcomp(dna)
        dna = dna[offset:].upper()
        aa = []
        table = self._CODON_TABLE
        for i in range(0, len(dna) - 2, 3):
            codon = dna[i:i+3]
            aa.append(table.get(codon, 'X'))
        return ''.join(aa)

    def _prepare_protein_for_hmmsearch(self, nucleotide_fasta: str, database: str,
                                       prep_mode: str, min_aa: int = 15) -> str:
        """Prepare a protein FASTA for hmmsearch from various input types.

        prep_mode:
          'protein'  - input is already AA; just strip '*' stop markers and decompress if gz
          'cds_dna'  - input is DNA CDS; translate frame +1 only
          'reads'    - unassembled reads; 6-frame translate, split on stops, keep >= min_aa aa
        Returns path to a (possibly temporary) protein FASTA file ready for hmmsearch.
        """
        suffix = '.faa'
        _tmp_fd, _tmp_path = tempfile.mkstemp(
            suffix=suffix,
            dir=str(self.temp_directory),
            prefix=f"hmmer_{database}_prep_",
        )
        os.close(_tmp_fd)

        frame_labels = ['+1', '+2', '+3', '-1', '-2', '-3']

        def _iter_fasta(path: str):
            """Yield (header, seq) pairs from a plain or gzipped FASTA."""
            opener = gzip.open(path, 'rt') if path.endswith(('.gz', '.gzip')) else open(path, 'r')
            with opener as fh:
                header, seq_parts = None, []
                for line in fh:
                    line = line.rstrip('\n')
                    if line.startswith('>'):
                        if header is not None:
                            yield header, ''.join(seq_parts)
                        header = line[1:]
                        seq_parts = []
                    else:
                        seq_parts.append(line)
                if header is not None:
                    yield header, ''.join(seq_parts)

        seq_count = 0
        with open(_tmp_path, 'w') as out_fh:
            for header, seq in _iter_fasta(nucleotide_fasta):
                seq_id = header.split()[0]

                if prep_mode == 'protein':
                    # Already protein — just strip stop-codon asterisks
                    clean = seq.replace('*', '').strip()
                    if len(clean) >= min_aa:
                        out_fh.write(f">{seq_id}\n{clean}\n")
                        seq_count += 1

                elif prep_mode == 'cds_dna':
                    # DNA CDS — translate frame +1 only
                    aa = self._translate_frame(seq, 0, False).rstrip('*').replace('*', '')
                    if len(aa) >= min_aa:
                        out_fh.write(f">{seq_id}\n{aa}\n")
                        seq_count += 1

                else:  # 'reads' — 6-frame translation
                    for fi, (offset, reverse) in enumerate([(0, False), (1, False), (2, False),
                                                             (0, True),  (1, True),  (2, True)]):
                        aa_full = self._translate_frame(seq, offset, reverse)
                        # Split on stop codons and keep fragments >= min_aa
                        fragments = aa_full.split('*')
                        for frag_i, frag in enumerate(fragments):
                            frag = frag.strip()
                            if len(frag) >= min_aa:
                                frag_id = f"{seq_id}|frame={frame_labels[fi]}|frag={frag_i}"
                                out_fh.write(f">{frag_id}\n{frag}\n")
                                seq_count += 1

        self.logger.info(f"  HMMER protein prep ({prep_mode}): wrote {seq_count} sequences to {_tmp_path}")
        return _tmp_path

    def run_hmmer(self, db_path: str, database: str, mode: str) -> Tuple[bool, Set[str]]:
        # Run HMMER search.
        if not db_path:
            return False, set()

        hmmer_cmd = 'nhmmer' if mode == 'dna' else 'hmmsearch'
        output_file = self.raw_dir / f"{database}_{hmmer_cmd}_results.tbl"
        domtbl_file = self.raw_dir / f"{database}_{hmmer_cmd}_domtbl.txt"
        tool_name = f"HMMER-{mode.upper()}"

        _tmp_fasta = None

        if mode == 'protein':
            # Determine prep mode from input type
            is_genes = getattr(self, 'is_genes_fasta', False)
            detected_gt = (getattr(self, 'detected_genes_type', None)
                           or getattr(self, 'genes_type', None))
            if is_genes and detected_gt == 'protein':
                prep_mode = 'protein'
            elif is_genes and detected_gt == 'dna':
                prep_mode = 'cds_dna'
            elif is_genes:
                prep_mode = 'protein'   # unknown genes type — assume protein
            else:
                prep_mode = 'reads'     # Paired-FASTQ / Single-FASTA

            self.logger.info(f"  HMMER protein mode: input prep_mode='{prep_mode}'")
            try:
                _tmp_path = self._prepare_protein_for_hmmsearch(
                    str(self.input_fasta), database, prep_mode)
                fasta_path_for_hmmer = _tmp_path
                _tmp_fasta = _tmp_path
            except Exception as e:
                self.logger.warning(f"  Protein prep for HMMER failed: {e}; falling back to raw input.")
                fasta_path_for_hmmer = str(self.input_fasta)
        else:
            # DNA mode (nhmmer) — decompress gz if needed, pass through otherwise
            fasta_path_for_hmmer = str(self.input_fasta)
            if str(self.input_fasta).endswith(('.gz', '.gzip')):
                try:
                    _tmp_fd, _tmp_path = tempfile.mkstemp(
                        suffix='.fna',
                        dir=str(self.temp_directory),
                        prefix=f"hmmer_{database}_",
                    )
                    os.close(_tmp_fd)
                    with gzip.open(str(self.input_fasta), 'rt') as gz_in, open(_tmp_path, 'w') as f_out:
                        for chunk in iter(lambda: gz_in.read(1 << 20), ''):
                            f_out.write(chunk)
                    fasta_path_for_hmmer = _tmp_path
                    _tmp_fasta = _tmp_path
                    self.logger.debug(f"  Decompressed input FASTA for nhmmer: {_tmp_path}")
                except Exception as e:
                    self.logger.warning(f"  Failed to decompress FASTA for nhmmer: {e}; attempting with compressed input.")

        cmd = [
            hmmer_cmd,
            '--tblout', str(output_file),
            '--domtblout', str(domtbl_file),
            '--cpu', str(self.threads),
            db_path,
            fasta_path_for_hmmer,
        ]

        # ── HMMER threshold handling ──────────────────────────────────────────────
        _tmode = getattr(self, 'hmmer_threshold_mode', 'evalue')
        if _tmode in ('tc', 'ga', 'nc'):
            # Use per-profile trusted/gathering/noise cutoffs embedded in the HMM database.
            # This is the gold standard for curated databases (e.g. ibbis biorisk.hmm).
            # E-value flags are NOT added; parse_hmmer_results will accept all reported hits.
            cmd.insert(3, f'--cut_{_tmode}')
            self.logger.info(f"  HMMER threshold mode: per-profile {_tmode.upper()} cutoffs (--cut_{_tmode})")
        else:
            # E-value mode: prefer --hmmer-evalue, then global -e, then log a warning
            _ev = getattr(self, 'hmmer_evalue', None)
            if _ev is None:
                _ev = getattr(self, 'evalue', None)
            if _ev is not None:
                cmd.insert(3, '-E')
                cmd.insert(4, str(_ev))
                self.logger.info(f"  HMMER threshold mode: E-value = {_ev}")
            else:
                self.logger.warning(
                    "  HMMER threshold mode: no explicit E-value set — using HMMER default (E=10, "
                    "very permissive). Consider --hmmer-evalue 1e-5 for Genes-FASTA, "
                    "1e-3 for reads, or --hmmer-threshold-mode tc if your HMM database "
                    "has per-profile trusted cutoffs."
                )

        # Append any user-specified extra parameters for HMMER/nhmmer/hmmsearch
        try:
            extra = None
            if (self.extra_tool_params or {}).get('hmmer'):
                extra = (self.extra_tool_params or {}).get('hmmer')
            elif (self.extra_tool_params or {}).get(hmmer_cmd):
                extra = (self.extra_tool_params or {}).get(hmmer_cmd)
            if extra:
                extra_args = shlex.split(str(extra))
                cmd.extend(extra_args)
        except Exception:
            pass

        # Load annotations for this database if available
        annotations_path = (self.databases.get(database) or {}).get('hmmer_annotations')
        if annotations_path and database not in self.hmmer_annotations:
            _ann, _meta = self._load_hmmer_annotations(annotations_path)
            self.hmmer_annotations[database] = _ann
            self.hmmer_annotations_meta[database] = _meta

        success = self.run_command(cmd, f"{database} - {tool_name}")
        detected = set()
        if success:
            detected, gene_reads = self.parse_hmmer_results(output_file, domtbl_file, database, tool_name)
            try:
                self.write_tool_stats(database, tool_name, gene_reads=gene_reads)
            except Exception:
                self.logger.debug("write_tool_stats failed for HMMER output; continuing")
            # Write annotations report if annotations were loaded
            if database in self.hmmer_annotations and self.hmmer_annotations[database]:
                self._write_hmmer_annotations_report(database, tool_name)

        # Clean up decompressed temp file
        if _tmp_fasta and not getattr(self, 'no_cleanup', False):
            try:
                os.remove(_tmp_fasta)
            except Exception:
                pass

        return success, detected

    def parse_blast_results(self, output_file: Path, database: str, tool_name: str,
                            count_only: bool = False, read_store=None) -> Set[str]:
        """Parse BLAST/DIAMOND tabular output and collect per-gene stats.

        Only alignments meeting identity and query-coverage thresholds are added.
        Gene detection is decided from coverage-at-depth, identity, and
        read-support thresholds.

        When count_only is True, per-gene read-name lists are replaced with
        integer counters. This is intended for GeneFior-Recompute on very large
        tabular outputs where FASTA reporting is not being requested.
        """
        detected_genes = set()
        gene_lengths = {}  # Store gene lengths
        count_only = bool(count_only and (not self.report_fasta or read_store is not None))
        if count_only:
            gene_reads = defaultdict(lambda: {'passing_count': 0, 'all_count': 0})
        else:
            gene_reads = defaultdict(lambda: {'passing': [], 'all': []})  # Track all reads per gene

        # Ensure structure exists for this database and tool (do not reset entire gene_stats)
        if database not in self.gene_stats:
            self.gene_stats[database] = defaultdict(lambda: defaultdict(GeneStats))
            self.detections.setdefault(database, defaultdict(lambda: defaultdict(bool)))
        elif self.gene_stats[database].get(tool_name) is None:
            self.gene_stats[database][tool_name] = defaultdict(GeneStats)

        if not output_file.exists():
            # No output present -> return empty detected set and empty gene_reads dict
            return set(), {}

        # Load all reads from input FASTA only when FASTA output is requested.
        # Recompute commonly runs without query FASTA and may process huge hit
        # tables, so avoid materialising reads unless the user explicitly needs
        # mapped-read FASTA files.
        if self.report_fasta and read_store is None and not hasattr(self, 'all_reads'):
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
                            read_name = line[1:].split(' ')[0].replace('\n','') # DIAMOND/BLAST read names split at space
                            seq_lines = []
                        else:
                            seq_lines.append(line.strip())
                    if read_name and seq_lines:
                        self.all_reads[read_name] = ''.join(seq_lines)
            except Exception as e:
                self.logger.error(f"Error reading FASTA file - Ignore if running Genefíor-Recompute without sequences")
        all_reads = getattr(self, 'all_reads', {})
        mapped_reads = 0
        passing_reads = 0
        malformed_lines = 0

        def process_query_group(read_name: str, hits: List[dict]):
            nonlocal mapped_reads, passing_reads
            if not hits:
                return

            qualifying_by_gene = defaultdict(list)
            raw_genes = set()
            for hit in hits:
                gene = hit['gene']
                raw_genes.add(gene)

                if hit['identity'] < self.query_min_identity:
                    continue
                if hit['query_coverage'] < self.query_min_coverage:
                    continue
                if getattr(self, 'evalue', None) is not None:
                    if hit['evalue'] is None or hit['evalue'] > float(self.evalue):
                        continue
                if getattr(self, 'min_bitscore', None) is not None:
                    if hit['bitscore'] is None or hit['bitscore'] < float(self.min_bitscore):
                        continue
                qualifying_by_gene[gene].append(hit)

            for gene in raw_genes:
                if count_only:
                    gene_reads[gene]['all_count'] += 1
                else:
                    gene_reads[gene]['all'].append(read_name)
                if read_store is not None:
                    read_store.add(
                        database, tool_name, gene, read_name,
                        passing=False,
                    )
            mapped_reads += len(raw_genes)

            best_score_by_gene: Dict[str, float] = {}
            for gene, gene_hits in qualifying_by_gene.items():
                if gene not in self.gene_stats[database][tool_name]:
                    self.gene_stats[database][tool_name][gene] = GeneStats(
                        gene_name=gene,
                        compact=count_only,
                    )
                stats = self.gene_stats[database][tool_name][gene]
                intervals = [
                    (
                        min(hit['sstart'], hit['send']) - 1,
                        max(hit['sstart'], hit['send']),
                    )
                    for hit in gene_hits
                ]
                identity_weights = [
                    max(1, abs(hit['send'] - hit['sstart']) + 1)
                    for hit in gene_hits
                ]
                aggregate_identity = (
                    sum(
                        hit['identity'] * weight
                        for hit, weight in zip(gene_hits, identity_weights)
                    )
                    / sum(identity_weights)
                )
                stats.add_intervals(
                    intervals,
                    aggregate_identity,
                    gene_lengths.get(gene, 0),
                )
                if count_only:
                    gene_reads[gene]['passing_count'] += 1
                else:
                    gene_reads[gene]['passing'].append(read_name)
                if read_store is not None:
                    read_store.add(
                        database, tool_name, gene, read_name,
                        passing=True,
                    )
                passing_reads += 1
                score = max(
                    (
                        hit['bitscore']
                        if hit['bitscore'] is not None
                        else hit['identity']
                    )
                    for hit in gene_hits
                )
                best_score_by_gene[gene] = score

            if not best_score_by_gene:
                return
            best_score = max(best_score_by_gene.values())
            tied_genes = {
                gene
                for gene, score in best_score_by_gene.items()
                if scores_tied(best_score, score, self.evidence_config)
            }
            for gene, score in best_score_by_gene.items():
                is_best = gene in tied_genes
                self.gene_stats[database][tool_name][gene].add_read_support(
                    mapped=gene in raw_genes,
                    passing=True,
                    best=is_best,
                    unique_best=is_best and len(tied_genes) == 1,
                    ambiguous_best=is_best and len(tied_genes) > 1,
                    high_confidence=is_best and len(tied_genes) == 1,
                    score=score,
                )

        try:
            with open(output_file, 'r') as f:
                current_read = None
                current_hits: List[dict] = []
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) < 14:
                        continue

                    try:
                        read_name = fields[0]  # qseqid
                        gene = fields[1]  # sseqid
                        identity = float(fields[2])  # pident
                        qstart = int(fields[6])  # query start
                        qend = int(fields[7])  # query end
                        sstart = int(fields[8])  # subject start
                        send = int(fields[9])  # subject end
                        qlen = int(fields[12])  # query length (added to output format)
                        slen = int(fields[13])  # subject length (added to output format)
                    except (TypeError, ValueError):
                        malformed_lines += 1
                        continue
                    # bitscore included in the standard outfmt
                    try:
                        bitscore = float(fields[11]) if len(fields) > 11 else None
                    except Exception:
                        bitscore = None
                    try:
                        evalue = float(fields[10]) if len(fields) > 10 else None
                    except Exception:
                        evalue = None

                    # Store gene length if available
                    if slen is not None:
                        gene_lengths[gene] = max(gene_lengths.get(gene, 0), slen)

                    # Query coverage: fraction of the query sequence spanned by this HSP.
                    #
                    # Formula: (abs(qend - qstart) + 1) / qlen  — used for ALL search modes.
                    #
                    # AA vs DNA coordinate handling
                    # ─────────────────────────────
                    # BLAST/DIAMOND always reports qstart/qend in the *query's own coordinate
                    # system*, regardless of the search type:
                    #
                    #   blastn          — qstart/qend in nucleotides; qlen in nucleotides
                    #   blastx          — qstart/qend in nucleotides; qlen in nucleotides
                    #                     (query is DNA even though the database is protein)
                    #   blastp / diamond-blastp — qstart/qend in amino acids; qlen in amino acids
                    #
                    # Because the numerator and denominator are always in the same unit, the
                    # span formula is self-consistent for every mode with no conversion needed.
                    # No is_translated_search flag is required here.
                    #
                
                    query_coverage = 0.0
                    try:
                        if qstart is not None and qend is not None and qlen and qlen > 0:
                            query_coverage = (abs(qend - qstart) + 1) / float(qlen) * 100.0
                    except Exception:
                        query_coverage = 0.0

                    if current_read is not None and read_name != current_read:
                        process_query_group(current_read, current_hits)
                        current_hits = []
                    current_read = read_name
                    current_hits.append({
                        'gene': gene,
                        'identity': identity,
                        'query_coverage': query_coverage,
                        'sstart': sstart,
                        'send': send,
                        'bitscore': bitscore,
                        'evalue': evalue,
                    })
                if current_read is not None:
                    process_query_group(current_read, current_hits)
        except KeyError:
            raise
        except Exception as e:
            self.logger.error(f"Error parsing {output_file}: {e}")
            raise RuntimeError(f"Failed to parse {output_file}") from e

        if malformed_lines:
            self.logger.warning(
                f"Skipped {malformed_lines} malformed rows while parsing {output_file}."
            )

        # Finalise statistics and classify exact/family/partial evidence.
        for gene in self.gene_stats[database][tool_name]:
            self.gene_stats[database][tool_name][gene].finalise()
        detected_genes = self._classify_tool_evidence(database, tool_name)

        if self.report_fasta and read_store is None:
            self.logger.info(f"Input reads loaded for FASTA output in {database.upper()} using {tool_name}: {len(all_reads)}")
        # Note: mapped_reads counts hits seen in the BLAST output file
        self.logger.info(f"Reads that returned a hit in {database.upper()} using {tool_name}: {mapped_reads}")
        self.logger.info(f"Reads passing thresholds in {database.upper()} using {tool_name}: {passing_reads}")
        if normalise_detection_system(
                getattr(self, "detection_system", DETECTION_SYSTEM_QUALIFIED)
        ) == DETECTION_SYSTEM_LEGACY_RELAXED:
            self.logger.info(
                f"Detected genes in {database.upper()} using {tool_name}: "
                f"{len(detected_genes)}"
            )
        else:
            evidence_count = sum(
                1 for call in self.evidence_calls.get(database, {}).values()
                if tool_name in call and call[tool_name].evidence_present
            )
            exact_count = sum(
                1 for calls in self.evidence_calls.get(database, {}).values()
                if tool_name in calls and calls[tool_name].exact_allele_detected
            )
            self.logger.info(
                f"Evidence genes in {database.upper()} using {tool_name}: "
                f"{evidence_count}; exact allele genes: {exact_count}"
            )

        # Output FASTA files of reads mapping to genes
        if self.report_fasta and read_store is None:
            self._write_fasta_outputs(database, tool_name, detected_genes, gene_reads, all_reads)

        return detected_genes, gene_reads

    def generate_coverage_from_bam(self, sorted_bam: Path, database: str, tool_name: str) -> None:
        """Generate sparse bedGraph coverage (non-zero positions) from a sorted/indexed BAM.

        Preference order:
        1. Use bedtools genomecov -ibam <bam> -bg to produce bedGraph (only non-zero regions).
        2. Fallback to samtools depth -a and stream-filter non-zero positions, converting to bedGraph.

        Writes file to the sample output directory named '<database>_<tool>_coverage.bedgraph' and
        includes a header row (chrom, start, end, depth) so the visualiser can read it as a tabular file.
        """
        try:
            if not sorted_bam.exists():
                self.logger.debug(f"Coverage generation skipped: BAM not found: {sorted_bam}")
                return
            out_file = self.output_dir / f"{database}_{tool_name.lower()}_coverage.bedgraph"
            out_file.parent.mkdir(parents=True, exist_ok=True)

            # Try bedtools genomecov first
            bedtools_cmd = ['bedtools', 'genomecov', '-ibam', str(sorted_bam), '-bg']
            try:
                self.logger.info(f"Generating bedGraph coverage (bedtools) for {database} using {tool_name}: {out_file}")
                with open(out_file, 'w') as fh:
                    # Write a header so visualiser can parse column names
                    fh.write('chrom\tstart\tend\tdepth\n')
                    proc = subprocess.run(bedtools_cmd, check=True, stdout=fh, stderr=subprocess.PIPE, text=True)
                    if proc.stderr:
                        self.logger.debug(f"bedtools genomecov stderr: {proc.stderr}")
                self.logger.info(f"Wrote bedGraph coverage file: {out_file}")
                return
            except FileNotFoundError:
                self.logger.debug('bedtools not found; falling back to samtools depth')
            except subprocess.CalledProcessError as e:
                self.logger.debug(f'bedtools genomecov failed: {e}; falling back to samtools depth')

            # Fallback: samtools depth -> convert to bedGraph by emitting only depth>0 positions
            self.logger.info(f"Generating bedGraph coverage (samtools fallback) for {database} using {tool_name}: {out_file}")
            with open(out_file, 'w') as fh:
                fh.write('chrom\tstart\tend\tdepth\n')
                proc = subprocess.Popen(['samtools', 'depth', '-a', str(sorted_bam)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                assert proc.stdout is not None
                for line in proc.stdout:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 3:
                        continue
                    chrom = parts[0]
                    try:
                        pos = int(parts[1])
                        depth = float(parts[2])
                    except Exception:
                        continue
                    if depth <= 0:
                        continue
                    # Convert 1-based samtools pos to 0-based bedGraph interval [pos-1, pos)
                    fh.write(f"{chrom}\t{pos-1}\t{pos}\t{depth}\n")
                _, stderr = proc.communicate()
                if stderr:
                    self.logger.debug(f"samtools depth stderr: {stderr}")
            self.logger.info(f"Wrote bedGraph coverage file: {out_file}")
        except Exception as e:
            self.logger.warning(f"Failed to generate coverage for {database} ({tool_name}): {e}")

    def parse_bam_results(self, bam_file: Path, database: str, tool_name: str,
                          count_only: bool = False, read_store=None) -> Set[str]:

        """Parse BAM (samtools view) and collect per-gene stats.

        Uses CIGAR to compute aligned positions and per-read identity; detection
        follows the same coverage-at-depth rules as BLAST parsing.

        When count_only is True, per-gene read-name lists and read sequences are
        replaced with counters. This keeps GeneFior-Recompute bounded on large
        BAM inputs when FASTA reporting is not requested.
        """
        detected_genes = set()
        if not bam_file.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_file}")

        # Ensure structure exists for this database and tool (do not reset entire gene_stats)
        if database not in self.gene_stats:
            self.gene_stats[database] = defaultdict(lambda: defaultdict(GeneStats))
            self.detections.setdefault(database, defaultdict(lambda: defaultdict(bool)))
        elif self.gene_stats[database].get(tool_name) is None:
            self.gene_stats[database][tool_name] = defaultdict(GeneStats)

        gene_lengths = {}  # Store gene lengths from BAM header
        count_only = bool(count_only and (not self.report_fasta or read_store is not None))
        if count_only:
            gene_reads = defaultdict(lambda: {
                'passing_count': 0, 'all_count': 0,
                'passing_r1_count': 0, 'passing_r2_count': 0,
            })
        else:
            gene_reads = defaultdict(lambda: {'passing': [], 'all': [], 'passing_r1': [], 'passing_r2': []})
        all_reads = {}
        mapped_reads = 0
        passing_reads = 0
        stderr = ''


        try:
            # start samtools in its own session; this helps ensure any child threads/processes can be killed together
            proc = subprocess.Popen(['samtools', 'view', '-h', str(bam_file)],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True, close_fds=True)

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
                is_secondary = bool(flag & 0x100)
                is_supplementary = bool(flag & 0x800)
                is_primary = not is_secondary and not is_supplementary
                if not is_primary:
                    continue
                try:
                    mapq = int(fields[4])
                except (TypeError, ValueError):
                    mapq = 0

                # Add /1 or /2 to read name based on SAM flag
                if flag & 0x40:
                    read_name += '/1'
                elif flag & 0x80:
                    read_name += '/2'

                gene = fields[2]
                try:
                    ref_start = int(fields[3]) - 1  # convert to 0-based
                except ValueError:
                    ref_start = 0
                cigar = fields[5]
                seq = fields[9]
                # store read sequence if available
                if self.report_fasta and read_store is None and read_name not in all_reads and seq and seq != '*':
                    all_reads[read_name] = seq

                gene_len = gene_lengths.get(gene, 0)
                if count_only:
                    gene_reads[gene]['all_count'] += 1
                else:
                    gene_reads[gene]['all'].append(read_name)
                if read_store is not None:
                    read_store.add(
                        database, tool_name, gene, read_name,
                        passing=False, sequence=(seq if seq and seq != '*' else None),
                    )
                mapped_reads +=1

                # try to get NM tag from optional fields
                nm = None
                for opt in fields[11:]:
                    if opt.startswith('NM:i:'):
                        try:
                            nm = int(opt.split(':')[-1])
                        except ValueError:
                            nm = None
                        break

                # parse CIGAR and compute aligned positions & alignment length
                ref_pos = ref_start
                aligned_positions = set()  # reference positions — feeds into gene coverage
                aligned_intervals = []  # compact half-open reference spans for recompute
                query_aligned_bases = 0  # query bases aligned — feeds into query coverage
                alignment_length = 0  # for identity calculation

                for count_str, op in cigar_re.findall(cigar):
                    length = int(count_str)
                    if op in ('M', '=', 'X'):
                        aligned_intervals.append((ref_pos, ref_pos + length))
                        if not count_only:
                            aligned_positions.update(range(ref_pos, ref_pos + length))
                        ref_pos += length
                        alignment_length += length
                        # M/=/X consume both query and reference: unambiguously "aligned"
                        query_aligned_bases += length
                    elif op == 'I':
                        # Insertion into reference: consumes query bases but no reference
                        # positions. We count these toward query_aligned_bases so that the
                        # BAM query-coverage formula matches the BLAST span-based approach
                        # (qend - qstart + 1 spans insertion positions in query coords too).
                        # Omitting I would make BAM systematically more conservative than
                        # BLAST/DIAMOND for insertion-bearing alignments.
                        alignment_length += length
                        query_aligned_bases += length
                    elif op == 'D':
                        # Deletion from reference: consumes reference positions but no query
                        # bases — does not advance query_aligned_bases.
                        ref_pos += length
                        alignment_length += length
                    elif op == 'N':
                        ref_pos += length
                    elif op in ('S', 'H'):
                        # Soft/hard clip: not part of the aligned region — excluded.
                        pass

                if alignment_length == 0:
                    identity = None
                elif nm is not None:
                    matches = max(0, alignment_length - nm)
                    identity = (matches / alignment_length) * 100.0
                else:
                    identity = None

                query_length = len(seq) if seq and seq != '*' else 0
                # Per-read query coverage: (M + I bases) / read length, consistent with
                # the BLAST span-based formula (qend - qstart + 1) / qlen.
                query_coverage = (query_aligned_bases / query_length) * 100.0 if query_length > 0 else 0.0

                # Per-read gate uses query_min_identity (same threshold as BLAST/DIAMOND
                # per-HSP filtering) so that all tools are consistent.  detection_min_identity
                # is reserved for the gene-level avg-identity gate applied after finalise().
                identity_passes = (
                    identity is not None
                    and identity >= self.query_min_identity
                )
                if (is_primary
                        and identity_passes
                        and query_coverage >= self.query_min_coverage):
                    if gene not in self.gene_stats[database][tool_name]:
                        self.gene_stats[database][tool_name][gene] = GeneStats(
                            gene_name=gene,
                            compact=count_only,
                        )
                    if count_only:
                        self.gene_stats[database][tool_name][gene].add_intervals(
                            aligned_intervals, identity, gene_len
                        )
                    else:
                        self.gene_stats[database][tool_name][gene].add_positions(
                            aligned_positions, identity, gene_len
                        )
                    if count_only:
                        gene_reads[gene]['passing_count'] += 1
                    else:
                        gene_reads[gene]['passing'].append(read_name)
                    if read_store is not None:
                        read_store.add(
                            database, tool_name, gene, read_name,
                            passing=True, sequence=(seq if seq and seq != '*' else None),
                        )
                    # MAPQ 255 means mapping quality is unavailable, not unique.
                    unique_mapping = 20 <= mapq < 255
                    self.gene_stats[database][tool_name][gene].add_read_support(
                        mapped=True,
                        passing=True,
                        best=True,
                        unique_best=unique_mapping,
                        ambiguous_best=not unique_mapping,
                        high_confidence=unique_mapping,
                        score=mapq,
                    )

                    # Track R1/R2 separately based on read name suffix
                    if read_name.endswith('_R1') or '/1' in read_name:
                        if count_only:
                            gene_reads[gene]['passing_r1_count'] += 1
                        else:
                            gene_reads[gene]['passing_r1'].append(read_name)
                    elif read_name.endswith('_R2') or '/2' in read_name:
                        if count_only:
                            gene_reads[gene]['passing_r2_count'] += 1
                        else:
                            gene_reads[gene]['passing_r2'].append(read_name)

                    passing_reads +=1

            # Close stdout and wait; then ensure process termination including any children
            try:
                if proc and getattr(proc, 'stdout', None):
                    proc.stdout.close()
            except Exception:
                pass
            try:
                if proc:
                    proc.wait()
            except Exception:
                pass

            # optionally capture and log stderr
            try:
                stderr = proc.stderr.read() if proc and proc.stderr else ''
                if stderr:
                    self.logger.debug(f"samtools stderr: {stderr}")
            except Exception:
                pass
            try:
                if proc and getattr(proc, 'stderr', None):
                    proc.stderr.close()
            except Exception:
                pass

            # Best-effort terminate any lingering processes in the group
            try:
                self._terminate_proc(proc, 'samtools')
            except Exception:
                pass

            if proc.returncode not in (0, None):
                raise RuntimeError(
                    f"samtools view failed for {bam_file} with exit code "
                    f"{proc.returncode}: {stderr.strip()}"
                )

        except Exception as e:
            self.logger.error(f"Error reading BAM via samtools: {e}")
            raise RuntimeError(f"Failed to parse BAM file {bam_file}") from e



        # Finalise statistics and classify exact/family/partial evidence.
        for gene in self.gene_stats[database][tool_name]:
            self.gene_stats[database][tool_name][gene].finalise()
        detected_genes = self._classify_tool_evidence(database, tool_name)

        #self.logger.info(f"Total reads processed in {database.upper()} using {tool_name}: {len(all_reads)}")
        ## Need a fix for this
        self.logger.info(f"Reads that returned a hit in {database.upper()} using {tool_name}: {mapped_reads}")
        self.logger.info(f"Reads passing thresholds in {database.upper()} using {tool_name}: {passing_reads}")
        if normalise_detection_system(
                getattr(self, "detection_system", DETECTION_SYSTEM_QUALIFIED)
        ) == DETECTION_SYSTEM_LEGACY_RELAXED:
            self.logger.info(
                f"Detected genes in {database.upper()} using {tool_name}: "
                f"{len(detected_genes)}"
            )
        else:
            evidence_count = sum(
                1 for call in self.evidence_calls.get(database, {}).values()
                if tool_name in call and call[tool_name].evidence_present
            )
            exact_count = sum(
                1 for calls in self.evidence_calls.get(database, {}).values()
                if tool_name in calls and calls[tool_name].exact_allele_detected
            )
            self.logger.info(
                f"Evidence genes in {database.upper()} using {tool_name}: "
                f"{evidence_count}; exact allele genes: {exact_count}"
            )

        # Output FASTA files of reads mapping to genes
        if self.report_fasta and read_store is None:
            self._write_fasta_outputs(database, tool_name, detected_genes, gene_reads, all_reads)


        return detected_genes, gene_reads


    def parse_hmmer_results(self, tbl_file: Path, domtbl_file: Path,
                            database: str, tool_name: str) -> Tuple[Set[str], dict]:
        """Parse HMMER output and extract detected HMM profile names with statistics.

        Uses two complementary output files:

        ``--tblout`` (tbl_file) — one row per (target, query) pair at the full-sequence
        level.  Used only for the final pass/fail E-value decision so we stay consistent
        with hmmsearch's own per-sequence reporting.

            Col 0: target name  (input FASTA sequence)
            Col 2: query name   (HMM profile = the "gene" being detected)
            Col 4: full-seq E-value
            Col 5: full-seq bit-score

        ``--domtblout`` (domtbl_file) — one row per domain hit.  Used for statistics
        because it contains the HMM profile coordinate information needed to compute
        meaningful profile coverage, sequence counts, etc.

            Col 0:  target name       (input sequence)
            Col 2:  tlen              (target sequence length)
            Col 3:  query name        (HMM profile)
            Col 5:  qlen              (profile length in residues)
            Col 6:  full-seq E-value
            Col 7:  full-seq bit-score
            Col 15: hmm_from          (alignment start on profile, 1-based)
            Col 16: hmm_to            (alignment end on profile, 1-based)

        Returns
        -------
        detected_genes : set[str]
            HMM profile names that passed the E-value threshold.
        gene_reads : dict
            gene_reads[profile]['all']     – all input sequences that matched
            gene_reads[profile]['passing'] – sequences that passed the E-value filter
        """
        detected_genes: Set[str] = set()
        gene_reads: dict = defaultdict(lambda: {'passing': set(), 'all': set()})

        # Ensure per-database/tool structures exist
        if database not in self.gene_stats:
            self.gene_stats[database] = defaultdict(lambda: defaultdict(GeneStats))
            self.detections.setdefault(database, defaultdict(lambda: defaultdict(bool)))

        annotations = self.hmmer_annotations.get(database, {})

        # ── Determine whether E-value post-filtering applies ─────────────────────
        # When using per-profile cutoffs (TC/GA/NC), every hit reported by HMMER
        # already passed the cutoff, so we accept all rows unconditionally.
        _tmode = getattr(self, 'hmmer_threshold_mode', 'evalue')
        _use_profile_cutoffs = (_tmode in ('tc', 'ga', 'nc'))

        # Resolve effective HMMER e-value: prefer hmmer_evalue, then global evalue
        _hmmer_ev = getattr(self, 'hmmer_evalue', None)
        if _hmmer_ev is None:
            _hmmer_ev = getattr(self, 'evalue', None)

        # ── Pass 1: tblout – decide which profiles pass the E-value threshold ──────
        # One representative row per (target, query) pair.
        detected_by_evalue: Set[str] = set()
        passing_target_profile_pairs: Set[Tuple[str, str]] = set()
        seq_best_evalue: dict = {}   # profile -> best (lowest) E-value seen
        if tbl_file.exists():
            try:
                with open(tbl_file, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        fields = line.strip().split()
                        if len(fields) < 6:
                            continue
                        target_name = fields[0]
                        gene = fields[2]   # HMM profile name
                        try:
                            evalue = float(fields[4])
                        except ValueError:
                            continue
                        # Track best e-value per profile
                        if gene not in seq_best_evalue or evalue < seq_best_evalue[gene]:
                            seq_best_evalue[gene] = evalue
                        # Accept hit: when using profile cutoffs all tblout rows already passed;
                        # otherwise apply the effective e-value threshold (or accept all if unset).
                        try:
                            if _use_profile_cutoffs or _hmmer_ev is None or evalue <= float(_hmmer_ev):
                                detected_by_evalue.add(gene)
                                passing_target_profile_pairs.add(
                                    (target_name, gene)
                                )
                        except Exception:
                            detected_by_evalue.add(gene)
                            passing_target_profile_pairs.add(
                                (target_name, gene)
                            )
            except Exception as e:
                self.logger.error(f"Error parsing {tbl_file}: {e}")
                raise RuntimeError(f"Failed to parse HMMER table {tbl_file}") from e

        # Fetch per-domain and gene-level thresholds from the workflow instance.
        # These mirror the same thresholds used by parse_blast_results so that
        # HMMER detection is gated consistently.
        #
        #   query_min_coverage  – minimum % of the INPUT SEQUENCE covered by the
        #                         domain alignment  (col ali_to - ali_from / tlen)
        #   query_min_identity  – minimum alignment accuracy (col 21, acc * 100) as
        #                         the closest HMMER equivalent to BLAST %identity
        #   detection_min_coverage – minimum % of the HMM PROFILE covered by the
        #                            union of all domain alignments (gene_coverage)
        _q_min_cov = float(getattr(self, 'query_min_coverage', 0.0))
        _q_min_id  = float(getattr(self, 'query_min_identity',  0.0))
        _d_min_cov = float(getattr(self, 'detection_min_coverage', 0.0))

        # ── Pass 2: domtblout – gather per-domain alignment stats ─────────────────
        # Multiple domain rows can exist per (target, query) pair; each row gives
        # per-domain alignment coordinates used for profile-coverage statistics.
        #
        # domtblout column layout (0-based):
        #   0  target_name   – input sequence name
        #   2  tlen          – input sequence length  (residues)
        #   3  query_name    – HMM profile name
        #   5  qlen          – profile length (residues)
        #   6  E-value       – full-sequence E-value
        #   7  score         – full-sequence bit-score
        #  15  hmm_from      – domain alignment start on the profile (1-based)
        #  16  hmm_to        – domain alignment end   on the profile (1-based)
        #  17  ali_from      – domain alignment start on the target  (1-based)
        #  18  ali_to        – domain alignment end   on the target  (1-based)
        #  21  acc           – mean posterior probability of aligned residues (0–1)
        if domtbl_file and domtbl_file.exists():
            try:
                with open(domtbl_file, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        fields = line.strip().split()
                        # Need at least 22 columns (index 0–21)
                        if len(fields) < 22:
                            continue
                        target_name = fields[0]
                        gene        = fields[3]   # HMM profile name
                        try:
                            tlen     = int(fields[2])    # input sequence length
                            qlen     = int(fields[5])    # profile length (residues)
                            score    = float(fields[7])  # full-seq bit-score
                            hmm_from = int(fields[15])   # profile alignment start
                            hmm_to   = int(fields[16])   # profile alignment end
                            ali_from = int(fields[17])   # target alignment start
                            ali_to   = int(fields[18])   # target alignment end
                            acc      = float(fields[21]) # alignment accuracy (0–1)
                        except (ValueError, IndexError):
                            continue

                        # ── Per-domain query-level filters ──────────────────────
                        # 1. Query-sequence coverage: how much of the input sequence
                        #    is spanned by this domain (analogous to HSP query cov).
                        query_seq_cov = (abs(ali_to - ali_from) + 1) / tlen * 100 if tlen > 0 else 0.0
                        passes_q_cov = (_q_min_cov <= 0.0 or query_seq_cov >= _q_min_cov)

                        # 2. Alignment accuracy as a proxy for %identity.
                        #    acc ranges 0–1; compare against query_min_identity/100.
                        passes_q_id = (_q_min_id <= 0.0 or (acc * 100.0) >= _q_min_id)

                        domain_passes = passes_q_cov and passes_q_id

                        # Always count in gene_reads['all'] for any hit (mirrors
                        # how BLAST 'all' includes any sequence that maps at all).
                        gene_reads[gene]['all'].add(target_name)

                        pair_passes = (
                            target_name, gene
                        ) in passing_target_profile_pairs
                        if not domain_passes or not pair_passes:
                            # Domain is too poor quality to contribute to stats
                            continue

                        # Initialise GeneStats on first qualifying hit
                        if gene not in self.gene_stats[database][tool_name]:
                            self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)

                        # Use HMM profile coordinates + bit-score for coverage/identity stats
                        self.gene_stats[database][tool_name][gene].add_hit(
                            hmm_from, hmm_to, score, qlen
                        )

                        # 'passing' = domain passed quality filters AND the full-sequence
                        # E-value threshold from Pass 1
                        gene_reads[gene]['passing'].add(target_name)

            except Exception as e:
                self.logger.error(f"Error parsing {domtbl_file}: {e}")
                raise RuntimeError(
                    f"Failed to parse HMMER domain table {domtbl_file}"
                ) from e

        else:
            # domtblout not available – fall back to tblout-only (no position data)
            self.logger.debug("  domtblout not found; using tblout only (no profile coverage or per-domain quality filtering)")
            if tbl_file.exists():
                try:
                    with open(tbl_file, 'r') as f:
                        for line in f:
                            if line.startswith('#'):
                                continue
                            fields = line.strip().split()
                            if len(fields) < 6:
                                continue
                            target_name = fields[0]
                            gene        = fields[2]
                            try:
                                score = float(fields[5])
                            except ValueError:
                                score = 0.0
                            if (target_name, gene) not in passing_target_profile_pairs:
                                gene_reads[gene]['all'].add(target_name)
                                continue
                            if gene not in self.gene_stats[database][tool_name]:
                                self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)
                            self.gene_stats[database][tool_name][gene].add_hit(1, 1, score)
                            gene_reads[gene]['all'].add(target_name)
                            gene_reads[gene]['passing'].add(target_name)
                except Exception as e:
                    self.logger.error(f"Error in tblout fallback: {e}")
                    raise RuntimeError(
                        f"Failed to parse HMMER table {tbl_file}"
                    ) from e

        # Finalise profile statistics and add unique target-sequence support.
        for gene, stats in self.gene_stats[database][tool_name].items():
            stats.mapped_read_support = len(gene_reads[gene]['all'])
            stats.passing_read_support = len(gene_reads[gene]['passing'])
            stats.best_read_support = stats.passing_read_support
            stats.unique_best_read_support = stats.passing_read_support
            stats.high_confidence_read_support = stats.passing_read_support
            stats.finalise()

        must_flag_overrides = {
            gene
            for gene in detected_by_evalue
            if getattr(self, 'hmmer_must_flag', True)
            and annotations.get(gene, {}).get('must_flag')
        }
        detected_genes = self._classify_tool_evidence(
            database,
            tool_name,
            must_flag_overrides=must_flag_overrides,
            eligible_genes=detected_by_evalue,
        )

        # Convert sets → lists for write_tool_stats
        gene_reads_lists = {
            gene: {
                'all':     list(v['all']),
                'passing': list(v['passing']),
            }
            for gene, v in gene_reads.items()
        }
        return detected_genes, gene_reads_lists

    def write_tool_stats(self, database: str, tool_name: str, gene_reads: dict = None):
        """Write detailed statistics for a specific tool to TSV.

        Output columns:
        - Gene: Gene name
        - Gene_Length: Length of the gene in the database (bp)
        - Num_Sequences_Mapped: Number of sequences that mapped to this gene with identity >= detection-min-identity
        - Num_Sequences_Passing_Thresholds: Number of sequences that passed all thresholds (identity, coverage, base_coverage, min_reads etc)
        - Gene_Coverage: Percentage of the gene covered by all qualifying alignments combined (%)
        - Base_Coverage: Average depth across the entire gene
        - Detection_Depth: Per-base depth required across the configured detection coverage
        - Detection_Depth_Coverage: Percentage of the gene covered at Detection_Depth or higher
        - Avg_Identity: Average identity across all qualifying sequences (%)
        - Detected: result selected by Detection_System
        """
        # Normalise gene_reads to a safe empty structure if caller did not provide one
        if gene_reads is None:
            gene_reads = defaultdict(lambda: {'passing': [], 'all': []})

        def _read_count(read_info: dict, key: str) -> int:
            """Return a read/hit count from either list-backed or count-only data."""
            if not isinstance(read_info, dict):
                return 0
            count_key = f"{key}_count"
            if count_key in read_info:
                try:
                    return int(read_info.get(count_key, 0))
                except Exception:
                    return 0
            value = read_info.get(key, [])
            if isinstance(value, int):
                return value
            try:
                return len(value)
            except Exception:
                return 0

        stats_file = self.stats_dir / f"{database}_{tool_name}_stats.tsv"
        detection_system = normalise_detection_system(
            getattr(self, "detection_system", DETECTION_SYSTEM_QUALIFIED)
        )

        gene_stats = self.gene_stats[database][tool_name]
        if not gene_stats:
            self.logger.warning(f"No statistics to write for {database} - {tool_name}")
            return

        with open(stats_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            metric_header = [
                'Gene', 'Gene_Length', 'Num_Sequences_Mapped',
                'Num_Sequences_Passing_Thresholds', 'Gene_Coverage',
                'Coverage_2x', 'Coverage_3x', 'Coverage_5x', 'Coverage_10x',
                'Base_Coverage', 'Detection_Depth',
                'Detection_Depth_Coverage', 'Median_Depth',
                'Depth_CV', 'Num_Internal_Gaps', 'Longest_Internal_Gap',
                'Num_Gaps', 'Longest_Gap', 'Longest_Covered_Run',
                'Avg_Identity',
            ]
            if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED:
                header = metric_header + ['Detected']
            else:
                header = metric_header + [
                'Mapped_Read_Support', 'Passing_Read_Support',
                'Best_Read_Support', 'Unique_Best_Read_Support',
                'Ambiguous_Best_Read_Support', 'Evidence_Status',
                'Evidence_Warnings', 'Family', 'Top_Database_Candidate',
                'Competing_Alleles', 'Always_Flagged', 'Evidence_Present',
                'Candidate_Allele_Detected', 'Exact_Allele_Detected',
                'Exact_Protein_Detected', 'Profile_Detected',
                'Detection_System', 'Detected',
                ]
            writer.writerow(header)

            # Sort genes alphabetically
            genes = sorted(gene_stats.keys())

            for gene in genes:
                stats = gene_stats[gene]
                try:
                    detected = self.detections[database][gene][tool_name]
                except (KeyError, TypeError):
                    detected = False
                if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED:
                    detection_depth = (
                        self.detection_min_depth
                        if not getattr(self, 'is_genes_fasta', False)
                        else 1
                    )
                else:
                    detection_depth = (
                        self.evidence_config.corroborating_depth
                        if not getattr(self, 'is_genes_fasta', False)
                        else 1
                    )
                detection_depth_coverage = stats.coverage_at_depth(
                    detection_depth
                )
                metric_row = [
                    gene,
                    stats.gene_length,
                    _read_count(gene_reads.get(gene, {}), 'all'), # 'all' reads mapping to gene
                    _read_count(gene_reads.get(gene, {}), 'passing'), # Just those that 'passed' thresholds
                    f"{stats.gene_coverage:.2f}",
                    f"{stats.coverage_2x:.2f}",
                    f"{stats.coverage_3x:.2f}",
                    f"{stats.coverage_5x:.2f}",
                    f"{stats.coverage_10x:.2f}",
                    f"{stats.base_depth:.2f}",
                    detection_depth,
                    f"{detection_depth_coverage:.2f}",
                    f"{stats.median_depth:.2f}",
                    f"{stats.depth_cv:.4f}",
                    stats.num_internal_gaps,
                    stats.longest_internal_gap,
                    stats.num_gaps,
                    stats.longest_gap,
                    stats.longest_covered_run,
                    f"{stats.avg_identity:.2f}",
                ]
                if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED:
                    row = metric_row + ['1' if detected else '0']
                    writer.writerow(row)
                    continue

                call = (
                    self.evidence_calls
                    .get(database, {})
                    .get(gene, {})
                    .get(tool_name)
                )
                row = metric_row + [
                    stats.mapped_read_support,
                    stats.passing_read_support,
                    stats.best_read_support,
                    stats.unique_best_read_support,
                    stats.ambiguous_best_read_support,
                    call.status if call else NOT_DETECTED,
                    call.warning_text() if call else '',
                    call.family if call else '',
                    call.best_allele if call else gene,
                    ';'.join(call.competing_alleles) if call else '',
                    '1' if gene in getattr(self, 'always_flag_genes', set()) else '0',
                    '1' if call and call.evidence_present else '0',
                    '1' if call and call.candidate_allele_detected else '0',
                    '1' if call and call.exact_allele_detected else '0',
                    '1' if call and call.status == EXACT_PROTEIN_DETECTED else '0',
                    '1' if call and call.profile_detected else '0',
                    detection_system,
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
        detection_system = normalise_detection_system(
            getattr(self, "detection_system", DETECTION_SYSTEM_QUALIFIED)
        )
        if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED:
            self._remove_qualified_reports(database)
        else:
            self.generate_evidence_matrix(database)
        if getattr(self, 'db_whole_genome', False):
            self.generate_genome_mapping_summary(database)

    def _remove_qualified_reports(self, database: str):
        """Remove stale qualified-only outputs from a legacy result directory."""
        for suffix in (
                'evidence_matrix.tsv',
                'evidence_summary.tsv',
                'allele_resolution.tsv'):
            path = self.output_dir / f"{database}_{suffix}"
            try:
                path.unlink()
                self.logger.info(f"Removed stale qualified report: {path}")
            except FileNotFoundError:
                pass

    def generate_genome_mapping_summary(self, database: str):
        """Write contig-level detail and aggregated multi-contig genome calls."""
        columns = [
            'Genome_ID', 'Contig_ID', 'Contig_Count', 'Reference_Length',
            'Tool', 'Coverage_1x', 'Coverage_2x', 'Coverage_3x',
            'Coverage_5x', 'Coverage_10x', 'Mean_Depth', 'Median_Depth',
            'Depth_CV', 'Mean_Identity', 'Mapped_Read_Support',
            'Passing_Read_Support', 'Best_Read_Support',
            'Unique_Best_Read_Support', 'Ambiguous_Best_Read_Support',
            'Num_Gaps', 'Longest_Gap', 'Num_Internal_Gaps',
            'Longest_Internal_Gap', 'Detection_Depth',
            'Detection_Depth_Coverage', 'Detection_Min_Coverage',
            'Detection_Min_Identity', 'Evidence_Status', 'Evidence_Warnings',
        ]
        database_stats = self.gene_stats.get(database, {})
        contig_file = self.output_dir / f'{database}_contig_mapping_summary.tsv'
        genome_file = self.output_dir / f'{database}_genome_mapping_summary.tsv'
        grouped = defaultdict(list)
        detection_depth = (
            self.evidence_config.corroborating_depth
            if not getattr(self, 'is_genes_fasta', False)
            else 1
        )

        with open(contig_file, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(columns)
            for tool_name in sorted(database_stats):
                for reference in sorted(database_stats[tool_name]):
                    stats = database_stats[tool_name][reference]
                    genome_id, contig_id = split_whole_genome_reference(reference)
                    grouped[(tool_name, genome_id)].append((contig_id, stats))
                    call = self.evidence_calls.get(database, {}).get(reference, {}).get(tool_name)
                    writer.writerow(self._mapping_summary_row(
                        genome_id, contig_id, 1, tool_name, stats,
                        detection_depth, call.status if call else NOT_DETECTED,
                        call.warning_text() if call else '',
                        self.detection_min_coverage, self.detection_min_identity,
                    ))

        with open(genome_file, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(columns + ['Aggregation_Level', 'Contig_List'])
            for (tool_name, genome_id), records in sorted(grouped.items()):
                metrics = self._aggregate_genome_mapping(records, detection_depth)
                writer.writerow(self._mapping_summary_row(
                    genome_id, ','.join(metrics.contigs), metrics.contig_count,
                    tool_name, metrics, detection_depth, metrics.status,
                    ';'.join(metrics.warnings), self.detection_min_coverage,
                    self.detection_min_identity,
                ) + ['GENOME', ','.join(metrics.contigs)])
        self.logger.info(f'Generated contig mapping summary: {contig_file}')
        self.logger.info(f'Generated whole-genome mapping summary: {genome_file}')

    @staticmethod
    def _mapping_summary_row(genome_id, contig_id, contig_count, tool_name,
                             stats, detection_depth, status, warnings,
                             detection_min_coverage, detection_min_identity):
        coverage = lambda depth: stats.coverage_at_depth(depth)
        return [
            genome_id, contig_id, contig_count, stats.gene_length, tool_name,
            f'{coverage(1):.2f}', f'{coverage(2):.2f}', f'{coverage(3):.2f}',
            f'{coverage(5):.2f}', f'{coverage(10):.2f}',
            f'{stats.base_depth:.2f}', f'{stats.median_depth:.2f}',
            f'{stats.depth_cv:.4f}', f'{stats.avg_identity:.2f}',
            stats.mapped_read_support, stats.passing_read_support,
            stats.best_read_support, stats.unique_best_read_support,
            stats.ambiguous_best_read_support, stats.num_gaps,
            stats.longest_gap, stats.num_internal_gaps,
            stats.longest_internal_gap, detection_depth,
            f'{coverage(detection_depth):.2f}',
            f'{detection_min_coverage:.2f}', f'{detection_min_identity:.2f}',
            status, warnings,
        ]

    def _aggregate_genome_mapping(self, records, detection_depth):
        """Aggregate contigs using reference-length-weighted base metrics."""
        total_length = sum(max(0, s.gene_length) for _, s in records)
        depth_hist = defaultdict(int)
        total_depth = total_depth_sq = 0
        identity_total = sequence_total = 0
        totals = defaultdict(int)
        longest_gap = longest_internal_gap = 0
        num_gaps = num_internal_gaps = 0
        discontinuous = False
        contigs = [name for name, _ in records]
        for _, stats in records:
            for start, end, depth in stats._iter_depth_segments():
                span = max(0, end - start)
                depth_hist[depth] += span
                total_depth += depth * span
                total_depth_sq += depth * depth * span
            identity_total += stats.identity_sum or (stats.avg_identity * stats.num_sequences)
            sequence_total += stats.num_sequences
            for field in ('mapped_read_support', 'passing_read_support',
                          'best_read_support', 'unique_best_read_support',
                          'ambiguous_best_read_support'):
                totals[field] += getattr(stats, field, 0)
            num_gaps += stats.num_gaps
            num_internal_gaps += stats.num_internal_gaps
            longest_gap = max(longest_gap, stats.longest_gap)
            longest_internal_gap = max(longest_internal_gap, stats.longest_internal_gap)
            max_gap = max(
                self.evidence_config.max_internal_gap_bp,
                int(stats.gene_length * self.evidence_config.max_internal_gap_fraction),
            )
            discontinuous |= (
                stats.num_internal_gaps > 1 or stats.longest_internal_gap > max_gap
            )

        if total_length:
            coverage = {
                depth: sum(span for level, span in depth_hist.items() if level >= depth)
                / total_length * 100
                for depth in (1, 2, 3, 5, 10)
            }
            mean_depth = total_depth / total_length
            variance = max(0.0, total_depth_sq / total_length - mean_depth ** 2)
            depth_cv = (variance ** 0.5 / mean_depth) if mean_depth else 0.0
            midpoint = (total_length + 1) / 2
            running = 0
            median_depth = 0.0
            for depth in sorted(depth_hist):
                running += depth_hist[depth]
                if running >= midpoint:
                    median_depth = float(depth)
                    break
        else:
            coverage = {depth: 0.0 for depth in (1, 2, 3, 5, 10)}
            mean_depth = median_depth = depth_cv = 0.0

        mean_identity = identity_total / sequence_total if sequence_total else 0.0
        detection_coverage = coverage.get(detection_depth, 0.0)
        passes = (
            detection_coverage >= self.detection_min_coverage
            and mean_identity >= self.detection_min_identity
            and totals['passing_read_support'] >= self.detection_min_num_reads
        )
        partial_signal = coverage[1] >= self.evidence_config.partial_min_coverage
        complete = (
            passes
            and coverage[1] >= 100.0 - 1e-9
            and detection_coverage >= 95.0
            and not discontinuous
        )
        near_complete = passes and detection_coverage >= 80.0
        if complete:
            status = WHOLE_GENOME_COMPLETE
        elif near_complete:
            status = WHOLE_GENOME_NEAR_COMPLETE
        elif partial_signal:
            status = WHOLE_GENOME_PARTIAL
        elif passes:
            status = MIXED_OR_MOSAIC
        else:
            status = NOT_DETECTED
        warnings = []
        if detection_coverage < self.detection_min_coverage:
            warnings.append('LOW_CORROBORATED_COVERAGE')
        if mean_identity < self.detection_min_identity:
            warnings.append('LOW_MEAN_IDENTITY')
        if totals['passing_read_support'] < self.detection_min_num_reads:
            warnings.append('LOW_READ_SUPPORT')
        if discontinuous:
            warnings.append('DISCONTINUOUS_COVERAGE')
        result = GeneStats(gene_name=','.join(contigs), gene_length=total_length)
        result._coverage_by_depth = coverage
        result.gene_coverage = coverage[1]
        result.coverage_2x = coverage[2]
        result.coverage_3x = coverage[3]
        result.coverage_5x = coverage[5]
        result.coverage_10x = coverage[10]
        result.base_depth = mean_depth
        result.median_depth = median_depth
        result.depth_cv = depth_cv
        result.avg_identity = mean_identity
        result.identity_sum = identity_total
        result.num_sequences = sequence_total
        result.num_gaps = num_gaps
        result.longest_gap = longest_gap
        result.num_internal_gaps = num_internal_gaps
        result.longest_internal_gap = longest_internal_gap
        for field, value in totals.items():
            setattr(result, field, value)
        result.detection_min_coverage = self.detection_min_coverage
        result.detection_min_identity = self.detection_min_identity
        result.status = status
        result.warnings = warnings
        result.contigs = contigs
        result.contig_count = len(contigs)
        return result

    def generate_evidence_matrix(self, database: str):
        """Write qualified per-tool calls and family-level allele resolution."""
        detection_system = normalise_detection_system(
            getattr(self, "detection_system", DETECTION_SYSTEM_QUALIFIED)
        )
        if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED:
            self._remove_qualified_reports(database)
            return
        database_calls = self.evidence_calls.get(database, {})
        all_tools = sorted({
            tool
            for gene_calls in database_calls.values()
            for tool in gene_calls
        })
        if not all_tools:
            return

        evidence_file = self.output_dir / f"{database}_evidence_matrix.tsv"
        genes = sorted(database_calls)
        with open(evidence_file, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow([
                'Gene',
                'Detection_System',
                *all_tools,
                'Best_Evidence_Status',
                'Evidence_Detections',
                'Candidate_Allele_Detections',
                'Exact_Allele_Detections',
                'Exact_Protein_Detections',
                'Profile_Detections',
                'Strict_Detections',
                'Always_Flagged',
                'Evidence_Warnings',
            ])
            for gene in genes:
                calls = database_calls[gene]
                statuses = [
                    calls[tool].status if tool in calls else NOT_DETECTED
                    for tool in all_tools
                ]
                exact_count = sum(
                    1 for tool in all_tools
                    if tool in calls and calls[tool].exact_allele_detected
                )
                exact_protein_count = sum(
                    1 for tool in all_tools
                    if tool in calls and calls[tool].status == EXACT_PROTEIN_DETECTED
                )
                candidate_count = sum(
                    1 for tool in all_tools
                    if tool in calls and calls[tool].candidate_allele_detected
                )
                profile_count = sum(
                    1 for tool in all_tools
                    if tool in calls and calls[tool].profile_detected
                )
                strict_count = sum(
                    1 for tool in all_tools
                    if tool in calls and calls[tool].exact_detected
                )
                evidence_count = sum(
                    1 for tool in all_tools
                    if tool in calls and calls[tool].evidence_present
                )
                warnings = sorted({
                    warning
                    for call in calls.values()
                    for warning in call.warnings
                })
                writer.writerow([
                    gene,
                    detection_system,
                    *statuses,
                    best_status(statuses),
                    evidence_count,
                    candidate_count,
                    exact_count,
                    exact_protein_count,
                    profile_count,
                    strict_count,
                    '1' if gene in getattr(self, 'always_flag_genes', set()) else '0',
                    ';'.join(warnings),
                ])
        self.logger.info(f"Generated detection evidence matrix: {evidence_file}")

        summary_file = self.output_dir / f"{database}_evidence_summary.tsv"
        with open(summary_file, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow([
                'Tool', 'Detection_System', 'Evidence_Genes',
                'Candidate_Allele_Genes', 'Exact_Allele_Genes',
                'Exact_Protein_Genes', 'Profile_Genes',
                'Exact_Or_Profile_Genes',
                'Always_Flagged_Genes',
            ])

            def _summary_counts(tool=None):
                selected = [
                    call
                    for gene_calls in database_calls.values()
                    for call_tool, call in gene_calls.items()
                    if tool is None or call_tool == tool
                ]
                evidence_genes = {
                    call.gene for call in selected if call.evidence_present
                }
                exact_genes = {
                    call.gene for call in selected
                    if call.exact_allele_detected
                }
                exact_protein_genes = {
                    call.gene for call in selected
                    if call.status == EXACT_PROTEIN_DETECTED
                }
                candidate_genes = {
                    call.gene for call in selected
                    if call.candidate_allele_detected
                }
                profile_genes = {
                    call.gene for call in selected if call.profile_detected
                }
                strict_genes = {
                    call.gene for call in selected if call.exact_detected
                }
                always_flagged_genes = {
                    call.gene for call in selected
                    if call.gene in getattr(self, 'always_flag_genes', set())
                }
                return (
                    len(evidence_genes),
                    len(candidate_genes),
                    len(exact_genes),
                    len(exact_protein_genes),
                    len(profile_genes),
                    len(strict_genes),
                    len(always_flagged_genes),
                )

            for tool in all_tools:
                writer.writerow([
                    tool, detection_system, *_summary_counts(tool)
                ])
            writer.writerow([
                'ALL_TOOLS', detection_system, *_summary_counts()
            ])
        self.logger.info(f"Generated evidence count summary: {summary_file}")

        family_file = self.output_dir / f"{database}_allele_resolution.tsv"
        with open(family_file, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow([
                'Family', 'Tool', 'Status', 'Top_Database_Candidate',
                'Candidate_Allele_Resolved', 'Exact_Allele_Resolved',
                'Exact_Protein_Resolved',
                'Competing_Alleles', 'Warnings',
            ])
            for tool in sorted(self.family_calls.get(database, {})):
                for family, summary in sorted(
                        self.family_calls[database][tool].items()):
                    writer.writerow([
                        family,
                        tool,
                        summary['status'],
                        summary['best_allele'],
                        '1' if summary['candidate_allele_resolved'] else '0',
                        '1' if summary['exact_allele_resolved'] else '0',
                        '1' if summary['exact_protein_resolved'] else '0',
                        ';'.join(summary['competing_alleles']),
                        ';'.join(summary['warnings']),
                    ])
        self.logger.info(f"Generated allele-resolution report: {family_file}")


    def run_workflow(self,options):
        results = {}

        if self.db_whole_genome:
            # Whole-genome references are nucleotide mapping targets, not
            # gene/allele databases. Keep the run mapping-oriented even when
            # Workflow is called directly rather than through the CLI.
            requested_tools = set(self.tools or [])
            allowed_tools = {'minimap2', 'bowtie2', 'bwa'}
            if 'blastn' in requested_tools or 'all' in requested_tools:
                allowed_tools.add('blastn')
            self.tools = sorted(
                tool for tool in requested_tools if tool in allowed_tools
            ) or sorted({'minimap2', 'bowtie2', 'bwa'})
            self.run_dna = True
            self.run_protein = False

        # If user indicated the inputs are full-length gene FASTA(s), adjust internal
        # flags: force Single-FASTA sequence type and (optionally) restrict which
        # tool classes to run based on the declared genes_type.
        # If the declared sequence type is Genes-FASTA, treat inputs as
        # full-length gene FASTA(s): force internal sequence_type to
        # 'Single-FASTA' for downstream logic and attempt auto-detection of
        # DNA vs protein to influence which tools are run.
        if getattr(self, 'sequence_type', None) == 'Genes-FASTA':
            try:
                self.logger.info("Genes-input mode enabled: treating input as full-length gene FASTA(s). Skipping read-mappers.")
                # Preserve original sequence_type but set a dedicated flag so
                # downstream checks can decide behavior without relying on a
                # mutated `sequence_type`.
                self.is_genes_fasta = True
                # If the caller provided an explicit genes_type use it; otherwise auto-detect
                if getattr(self, 'genes_type', None):
                    detected = getattr(self, 'genes_type')
                    self.logger.info(f"Using user-specified genes_type: {detected}")
                else:
                    detected = self._detect_genes_type()
                if detected == 'protein':
                    self.logger.info("Detected gene FASTA type: protein. Will run protein searches (BLASTp/DIAMOND blastp).")
                    self.run_dna = False
                    self.run_protein = True
                elif detected == 'dna':
                    self.logger.info("Detected gene FASTA type: DNA. Will run BLASTn and translated protein searches (BLASTx/DIAMOND blastx).")
                    self.run_dna = True
                    self.run_protein = True
                else:
                    self.logger.info("Could not auto-detect gene FASTA type; leaving run_dna/run_protein as configured.")
                self.detected_genes_type = detected
            except Exception:
                # Non-fatal: if detection fails we continue with configured flags
                self.detected_genes_type = None

        # Iterate over each database in the provided `databases` dictionary
        for db_name, db_paths in self.databases.items():
            results[db_name] = {}
            self.logger.info(f"\n### Processing {db_name.capitalize()} Database ###")

            # BLASTn for DNA queries
            if self.run_dna and self._tool_enabled('blastn') and db_paths.get('blastn'):
                results[db_name]['BLASTn-DNA'] = self.run_blast(
                    db_paths['blastn'], db_name, 'blastn')

            # Protein searches: choose BLASTp for protein queries, BLASTx for nucleotide->protein
            if self.run_protein:
                # If inputs are Genes-FASTA and detected/declared as protein prefer
                # BLASTp when available; otherwise fall back to BLASTx.
                if getattr(self, 'is_genes_fasta', False) and getattr(self, 'detected_genes_type', None) == 'protein':
                    # prefer blastp if available
                    if self._tool_enabled('blastp') and db_paths.get('blastp'):
                        results[db_name]['BLASTp-AA'] = self.run_blast(db_paths['blastp'], db_name, 'blastp')
                    elif self._tool_enabled('blastx') and db_paths.get('blastx'):
                        self.logger.warning(f"No BLASTp available for {db_name}; falling back to BLASTx (may be suboptimal for protein queries)")
                        results[db_name]['BLASTx-AA'] = self.run_blast(db_paths['blastx'], db_name, 'blastx')
                else:
                    if self._tool_enabled('blastx') and db_paths.get('blastx'):
                        results[db_name]['BLASTx-AA'] = self.run_blast(db_paths['blastx'], db_name, 'blastx')

            # DIAMOND: choose blastp for protein queries, blastx for nucleotide queries
            if self.run_protein and self._tool_enabled('diamond') and db_paths.get('diamond'):
                if getattr(self, 'is_genes_fasta', False) and getattr(self, 'detected_genes_type', None) == 'protein':
                    results[db_name]['DIAMOND-AA'] = self.run_diamond(db_paths['diamond'], db_name, query_mode='blastp')
                else:
                    results[db_name]['DIAMOND-AA'] = self.run_diamond(db_paths['diamond'], db_name, query_mode='blastx')

            # Skip read-mappers entirely when the user has indicated the input
            # are full-length gene FASTA(s). In that case we only run BLAST/DIAMOND/HMMER
            # style searches appropriate for full-length sequences.
            if (not getattr(self, 'is_genes_fasta', False)) and self.run_dna and self._tool_enabled('bowtie2') and db_paths.get('bowtie2'):
                results[db_name]['Bowtie2-DNA'] = self.run_bowtie2(
                    db_paths['bowtie2'], db_name)

            if (not getattr(self, 'is_genes_fasta', False)) and self.run_dna and self._tool_enabled('bwa') and db_paths.get('bwa'):
                results[db_name]['BWA-DNA'] = self.run_bwa(
                    db_paths['bwa'], db_name)

            if (not getattr(self, 'is_genes_fasta', False)) and self._tool_enabled('minimap2') and db_paths.get('minimap2'):
                results[db_name]['Minimap2-DNA'] = self.run_minimap2(
                    db_paths['minimap2'], db_name, options.minimap2_preset)

            # HMMER protein search (hmmsearch) – runs against protein HMM profiles.
            # Requires protein sequences as input; suited for biorisk detection and
            # any other protein-family HMM database.
            if self.run_protein and self._tool_enabled('hmmer_protein') and db_paths.get('hmmer_protein'):
                results[db_name]['HMMER-PROTEIN'] = self.run_hmmer(
                    db_paths['hmmer_protein'], db_name, 'protein')

            # HMMER DNA search (nhmmer) – runs against nucleotide HMM profiles.
            if self.run_dna and self._tool_enabled('hmmer_dna') and db_paths.get('hmmer_dna'):
                results[db_name]['HMMER-DNA'] = self.run_hmmer(
                    db_paths['hmmer_dna'], db_name, 'dna')

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

   

        # Final summary
        self.logger.info("\n" + "=" * 70)
        self.logger.info("PIPELINE SUMMARY")
        self.logger.info("=" * 70)
        detection_system = normalise_detection_system(
            getattr(self, "detection_system", DETECTION_SYSTEM_QUALIFIED)
        )

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
                    evidence_tool = {
                        'BLASTn-DNA': 'BLASTn',
                        'BLASTx-AA': 'BLASTx',
                        'BLASTp-AA': 'BLASTp',
                        'Bowtie2-DNA': 'Bowtie2',
                        'BWA-DNA': 'BWA',
                        'Minimap2-DNA': 'Minimap2',
                    }.get(tool, tool)
                    if tool == 'DIAMOND-AA':
                        evidence_tool = (
                            'DIAMOND-BLASTP'
                            if (getattr(self, 'is_genes_fasta', False)
                                and getattr(
                                    self, 'detected_genes_type', None
                                ) == 'protein')
                            else 'DIAMOND-BLASTX'
                        )
                    if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED:
                        detected_count = len(genes) if isinstance(
                            genes, (set, list, tuple, dict)
                        ) else 0
                        self.logger.info(
                            f"  {status} {evidence_tool:.<30} "
                            f"{detected_count} detected"
                        )
                        continue
                    tool_calls = [
                        gene_calls[evidence_tool]
                        for gene_calls in self.evidence_calls
                        .get(db_name, {}).values()
                        if evidence_tool in gene_calls
                    ]
                    evidence_count = sum(
                        call.evidence_present for call in tool_calls
                    )
                    exact_count = sum(
                        call.exact_allele_detected for call in tool_calls
                    )
                    candidate_count = sum(
                        call.candidate_allele_detected for call in tool_calls
                    )
                    profile_count = sum(
                        call.profile_detected for call in tool_calls
                    )
                    # Append MUST-FLAG count for HMMER tools when annotations are present
                    must_flag_suffix = ""
                    if 'HMMER' in tool and isinstance(genes, set):
                        try:
                            _ann = self.hmmer_annotations.get(db_name, {})
                            _meta = self.hmmer_annotations_meta.get(db_name, {})
                            if _meta.get('has_must_flag'):
                                _mf = sum(1 for g in genes if _ann.get(g, {}).get('must_flag'))
                                if _mf > 0:
                                    must_flag_suffix = f"  ⚑  {_mf} must-flag"
                        except Exception:
                            pass
                    self.logger.info(
                        f"  {status} {evidence_tool:.<30} "
                        f"{evidence_count} evidence; "
                        f"{candidate_count} candidate alleles; "
                        f"{exact_count} exact alleles; "
                        f"{profile_count} profiles"
                        f"{must_flag_suffix}"
                    )
                if detection_system == DETECTION_SYSTEM_LEGACY_RELAXED:
                    detected_genes = {
                        gene
                        for gene, tool_calls in self.detections
                        .get(db_name, {}).items()
                        if any(tool_calls.values())
                    }
                    self.logger.info(
                        f"  DETECTED: {len(detected_genes)} genes"
                    )
                    continue

                evidence_calls = [
                    call
                    for gene_calls in self.evidence_calls.get(db_name, {}).values()
                    for call in gene_calls.values()
                    if call.evidence_present
                ]
                evidence_genes = {
                    call.gene for call in evidence_calls
                }
                exact_allele_calls = sum(
                    call.exact_allele_detected for call in evidence_calls
                )
                exact_allele_genes = {
                    call.gene for call in evidence_calls
                    if call.exact_allele_detected
                }
                candidate_allele_calls = sum(
                    call.candidate_allele_detected for call in evidence_calls
                )
                candidate_allele_genes = {
                    call.gene for call in evidence_calls
                    if call.candidate_allele_detected
                }
                profile_calls = sum(
                    call.profile_detected for call in evidence_calls
                )
                if evidence_calls:
                    self.logger.info(
                        f"  Evidence: {len(evidence_genes)} genes "
                        f"({len(evidence_calls)} gene/tool calls); "
                        f"candidate alleles: {len(candidate_allele_genes)} genes "
                        f"({candidate_allele_calls} calls); "
                        f"exact alleles: {len(exact_allele_genes)} genes "
                        f"({exact_allele_calls} calls); "
                        f"profiles: {profile_calls} calls"
                    )
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
