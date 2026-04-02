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
from typing import Any, Dict, List, Set, Tuple, Iterator
import re
from typing import Optional


try:
    from .gene_stats import GeneStats
    from .constants import *
except (ModuleNotFoundError, ImportError) as error:
    from gene_stats import GeneStats
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
                 detection_min_base_depth: float = 1.0,
                 detection_min_num_reads: int = 1,

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
                 preserve_chunks: bool = False):
                 
                 
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
        self.blastx_task = blastx_task
       # self.evalue = evalue
        self.query_min_coverage = query_min_coverage
        self.query_min_identity = query_min_identity
        self.detection_min_coverage = detection_min_coverage
        self.detection_min_identity = detection_min_identity
        self.detection_min_base_depth = detection_min_base_depth
        self.detection_min_num_reads = detection_min_num_reads

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
                    except Exception:
                        continue

                    query_coverage = ((abs(qend - qstart) + 1) / qlen) * 100 if qlen else 0.0

                    if identity >= getattr(self, 'query_min_identity', 0.0) and query_coverage >= getattr(self, 'query_min_coverage', 0.0):
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
                return False

            self.logger.info(f"Filtered {tool_name} output {output_file}: kept {kept}/{total} hits (identity>={self.detection_min_identity}, qcov>={self.query_min_coverage}%)")
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
                            read_type = determine_read_type(read_name)
                            fasta_out.write(f">{read_name}{read_type}\n{seq}\n")
                            count += 1
                if getattr(self, "verbose", True):
                    r1_count, r2_count = count_read_types(read_names)
                    self.logger.info(f"  FASTA file: {fasta_path} ({count} reads: {r1_count} R1, {r2_count} R2)")


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
                            read_type = determine_read_type(read_name)
                            fasta_out.write(f">{read_name}{read_type}\n{seq}\n")
                            count += 1

                if getattr(self, "verbose", True):
                    r1_count, r2_count = count_read_types(read_names)
                    self.logger.info(f"  FASTA file: {fasta_path} ({count} reads: {r1_count} R1, {r2_count} R2)")

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
            # Request BLAST to write tabular output to stdout so we can filter in Python
            cmd = [
                blast_cmd,
                '-query', query_arg,
                '-db', db_path,
                '-out', '-',
                '-outfmt', outfmt_fields,
                '-num_threads', str(self.threads),
                '-qcov_hsp_perc', str(int(self.query_min_coverage)),
                '-perc_identity', str(self.query_min_identity)
            ]
        elif mode == 'protein':
            blast_cmd = 'blastx'
            query_arg = '-' if gz_input else fasta_path_str
            cmd = [
                blast_cmd,
                '-query', query_arg,
                '-db', db_path,
                '-task', self.blastx_task,
                '-out', '-',
                '-outfmt', outfmt_fields,
                '-num_threads', str(self.threads),
                '-qcov_hsp_perc', str(int(self.query_min_coverage))
            ]
        else:
            self.logger.error(f"Invalid BLAST'ing' mode: {mode}")
            return False, set()

        success = False
        # If input FASTA file is large, split into chunks and run BLAST per-chunk to reduce memory/IO
        try:
            file_size = None
            try:
                file_size = os.path.getsize(str(self.input_fasta))
            except Exception:
                file_size = None

            # Determine if input FASTA is a regular file that can be chunked
            input_path_str = str(self.input_fasta) if hasattr(self, 'input_fasta') else None
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
                # # Determine concurrency now (we don't yet know num_chunks) -- cap by chunk_jobs and thread-based cap
                # if self.chunk_jobs and self.chunk_jobs > 0:
                #     if max_concurrency_by_threads is not None:
                #         concurrency = min(self.chunk_jobs, max_concurrency_by_threads)
                #     else:
                #         concurrency = self.chunk_jobs
                # else:
                #     concurrency = max_concurrency_by_threads if max_concurrency_by_threads is not None else 1
                #
                # concurrency = max(1, int(concurrency))
                #
                # # Compute threads per concurrent job distribution
                # if self.chunk_threads_per_job:
                #     threads_per_job_list = [max(1, int(self.chunk_threads_per_job))] * concurrency
                # else:
                #     base = max(1, int(self.threads) // concurrency)
                #     remainder = int(self.threads) - (base * concurrency)
                #     threads_per_job_list = [base + (1 if i < remainder else 0) for i in range(concurrency)]

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
                        if mode == 'dna':
                            per_cmd = [
                                'blastn',
                                '-query', str(chunk),
                                '-db', db_path,
                                '-out', '-',
                                '-outfmt', outfmt_fields,
                                '-num_threads', str(threads_for_this),
                                '-qcov_hsp_perc', str(int(self.query_min_coverage)),
                                '-perc_identity', str(self.query_min_identity)
                            ]
                        else:
                            per_cmd = [
                                'blastx',
                                '-query', str(chunk),
                                '-db', db_path,
                                '-task', self.blastx_task,
                                '-out', '-',
                                '-outfmt', outfmt_fields,
                                '-num_threads', str(threads_for_this),
                                '-qcov_hsp_perc', str(int(self.query_min_coverage))
                            ]

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
                '--id', str(self.query_min_identity),
                '--query-cover', str(self.query_min_coverage),
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
            '-o', str(sam_file),
            db_path,
            ] + flags + [
        ]

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
        """Parse BLAST/DIAMOND tabular output and collect per-gene stats.

        Only alignments meeting identity and query-coverage thresholds are added.
        Gene detection is decided from combined gene coverage and base-depth thresholds.
        """
        detected_genes = set()
        gene_lengths = {}  # Store gene lengths
        gene_reads = defaultdict(lambda: {'passing': [], 'all': []})  # Track all reads per gene

        # Ensure structure exists for this database and tool (do not reset entire gene_stats)
        if database not in self.gene_stats:
            self.gene_stats[database] = defaultdict(lambda: defaultdict(GeneStats))
            self.detections.setdefault(database, defaultdict(lambda: defaultdict(bool)))
        elif self.gene_stats[database].get(tool_name) is None:
            self.gene_stats[database][tool_name] = defaultdict(GeneStats)

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
                            read_name = line[1:].split(' ')[0].replace('\n','') # DIAMOND/BLAST read names split at space
                            seq_lines = []
                        else:
                            seq_lines.append(line.strip())
                    if read_name and seq_lines:
                        self.all_reads[read_name] = ''.join(seq_lines)
            except Exception as e:
                self.logger.error(f"Error reading FASTA file - Ignore if running Genefíor-Recompute without sequences")
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

                    try:
                        alignment_len = int(fields[3])
                    except Exception:
                        alignment_len = None


                    # Store gene length if available
                    if slen is not None:
                        gene_lengths[gene] = max(gene_lengths.get(gene, 0), slen)

                    # Determine whether this is a translated search (blastx/diamond-blastx)
                    tool_name_l = tool_name.lower() if isinstance(tool_name, str) else ''
                    is_translated_search = ('blastx' in tool_name_l) or ('diamond' in tool_name_l)

                    # Compute query coverage. Prefer alignment_len/qlen when available.
                    # For translated searches alignment_len is in amino-acids while qlen is in nucleotides;
                    # multiply alignment_len by 3 to convert to nucleotides before comparison.
                    query_coverage = 0.0
                    try:
                        if alignment_len is not None and qlen and qlen > 0:
                            if is_translated_search:
                                aligned_nt = alignment_len * 3
                                query_coverage = (aligned_nt / float(qlen)) * 100.0
                            else:
                                query_coverage = (alignment_len / float(qlen)) * 100.0
                        elif qstart is not None and qend is not None and qlen and qlen > 0:
                            aln_span = (abs(qend - qstart) + 1)
                            query_coverage = (aln_span / float(qlen)) * 100.0
                    except Exception:
                        query_coverage = 0.0

                    # Track all reads mapping to this gene
                    gene_reads[gene]['all'].append(read_name)
                    mapped_reads += 1

                    if identity >= self.query_min_identity and query_coverage >= self.query_min_coverage:
                        # Initialise stats if first hit for this gene
                        if gene not in self.gene_stats[database][tool_name]:
                            self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)

                        # Add hit to statistics (use subject coordinates if available)
                        ss = sstart if sstart is not None else 1
                        se = send if send is not None else ss
                        self.gene_stats[database][tool_name][gene].add_hit(
                            ss, se, identity, gene_lengths.get(gene, 0)
                        )

                        # Track reads that pass thresholds
                        gene_reads[gene]['passing'].append(read_name)
                        passing_reads += 1
        except KeyError:
            raise
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
        # Note: mapped_reads counts hits seen in the BLAST output file
        self.logger.info(f"Reads that returned a hit in {database.upper()} using {tool_name}: {mapped_reads}")
        self.logger.info(f"Reads passing thresholds in {database.upper()} using {tool_name}: {passing_reads}")
        self.logger.info(f"Detected {len(detected_genes)} genes in {database.upper()} using {tool_name}")

        # Output FASTA files of reads mapping to genes
        if self.report_fasta:
            self._write_fasta_outputs(database, tool_name, detected_genes, gene_reads, all_reads)

        return detected_genes, gene_reads

    def parse_bam_results(self, bam_file: Path, database: str, tool_name: str) -> Set[str]:
        """Parse BAM (samtools view) and collect per-gene stats.

        Uses CIGAR to compute aligned positions and per-read identity; detection
        follows the same coverage/base-depth rules as BLAST parsing.
        """
        detected_genes = set()
        if not bam_file.exists():
            self.logger.error(f"BAM file not found: {bam_file}")
            return detected_genes

        # Ensure structure exists for this database and tool (do not reset entire gene_stats)
        if database not in self.gene_stats:
            self.gene_stats[database] = defaultdict(lambda: defaultdict(GeneStats))
            self.detections.setdefault(database, defaultdict(lambda: defaultdict(bool)))
        elif self.gene_stats[database].get(tool_name) is None:
            self.gene_stats[database][tool_name] = defaultdict(GeneStats)

        gene_lengths = {}  # Store gene lengths from BAM header
        gene_reads = defaultdict(lambda: {'passing': [], 'all': [], 'passing_r1': [], 'passing_r2': []})
        all_reads = {}
        mapped_reads = 0
        passing_reads = 0


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
                aligned_positions = set()  # reference positions — feeds into gene coverage
                query_aligned_bases = 0  # query bases aligned — feeds into query coverage
                alignment_length = 0  # for identity calculation

                for count_str, op in cigar_re.findall(cigar):
                    length = int(count_str)
                    if op in ('M', '=', 'X'):
                        aligned_positions.update(range(ref_pos, ref_pos + length))
                        ref_pos += length
                        alignment_length += length
                        query_aligned_bases += length  # these bases exist in both the read AND the gene
                    elif op == 'I':  # insertion to reference
                        alignment_length += length
                        # consumes query bases, but they have no reference counterpart
                        # — do NOT count toward query_aligned_bases
                    elif op == 'D':  # deletion from reference
                        ref_pos += length
                        alignment_length += length
                        # consumes reference positions but no query bases
                        # — do NOT count toward query_aligned_bases
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
                #query_coverage = (len(aligned_positions) / query_length) * 100.0 if query_length > 0 else 0.0
                # Per-read query coverage — what proportion of this read overlaps the gene
                query_coverage = (query_aligned_bases / query_length) * 100.0 if query_length > 0 else 0.0

                if identity >= self.detection_min_identity and query_coverage >= self.query_min_coverage:
                    if gene not in self.gene_stats[database][tool_name]:
                        self.gene_stats[database][tool_name][gene] = GeneStats(gene_name=gene)
                    self.gene_stats[database][tool_name][gene].add_positions(aligned_positions, identity, gene_len)
                    gene_reads[gene]['passing'].append(read_name)

                    # Track R1/R2 separately based on read name suffix
                    if read_name.endswith('_R1') or '/1' in read_name:
                        gene_reads[gene]['passing_r1'].append(read_name)
                    elif read_name.endswith('_R2') or '/2' in read_name:
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

        except Exception as e:
            self.logger.error(f"Error reading BAM via samtools: {e}")



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
        if self.report_fasta:
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
        - Gene: Gene name
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

