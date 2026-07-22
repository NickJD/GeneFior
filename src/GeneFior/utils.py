import os, sys
import glob
import re
import logging
import gzip
import shutil
import json
import shlex
import time
import tempfile
import subprocess
import csv


def diamond_mode_label(sequence_type, genes_type=None):
    """Return the concrete DIAMOND mode used for the supplied input type."""
    if sequence_type == 'Genes-FASTA' and genes_type in ('aa', 'protein'):
        return 'DIAMOND-BLASTP'
    return 'DIAMOND-BLASTX'


def sample_temp_directory(base_temp, root_output, sample_output, sample_name,
                          user_specified=False):
    """Return an isolated temporary directory for one batch sample."""
    if not user_specified:
        return os.path.abspath(sample_output)
    safe_sample = re.sub(r'[^A-Za-z0-9._-]+', '_', str(sample_name)).strip('._')
    return os.path.join(
        os.path.abspath(base_temp or root_output),
        safe_sample or 'sample',
    )


def workflow_has_success(results):
    """Return True when at least one database/tool result completed successfully."""
    for db_results in (results or {}).values():
        for value in (db_results or {}).values():
            if isinstance(value, tuple) and value:
                success = bool(value[0])
            elif isinstance(value, (set, list, dict)):
                success = True
            else:
                success = bool(value)
            if success:
                return True
    return False


def load_gene_id_file(path, logger=None):
    """Load database gene IDs from a one-column text, CSV, or TSV file.

    If the first row contains a recognised header (Gene, Gene_ID, GeneID, or
    ID), that column is used. Otherwise the first column is treated as the gene
    ID column. Empty rows and comment rows beginning with # are ignored.
    """
    if not path:
        return set()

    file_path = os.path.abspath(os.path.expanduser(str(path)))
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Gene ID list not found: {file_path}")

    with open(file_path, 'r', newline='') as handle:
        raw_lines = [
            line for line in handle
            if line.strip() and not line.lstrip().startswith('#')
        ]

    if not raw_lines:
        if logger:
            logger.warning(f"Gene ID list is empty: {file_path}")
        return set()

    first_line = raw_lines[0]
    delimiter = '\t' if '\t' in first_line else ','
    rows = list(csv.reader(raw_lines, delimiter=delimiter))
    header_names = {'gene', 'gene_id', 'gene id', 'geneid', 'id'}
    header = [cell.strip().lower() for cell in rows[0]]
    id_column = 0
    start_index = 0
    for index, name in enumerate(header):
        if name in header_names:
            id_column = index
            start_index = 1
            break

    gene_ids = set()
    for row in rows[start_index:]:
        if id_column >= len(row):
            continue
        gene_id = row[id_column].strip()
        if gene_id:
            gene_ids.add(gene_id)

    if logger:
        logger.info(f"Loaded {len(gene_ids)} always-flag gene ID(s) from {file_path}")
    return gene_ids


def run_fastp_for_paired_reads(options, r1_path, r2_path, logger):
    """Trim paired FASTQ reads with fastp and return the generated read paths."""
    fastp_path = shutil.which('fastp')
    if not fastp_path:
        message = (
            "fastp trimming was requested, but the 'fastp' executable was not found "
            "in PATH. Install fastp or rerun without --fastp-trim."
        )
        logger.error(message)
        raise RuntimeError(message)

    output_dir = os.path.abspath(options.output)
    temp_root = os.path.abspath(getattr(options, 'temp_directory', None) or output_dir)
    trim_dir = os.path.join(temp_root, 'fastp_trimmed_reads')
    if os.path.normcase(temp_root) != os.path.normcase(output_dir):
        sample_name = re.sub(r'[^A-Za-z0-9._-]+', '_', os.path.basename(output_dir))
        trim_dir = os.path.join(trim_dir, sample_name or 'sample')
    report_dir = os.path.join(output_dir, 'fastp')
    os.makedirs(trim_dir, exist_ok=True)
    os.makedirs(report_dir, exist_ok=True)

    def _read_stem(path):
        name = os.path.basename(path)
        for suffix in ('.fastq.gz', '.fq.gz', '.fastq', '.fq', '.gz'):
            if name.lower().endswith(suffix):
                return name[:-len(suffix)]
        return name

    trimmed_r1 = os.path.join(trim_dir, f"R1_{_read_stem(r1_path)}_fastp_trimmed.fastq.gz")
    trimmed_r2 = os.path.join(trim_dir, f"R2_{_read_stem(r2_path)}_fastp_trimmed.fastq.gz")
    json_report = os.path.join(report_dir, 'fastp.json')
    html_report = os.path.join(report_dir, 'fastp.html')

    generated_paths = (trimmed_r1, trimmed_r2, json_report, html_report)
    for path in generated_paths:
        try:
            if os.path.isfile(path):
                os.remove(path)
        except OSError as exc:
            message = f"Could not replace existing fastp output '{path}': {exc}"
            logger.error(message)
            raise RuntimeError(message) from exc

    try:
        requested_threads = max(1, int(getattr(options, 'threads', 4)))
    except (TypeError, ValueError):
        requested_threads = 4
    fastp_threads = min(requested_threads, 16)
    if fastp_threads != requested_threads:
        logger.info(
            f"fastp will use {fastp_threads} threads (requested {requested_threads}; "
            "fastp supports at most 16 worker threads)"
        )

    command = [
        fastp_path,
        '--in1', r1_path,
        '--in2', r2_path,
        '--out1', trimmed_r1,
        '--out2', trimmed_r2,
        '--thread', str(fastp_threads),
        '--detect_adapter_for_pe',
        '--json', json_report,
        '--html', html_report,
    ]

    logger.info("fastp auto trimming enabled for paired FASTQ input")
    logger.info(f"Running fastp: {shlex.join(command)}")
    try:
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except OSError as exc:
        message = f"Could not launch fastp: {exc}"
        logger.error(message)
        raise RuntimeError(message) from exc

    if completed.returncode != 0:
        for path in generated_paths:
            try:
                if os.path.isfile(path):
                    os.remove(path)
            except OSError:
                pass
        detail = (completed.stderr or completed.stdout or '').strip()
        message = f"fastp failed with exit code {completed.returncode}"
        if detail:
            message = f"{message}: {detail}"
        logger.error(message)
        raise RuntimeError(message)

    missing_outputs = [
        path for path in (trimmed_r1, trimmed_r2)
        if not os.path.isfile(path) or os.path.getsize(path) == 0
    ]
    if missing_outputs:
        for path in generated_paths:
            try:
                if os.path.isfile(path):
                    os.remove(path)
            except OSError:
                pass
        message = (
            "fastp completed without producing valid paired FASTQ output(s): "
            + ', '.join(missing_outputs)
        )
        logger.error(message)
        raise RuntimeError(message)

    options.fastp_report_json = json_report
    options.fastp_report_html = html_report
    options.input_fastq_trimmed = (trimmed_r1, trimmed_r2)
    logger.info("Downstream analysis will use fastp-trimmed paired reads")
    for label, path in (('JSON', json_report), ('HTML', html_report)):
        if os.path.isfile(path):
            logger.info(f"fastp {label} report: {path}")
        else:
            logger.warning(f"fastp did not create its expected {label} report: {path}")
    return trimmed_r1, trimmed_r2


def prepare_fastq_for_alignment(r1_path, r2_path, temp_dir, logger, force: bool = False):
    # Check if FASTQ read IDs need suffixes and create temporary modified files if needed.

    # Check if suffixes are needed
    needs_suffix = check_read_id_uniqueness(r1_path, r2_path, logger, sample_size=10000)

    if not needs_suffix and not force:
        logger.info("FASTQ read IDs are already unique, no modification needed")
        return r1_path, r2_path, False
    if not needs_suffix and force:
        logger.info("FASTQ read IDs appear unique but forced modification requested — creating modified FASTQ files with /1 and /2 suffixes")

    logger.warning(
        "FASTQ read IDs are not unique after first space. Creating/modifying FASTQ files with /1 and /2 suffixes...")

    # Decide where to place (or look for) modified files. If a temp_dir is provided, prefer it so the caller
    # can control lifecycle. Otherwise place modified files alongside the originals so they can be reused
    # across runs.
    if temp_dir:
        work_dir = temp_dir
    else:
        work_dir = os.path.dirname(os.path.abspath(r1_path))

    os.makedirs(work_dir, exist_ok=True)

    # Create expected output paths (consistent naming so they can be detected later)
    # Normalise basename and strip common suffixes (including gzipped extensions)
    def _strip_fastq_suffix(fn):
        bn = os.path.basename(fn)
        for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq', '.gz'):
            if bn.endswith(ext):
                bn = bn[:-len(ext)]
                break
        return bn

    base_r1 = _strip_fastq_suffix(r1_path)
    base_r2 = _strip_fastq_suffix(r2_path)
    r1_modified = os.path.join(work_dir, f"{base_r1}_R1_modified.fastq.gz")
    r2_modified = os.path.join(work_dir, f"{base_r2}_R2_modified.fastq.gz")

    def _file_looks_modified(path, expected_suffix):
        # Quick validation: open the file and inspect the first header line to ensure the expected suffix exists
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            return False
        opener = gzip.open if path.endswith('.gz') else open
        mode = 'rt' if path.endswith('.gz') else 'r'
        try:
            with opener(path, mode) as fh:
                # Read until first header line
                for i in range(4):
                    line = fh.readline()
                    if not line:
                        break
                    if line.startswith('@'):
                        read_id = line[1:].strip().split()[0]
                        # Accept suffix appended (e.g. ID/1 or ID/2) or suffix present anywhere in the base id
                        if read_id.endswith(expected_suffix) or expected_suffix in read_id:
                            return True
                        return False
        except Exception:
            return False
        return False

    # If candidate modified files already exist and look valid, reuse them
    if _file_looks_modified(r1_modified, '/1') and _file_looks_modified(r2_modified, '/2'):
        logger.info(f"Found existing modified FASTQ files, reusing: {r1_modified}, {r2_modified}")
        return r1_modified, r2_modified, False

    # Otherwise create them in the chosen work_dir
    logger.info(f"Creating modified FASTQ files: {r1_path} -> {r1_modified} ; {r2_path} -> {r2_modified}")

    def modify_fastq_ids(input_path, output_path, suffix):
        # Add suffix to FASTQ read IDs
        # Try an optimised shell pipeline using awk + pigz when available. This
        # avoids Python-level per-line processing for very large files.
        pigz_path = shutil.which('pigz')
        gzip_path = shutil.which('gzip')
        compressor = pigz_path or gzip_path
        decompressor = 'gzip -dc' if input_path.endswith('.gz') else 'cat'
        # AWK program: for every header line (NR%4==1) strip leading @, split at
        # first space, append suffix to id and reprint, otherwise print as-is.
        awk_prog = (
            "{\n"
            " if (NR%4==1) {\n"
            "   line=substr($0,2); split(line,a,\" \"); id=a[1]; rest=\"\"; pos=index($0,\" \");\n"
            "   if (pos>0) rest=substr($0,pos);\n"
            "   printf(\"@%s%s%s\\n\", id, suffix, rest);\n"
            " } else { print $0 }\n"
            "}\n"
        )

        if compressor:
            # Build shell command and run it. Use -v to pass suffix safely to awk.
            cmd = f"{decompressor} {shlex.quote(input_path)} | awk -v suffix='{suffix}' '{awk_prog}' | {shlex.quote(compressor)} -c > {shlex.quote(output_path)}"
            try:
                subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
                return
            except Exception:
                # Fall back to Python implementation on any failure
                pass

        # Fallback Python implementation (safe, portable)
        opener = gzip.open if input_path.endswith('.gz') else open
        input_mode = 'rt' if input_path.endswith('.gz') else 'r'

        with opener(input_path, input_mode) as infile, gzip.open(output_path, 'wt') as outfile:
            line_num = 0
            for line in infile:
                line_num += 1
                if line_num % 4 == 1:  # Header line
                    # Modify read ID
                    read_id = line[1:].strip()  # Remove @ and whitespace
                    if ' ' in read_id:
                        parts = read_id.split(' ', 1)
                        modified_line = f"@{parts[0]}{suffix} {parts[1]}\n"
                    else:
                        modified_line = f"@{read_id}{suffix}\n"
                    outfile.write(modified_line)
                else:
                    outfile.write(line)

    # Modify both files
    try:
        modify_fastq_ids(r1_path, r1_modified, '/1')
        modify_fastq_ids(r2_path, r2_modified, '/2')
    except Exception as e:
        logger.error(f"Failed to create modified FASTQ files: {e}")
        raise

    return r1_modified, r2_modified, True

def check_read_id_uniqueness(r1_path, r2_path, logger, sample_size=10000):
    """
    Check if read IDs remain unique after truncation at first space.

    Returns True if suffixes are required (i.e. some base IDs collide), False
    otherwise.
    """
    # More robust: collect a sample of base IDs from R1 and R2 separately and
    # check whether any base ID appears in both files. This directly tests the
    # paired-read collision case which requires /1 and /2 suffixes.
    def _iter_base_ids(path, limit):
        opener = gzip.open if path.endswith(('.gz', '.gzip')) else open
        mode = 'rt' if path.endswith(('.gz', '.gzip')) else 'r'
        count = 0
        try:
            with opener(path, mode) as fh:
                for line in fh:
                    if not line:
                        break
                    if line.startswith('@'):
                        header = line[1:].strip()
                        base = header.split()[0] if ' ' in header else header
                        yield base
                        count += 1
                        if limit and count >= limit:
                            break
        except Exception as e:
            logger.debug(f"Failed to read headers from {path}: {e}")
            return

    # Read up to sample_size headers from each file (split evenly if possible)
    half = max(1, sample_size // 2)
    r1_ids = set()
    r2_ids = set()
    try:
        for i, bid in enumerate(_iter_base_ids(r1_path, half)):
            if not bid:
                continue
            r1_ids.add(bid)
    except Exception:
        logger.debug("Error sampling R1 headers for uniqueness check")

    try:
        for i, bid in enumerate(_iter_base_ids(r2_path, sample_size - half)):
            if not bid:
                continue
            r2_ids.add(bid)
    except Exception:
        logger.debug("Error sampling R2 headers for uniqueness check")

    # If there's any intersection between sets, we need suffixes
    needs_suffix = len(r1_ids & r2_ids) > 0
    logger.debug(f"check_read_id_uniqueness: sampled {len(r1_ids)} R1 ids and {len(r2_ids)} R2 ids; intersection={len(r1_ids & r2_ids)}")
    return needs_suffix

def requires_fasta_conversion(tools):
    # BLAST/DIAMOND tools need FASTA; hmmer_protein (hmmsearch) also searches a FASTA target file.
    return any(tool in ('blastn', 'blastx', 'diamond', 'all', 'hmmer_protein') for tool in (tools or []))


# Module-level canonical constants for extensions and suffixes
FASTA_EXTS = ('.fa', '.fasta', '.fna', '.ffn', '.faa', '.fas')
FASTQ_EXTS = ('.fastq', '.fq')
GZ_VARIANTS = ('.gz', '.gzip')
R1_SUFFIXES = ['_r1', '_1', '.r1', '.1', '-r1', '-1']
R2_SUFFIXES = ['_r2', '_2', '.r2', '.2', '-r2', '-2']


def get_max_fasta_chunk_bytes(mb_value):
    """Convert MiB value to bytes, with safe fallback."""
    try:
        return int(float(mb_value) * 1024.0 * 1024.0)
    except Exception:
        return 200 * 1024 * 1024


def add_file_handler_for_path(logger, path, level=logging.INFO):
    """Create a FileHandler for `path`, attach it to `logger`, and return the handler.

    The caller is responsible for removing/closing the handler via remove_file_handler.
    """
    try:
        parent = os.path.dirname(str(path))
        if parent:
            os.makedirs(parent, exist_ok=True)
    except Exception:
        pass
    try:
        fh = logging.FileHandler(str(path), mode='a')
        fh.setLevel(level)
        fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(fh)
        return fh
    except Exception:
        return None


def remove_file_handler(logger, handler):
    try:
        if handler:
            try:
                logger.removeHandler(handler)
            except Exception:
                pass
            try:
                handler.close()
            except Exception:
                pass
    except Exception:
        pass

def FASTQ_to_FASTA(options, logger):
    # Convert paired FASTQ -> single combined FASTA when needed.
    # Prefers seqtk + compressor (pigz/gzip) for speed; falls back to the Python streamer.

    def find_pair(input_spec):
        if ',' in input_spec:
            r1, r2 = map(str.strip, input_spec.split(',', 1))
            return r1, r2
        base = input_spec
        candidates = [base, base + '_R1.fastq', base + '_R1.fq', base + '_1.fastq', base + '_1.fq']
        r1_path = None
        for c in candidates:
            if os.path.isfile(c):
                r1_path = c
                break
        if not r1_path:
            logger.error("Could not locate R1 FASTQ. Provide `R1.fastq,R2.fastq` as `-i`.")
            sys.exit(1)
        if '_R1.' in r1_path:
            r2_path = r1_path.replace('_R1.', '_R2.')
        elif '_1.' in r1_path:
            r2_path = r1_path.replace('_1.', '_2.')
        else:
            r2_path = r1_path.replace('_R1', '_R2')
        return r1_path, r2_path

    r1_path, r2_path = find_pair(options.input)
    if not os.path.isfile(r1_path) or not os.path.isfile(r2_path):
        logger.error(f"Paired FASTQ files not found or not regular files: {r1_path}, {r2_path}")
        sys.exit(1)

    conv_dir = os.path.join(options.output, 'paired_fastq_fasta')
    os.makedirs(conv_dir, exist_ok=True)
    combined_fasta = os.path.join(conv_dir, 'fastq_to_fasta_combined.fasta.gz') # Could need to modify in future
    # Do not automatically reuse a combined FASTA just because it exists; prefer
    # to validate using metadata (later) so we don't silently reuse an outdated
    # combined FASTA created from different inputs or settings.

    # Prefer seqtk + compressor when available (fast C tools). Python streamer
    # is used only as a fallback or when external tools fail.
    def fastq_to_fasta_py(fastq_path, out_handle, suffix=None):
        opener = gzip.open if fastq_path.endswith(('.gz', '.gzip')) else open
        mode = 'rt' if fastq_path.endswith(('.gz', '.gzip')) else 'r'
        try:
            with opener(fastq_path, mode) as fh:
                while True:
                    header = fh.readline()
                    if not header:
                        break
                    seq = fh.readline()
                    plus = fh.readline()
                    qual = fh.readline()
                    if not seq:
                        break
                    # Normalise
                    header = header.rstrip('\n')
                    seq = seq.rstrip('\n')
                    read_id = header[1:].strip() if header.startswith('@') else header.strip()
                    if suffix:
                        if ' ' in read_id:
                            parts = read_id.split(' ', 1)
                            out_header = f">{parts[0]}{suffix} {parts[1]}\n"
                        else:
                            out_header = f">{read_id}{suffix}\n"
                    else:
                        if ' ' in read_id:
                            parts = read_id.split(' ', 1)
                            out_header = f">{parts[0]} {parts[1]}\n"
                        else:
                            out_header = f">{read_id}\n"
                    out_handle.write(out_header.encode('utf-8'))
                    out_handle.write((seq + '\n').encode('utf-8'))
        except Exception as e:
            logger.error(f"FASTQ->FASTA conversion failed for {fastq_path}: {e}")
            raise

    # Check if we need suffixes
    needs_suffix = check_read_id_uniqueness(r1_path, r2_path, logger)

    # Define suffixes to use
    r1_suffix = '/1' if needs_suffix else None
    r2_suffix = '/2' if needs_suffix else None

    if needs_suffix:
        logger.warning(
            f"Read IDs are not unique after first space. Appending {r1_suffix} and {r2_suffix} to read names.")

    # If combined FASTA already exists, validate it using a small metadata file
    meta_path = combined_fasta + '.meta.json'
    force_regen = bool(getattr(options, 'force_regenerate_fasta', False))
    def _file_meta_matches(meta_file, r1, r2, r1_suf, r2_suf):
        try:
            if not os.path.isfile(meta_file):
                return False
            with open(meta_file, 'r') as mh:
                meta = json.load(mh)
            # Compare paths, sizes and mtimes
            if meta.get('r1_path') != os.path.abspath(r1):
                return False
            if meta.get('r2_path') != os.path.abspath(r2):
                return False
            if meta.get('r1_size') != os.path.getsize(r1):
                return False
            if meta.get('r2_size') != os.path.getsize(r2):
                return False
            if meta.get('r1_mtime') != os.path.getmtime(r1):
                return False
            if meta.get('r2_mtime') != os.path.getmtime(r2):
                return False
            if meta.get('r1_suffix') != r1_suf or meta.get('r2_suffix') != r2_suf:
                return False
            return True
        except Exception:
            return False

    # If metadata matches, reuse. If metadata is missing but the combined FASTA exists
    # and is newer than the input FASTQ files, assume it was created from these inputs
    # and reuse it (write metadata so future runs are explicit). This avoids
    # unnecessary re-conversion when users re-run the pipeline in the same folder.
    if os.path.isfile(combined_fasta) and os.path.getsize(combined_fasta) > 0:
        if not force_regen and _file_meta_matches(meta_path, r1_path, r2_path, r1_suffix, r2_suffix):
            logger.info(f"Found existing validated combined FASTA at `{combined_fasta}`; using it (skipping conversion)")
            options.input_fasta = combined_fasta
            options.input_fastq = (r1_path, r2_path)
            return
        # metadata missing or doesn't match; attempt a conservative heuristic:
        try:
            combined_mtime = os.path.getmtime(combined_fasta)
            r1_mtime = os.path.getmtime(r1_path) if os.path.exists(r1_path) else 0
            r2_mtime = os.path.getmtime(r2_path) if os.path.exists(r2_path) else 0
            # If the combined FASTA is newer than both FASTQ inputs, assume it was produced
            # from them and reuse it. This is a pragmatic heuristic to avoid unnecessary work.
            if not force_regen and combined_mtime >= max(r1_mtime, r2_mtime):
                logger.info(f"Found existing combined FASTA `{combined_fasta}` newer than inputs; reusing it and writing metadata")
                # attempt to write metadata for future runs
                try:
                    meta = {
                        'r1_path': os.path.abspath(r1_path),
                        'r2_path': os.path.abspath(r2_path),
                        'r1_size': os.path.getsize(r1_path) if os.path.exists(r1_path) else None,
                        'r2_size': os.path.getsize(r2_path) if os.path.exists(r2_path) else None,
                        'r1_mtime': r1_mtime,
                        'r2_mtime': r2_mtime,
                        'r1_suffix': r1_suffix,
                        'r2_suffix': r2_suffix,
                        'created': time.time()
                    }
                    with open(meta_path, 'w') as mh:
                        json.dump(meta, mh)
                except Exception:
                    logger.debug("Failed to write metadata for existing combined FASTA; continuing to reuse without metadata")
                options.input_fasta = combined_fasta
                options.input_fastq = (r1_path, r2_path)
                return
        except Exception:
            # If any check fails, fall through to conversion logic
            pass
    if force_regen:
        logger.info("Force regeneration of combined FASTA requested (options.force_regenerate_fasta=True)")

    # Try optimised path: prefer seqtk + pigz when available (fast C tools).
    # If not available or an error occurs, fall back to the native Python streamer.
    seqtk_path = shutil.which('seqtk')
    pigz_path = shutil.which('pigz')
    gzip_path = shutil.which('gzip')
    logger.info(f"seqtk found: {bool(seqtk_path)} ({seqtk_path}) ; pigz found: {bool(pigz_path)} ({pigz_path}) ; gzip found: {bool(gzip_path)} ({gzip_path})")
    used_optimised = False
    # Prefer seqtk + pigz, but fall back to seqtk + gzip if pigz is not available
    compressor = pigz_path or gzip_path
    if seqtk_path and compressor:
        logger.info("Attempting optimised seqtk+pigz FASTQ->FASTA conversion")
        try:
            # If no header suffix modification needed we can run a simple shell pipeline
            if not r1_suffix and not r2_suffix:
                cmd = f"(seqtk seq -A {shlex.quote(r1_path)}; seqtk seq -A {shlex.quote(r2_path)}) | {shlex.quote(compressor)} -c > {shlex.quote(combined_fasta)}"
                subprocess.run(cmd, shell=True, check=True)
                used_optimised = True
                logger.info(f"Optimised seqtk+{os.path.basename(compressor)} conversion succeeded (no header modification needed)")
            else:
                # Need to modify headers; stream seqtk outputs through a small modifier into pigz
                with open(combined_fasta, 'wb') as outf:
                    comp_proc = subprocess.Popen([compressor, '-c'], stdin=subprocess.PIPE, stdout=outf)
                try:
                    for fastq_path, suffix in ((r1_path, r1_suffix), (r2_path, r2_suffix)):
                        # Spawn seqtk for this fastq and stream its stdout
                        proc = subprocess.Popen([seqtk_path, 'seq', '-A', fastq_path], stdout=subprocess.PIPE)
                        assert proc.stdout is not None
                        buf = b''
                        for chunk in iter(lambda: proc.stdout.read(8192), b''):
                            if not chunk:
                                break
                            buf += chunk
                            parts = buf.split(b'\n')
                            buf = parts.pop()  # last partial
                            for part in parts:
                                if part.startswith(b'>'):
                                    header = part[1:].decode('utf-8', errors='ignore').rstrip()
                                    if suffix:
                                        if ' ' in header:
                                            hparts = header.split(' ', 1)
                                            out_line = f">{hparts[0]}{suffix} {hparts[1]}\n"
                                        else:
                                            out_line = f">{header}{suffix}\n"
                                    else:
                                        out_line = f">{header}\n"
                                    comp_proc.stdin.write(out_line.encode('utf-8'))
                                else:
                                    comp_proc.stdin.write(part + b'\n')
                        # Flush any remaining buffer
                        if buf:
                            part = buf
                            if part.startswith(b'>'):
                                header = part[1:].decode('utf-8', errors='ignore').rstrip()
                                if suffix:
                                    if ' ' in header:
                                        hparts = header.split(' ', 1)
                                        out_line = f">{hparts[0]}{suffix} {hparts[1]}\n"
                                    else:
                                        out_line = f">{header}{suffix}\n"
                                else:
                                    out_line = f">{header}\n"
                                comp_proc.stdin.write(out_line.encode('utf-8'))
                            else:
                                comp_proc.stdin.write(part + b'\n')
                        proc.stdout.close()
                        proc.wait()
                    # close pigz stdin to finish compression
                    comp_proc.stdin.close()
                    comp_proc.wait()
                    used_optimised = True
                    logger.info(f"optimised seqtk+{os.path.basename(compressor)} conversion succeeded (with header modification)")
                except Exception:
                        try:
                            comp_proc.kill()
                        except Exception:
                            pass
                        raise
        except Exception as e:
            logger.warning(f"seqtk+{os.path.basename(compressor) if compressor else 'compressor'} optimised path failed, falling back to Python converter: {e}")

    if not used_optimised:
        logger.info("Using native Python FASTQ->FASTA converter (seqtk+pigz not used)")
        # Convert both FASTQ files and append into a single FASTA (native Python streamer)
        with gzip.open(combined_fasta, 'wb') as out:
            fastq_to_fasta_py(r1_path, out, suffix=r1_suffix)
            fastq_to_fasta_py(r2_path, out, suffix=r2_suffix)

    # Write metadata to enable safe reuse in future runs
    try:
        meta = {
            'r1_path': os.path.abspath(r1_path),
            'r2_path': os.path.abspath(r2_path),
            'r1_size': os.path.getsize(r1_path),
            'r2_size': os.path.getsize(r2_path),
            'r1_mtime': os.path.getmtime(r1_path),
            'r2_mtime': os.path.getmtime(r2_path),
            'r1_suffix': r1_suffix,
            'r2_suffix': r2_suffix,
            'created': time.time()
        }
        with open(meta_path, 'w') as mh:
            json.dump(meta, mh)
    except Exception:
        # Non-fatal — conversion succeeded but we couldn't write metadata
        logger.debug("Failed to write FASTA metadata file; future runs may re-create the combined FASTA")

    logger.info(f"Combined FASTA created at {combined_fasta}")
    options.input_fasta = combined_fasta
    options.input_fastq = (r1_path, r2_path)

def copy_to_temp_directory(source_path, temp_dir, logger):
    # Copy a file to the temporary directory for faster I/O.
    if not temp_dir:
        return source_path

    if not os.path.isfile(source_path):
        logger.error(f"Source file does not exist or is not a regular file: {source_path}")
        return source_path

    try:
        # Create temp directory if it doesn't exist
        os.makedirs(temp_dir, exist_ok=True)

        # Get filename and create destination path
        filename = os.path.basename(source_path)
        dest_path = os.path.join(temp_dir, filename)

        # Check if file already exists in temp directory
        if os.path.exists(dest_path):
            # Verify it's the same file (by size as a quick check)
            if os.path.getsize(source_path) == os.path.getsize(dest_path):
                logger.info(f"File already exists in temp directory: {dest_path}")
                return dest_path
            else:
                logger.warning(f"Existing temp file has different size, overwriting: {dest_path}")

        # Copy file to temp directory
        logger.info(f"Copying {source_path} to temp directory {temp_dir}...")
        file_size_mb = os.path.getsize(source_path) / (1024 * 1024)
        logger.info(f"File size: {file_size_mb:.2f} MB")

        shutil.copy2(source_path, dest_path)
        logger.info(f"Successfully copied to: {dest_path}")

        return dest_path

    except Exception as e:
        logger.error(f"Failed to copy file to temp directory: {e}")
        logger.warning(f"Falling back to original file: {source_path}")
        return source_path


def cleanup_temp_files(temp_dir, files_to_remove, logger):
    # Remove temporary files created in the temp directory
    if not temp_dir:
        return
    cleaned_count = 0
    # Normalise temp_dir to absolute path for safe containment checks
    try:
        abs_temp_dir = os.path.abspath(temp_dir)
    except Exception:
        abs_temp_dir = temp_dir

    for file_path in files_to_remove:
        try:
            if not file_path:
                continue
            fp = str(file_path)
            if not os.path.exists(fp):
                continue
            # Ensure file is within the temp directory (avoid substring pitfalls)
            abs_fp = os.path.abspath(fp)
            if not abs_fp.startswith(abs_temp_dir.rstrip(os.path.sep) + os.path.sep) and abs_fp != abs_temp_dir:
                # Not inside temp dir; skip
                continue
            try:
                os.remove(abs_fp)
                logger.info(f"Removed temporary file: {abs_fp}")
                cleaned_count += 1
            except Exception as e:
                logger.warning(f"Failed to remove temporary file {abs_fp}: {e}")
        except Exception:
            continue

    if cleaned_count > 0:
        logger.info(f"Cleaned up {cleaned_count} temporary file(s)")


def validate_paired_fastq(options, logger):
    # Validate paired FASTQ input files.
    # Parse the input specification
    raw_inputs = [p.strip() for p in options.input.split(',')] if isinstance(options.input, str) else list(
        options.input)

    if len(raw_inputs) != 2:
        logger.error("For Paired-FASTQ input, please provide exactly two FASTQ files separated by a comma")
        print("Error: For Paired-FASTQ input, please provide exactly two FASTQ files separated by a comma",
              file=sys.stderr)
        sys.exit(1)

    r1_path, r2_path = raw_inputs

    # Check if R1 exists and is a file
    if not os.path.isfile(r1_path):
        logger.error(f"R1 input file '{r1_path}' not found or is not a regular file")
        print(f"Error: Input file '{r1_path}' not found", file=sys.stderr)
        sys.exit(1)

    # Check if R2 exists and is a file
    if not os.path.isfile(r2_path):
        logger.error(f"R2 input file '{r2_path}' not found or is not a regular file")
        print(f"Error: Input file '{r2_path}' not found", file=sys.stderr)
        sys.exit(1)

    # Check if they're the same file
    if os.path.abspath(r1_path) == os.path.abspath(r2_path):
        logger.error("R1 and R2 FASTQ files cannot be the same file")
        print("Error: R1 and R2 FASTQ files cannot be the same file", file=sys.stderr)
        sys.exit(1)

    logger.info(f"Validated paired FASTQ files:")
    logger.info(f"  R1: {r1_path} ({os.path.getsize(r1_path) / (1024 * 1024):.2f} MB)")
    logger.info(f"  R2: {r2_path} ({os.path.getsize(r2_path) / (1024 * 1024):.2f} MB)")

    return r1_path, r2_path


def discover_samples_from_input_dir(input_dir, sequence_type, logger):
    """Discover samples in a flat input directory.
    Returns list of (sample_name, input_spec) where input_spec is either a FASTA path
    or a comma-separated R1,R2 FASTQ specification compatible with existing CLI.
    """
    samples = []
    if not input_dir:
        return samples
    if not os.path.isdir(input_dir):
        logger.error(f"Input directory not found: {input_dir}")
        return samples

    fasta_exts = FASTA_EXTS
    fastq_exts = FASTQ_EXTS

    def _strip_name_no_gz(fn, exts):
        """Return the basename without known extension(s) or trailing .gz/gzip.
        Preserves original case except for removed suffixes.
        """
        bn = os.path.basename(fn)
        lbn = bn.lower()
        # remove .gz or .gzip if present
        if lbn.endswith('.gz'):
            bn_no_gz = bn[:-3]
            lbn_no_gz = lbn[:-3]
        elif lbn.endswith('.gzip'):
            bn_no_gz = bn[:-5]
            lbn_no_gz = lbn[:-5]
        else:
            bn_no_gz = bn
            lbn_no_gz = lbn

        for ext in exts:
            if lbn_no_gz.endswith(ext):
                return bn_no_gz[:-len(ext)]
        # fallback to splitting last extension
        return os.path.splitext(bn_no_gz)[0]

    entries = sorted(os.listdir(input_dir))

    if sequence_type == 'Single-FASTA':
        for fn in entries:
            fp = os.path.join(input_dir, fn)
            if not os.path.isfile(fp):
                continue
            ln = fn.lower()
            if any(ln.endswith(ext) for ext in fasta_exts) or any(ln.endswith(ext + '.gz') or ln.endswith(ext + '.gzip') for ext in fasta_exts):
                name = _strip_name_no_gz(fn, fasta_exts)
                samples.append((name, fp))
        return samples

    # Paired-FASTQ: find R1/R2 pairs by common suffixes
    # Build map of filenames (ignore hidden/system files like .DS_Store)
    files = [f for f in entries if os.path.isfile(os.path.join(input_dir, f)) and not f.startswith('.') and f.lower() not in ('.ds_store', 'thumbs.db')]
    # Map lowercased filename -> original filename for case-insensitive lookups
    lc_to_orig = {f.lower(): f for f in files}

    used = set()
    for f in files:
        if f in used:
            continue
        lc = f.lower()
        # candidate R1 patterns (common naming conventions)
        candidate_suffixes = ['_r1', '_1', '.r1', '.1', '-r1', '-1']
        found = False
        for suf in candidate_suffixes:
            for ext in fastq_exts:
                for gz in ['', '.gz', '.gzip']:
                    pattern = suf + ext + gz
                    if lc.endswith(pattern):
                        # Work in lowercased space for reliable prefix computation
                        prefix_lc = lc[:-len(pattern)]
                        # try to find corresponding R2 (lowercased lookup)
                        for r2suf in ['_r2', '_2', '.r2', '.2', '-r2', '-2']:
                            r2_lcname = prefix_lc + r2suf + ext + gz
                            r2orig = lc_to_orig.get(r2_lcname)
                            if r2orig:
                                r1p = os.path.join(input_dir, f)
                                r2p = os.path.join(input_dir, r2orig)
                                # sample name from prefix (preserve original case)
                                prefix_orig = f[:-len(pattern)] if len(pattern) > 0 else f
                                sample_name = prefix_orig if prefix_orig else _strip_name_no_gz(f, fastq_exts)
                                # Normalise sample name (remove trailing separators)
                                sample_name = sample_name.rstrip('._-')
                                samples.append((sample_name, f"{r1p},{r2p}"))
                                used.add(f)
                                used.add(r2orig)
                                found = True
                                break
                        if found:
                            break
                if found:
                    break
            if found:
                break
        # If not found via strict suffix matching, attempt a more flexible regex-based
        # pairing that handles names like sample_1_trimmed.fastq.gz -> sample_2_trimmed.fastq.gz
        if not found:
            try:
                # Determine the full extension (preserve original ext to build candidate names)
                ext_full = None
                for cand_ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
                    if lc.endswith(cand_ext):
                        ext_full = f[-len(cand_ext):]
                        base_orig = f[:-len(cand_ext)]
                        base_lc = lc[:-len(cand_ext)]
                        break
                if not ext_full:
                    # Unknown extension, skip
                    continue

                # Regex: capture prefix, an R1/1 token (with optional surrounding separators), and optional tail
                # Use greedy prefix so we capture the last occurrence of the R1 token
                m = re.match(r"(?i)^(?P<prefix>.*)(?P<suf>(?:[_\-]r?1|r?1))(?P<tail>[_\-\w]*)$", base_lc)
                if m:
                    prefix_lc = m.group('prefix')
                    tail_lc = m.group('tail') or ''
                    # Try alternative R2 suffixes
                    for r2suf in ['_r2', '_2', '-r2', '-2', 'r2', '2']:
                        r2_lcname = prefix_lc + r2suf + tail_lc + ext_full
                        r2orig = lc_to_orig.get(r2_lcname)
                        if r2orig:
                            r1p = os.path.join(input_dir, f)
                            r2p = os.path.join(input_dir, r2orig)
                            # preserve original-case sample name using original base
                            prefix_orig = base_orig[:len(base_orig) - len(m.group('suf') + (m.group('tail') or ''))]
                            sample_name = prefix_orig if prefix_orig else _strip_name_no_gz(f, fastq_exts)
                            sample_name = sample_name.rstrip('._-')
                            samples.append((sample_name, f"{r1p},{r2p}"))
                            used.add(f)
                            used.add(r2orig)
                            found = True
                            break
                        # also try with .gz variant keys if ext_full didn't include .gz originally
                    # end for r2suf
            except Exception:
                # Best-effort fallback; on any error, ignore and continue
                pass

    # As a fallback, attempt to pair by common prefix before the first dot/underscore
    if not samples:
        prefixes = {}
        for f in files:
            # strip gz and fastq extensions so base tokenisation is correct
            base = _strip_name_no_gz(f, fastq_exts)
            # split on common delimiters
            for sep in ('_','.'):
                if sep in base:
                    base = base.split(sep)[0]
                    break
            prefixes.setdefault(base, []).append(f)
        for base, flist in prefixes.items():
            if len(flist) >= 2:
                # attempt to select one R1 and one R2 by reading names
                r1 = None
                r2 = None
                for ff in flist:
                    if 'r1' in ff.lower() or ff.lower().endswith('_1') or ff.lower().endswith('.1'):
                        r1 = ff
                    elif 'r2' in ff.lower() or ff.lower().endswith('_2') or ff.lower().endswith('.2'):
                        r2 = ff
                if r1 and r2:
                    samples.append((base, os.path.join(input_dir, r1) + ',' + os.path.join(input_dir, r2)))

    return samples


def discover_samples_from_subdirs(parent_dir, sequence_type, logger, exclude_paths=None):
    """Discover samples where each subdirectory contains one sample (one FASTA or a paired FASTQ set).
    Returns list of (sample_name, input_spec).

    exclude_paths: optional iterable of directory paths (strings). Any subdirectory whose absolute
    path matches one of the exclude_paths will be skipped. This is intended to prevent the
    discovery routine from scanning an output directory that lives inside the parent_dir.
    """
    samples = []
    if not parent_dir or not os.path.isdir(parent_dir):
        logger.error(f"Input subdirs path not found: {parent_dir}")
        return samples
    # Normalise exclude paths to absolute strings for robust comparison
    exclude_abs = set()
    try:
        if exclude_paths:
            for p in exclude_paths:
                if not p:
                    continue
                exclude_abs.add(os.path.abspath(str(p)))
    except Exception:
        exclude_abs = set()

    for entry in sorted(os.listdir(parent_dir)):
        sub = os.path.join(parent_dir, entry)
        # Skip any subdir that matches an excluded path (e.g. the pipeline output directory)
        try:
            if os.path.abspath(sub) in exclude_abs:
                logger.info(f"Skipping excluded directory during discovery: {sub}")
                continue
        except Exception:
            pass
        if not os.path.isdir(sub):
            continue
        # Attempt to find FASTA first
        fasta_found = []
        for ext in ('.fa', '.fasta', '.faa', '.fna'):
            # accept plain and gz/gzip-compressed variants
            for pattern in (f'*{ext}', f'*{ext}.gz', f'*{ext}.gzip'):
                for fn in glob.glob(os.path.join(sub, pattern)):
                    if os.path.isfile(fn):
                        fasta_found.append(fn)
        if fasta_found:
            samples.append((entry, fasta_found[0]))
            continue

        # Otherwise attempt to find paired fastq in subdir
        # Reuse discover_samples_from_input_dir on subdir (Paired-FASTQ mode)
        sub_samples = discover_samples_from_input_dir(sub, sequence_type, logger)
        if sub_samples:
            # pick first sample found
            name, spec = sub_samples[0]
            samples.append((entry, spec))
            continue

        logger.warning(f"No FASTA or paired FASTQ found in subdir: {sub}")
    # Return discovered samples (one per subdir)
    return samples


def combine_detection_matrices(output_root, sample_names, logger):
    """Combine per-sample detection_matrix files into per-database combined matrices.

    Produces two outputs per database found in sample outputs:

    1) <output_root>/<database>_combined_detection_matrix.tsv
       - Backwards-compatible binary presence/absence matrix: rows are genes, columns are samples (1/0), plus Total_Samples.

    2) <output_root>/<database>_combined_detection_matrix_tools.tsv
       - Informative matrix where each cell lists which tools detected the gene in that sample (comma-separated).
         Also includes Total_Samples (number of samples where gene was detected by at least one tool).

    This function reads each per-sample <database>_detection_matrix.tsv (generated by the per-sample pipeline)
    and infers which tools reported each gene for that sample by reading tool columns in the header.
    """
    import csv
    from pathlib import Path

    out_root = Path(output_root)
    db_sample_files = {}

    for s in sample_names:
        sdir = out_root / s
        if not sdir.is_dir():
            continue
        for p in sdir.glob('*_detection_matrix.tsv'):
            db = p.name.replace('_detection_matrix.tsv', '')
            db_sample_files.setdefault(db, {})[s] = p

    if not db_sample_files:
        logger.info("No per-sample detection_matrix files found to combine.")
        return

    created_files = []
    for db, sample_map in db_sample_files.items():
        # Collect all genes across samples and store per-sample detected tools
        all_genes = set()
        # sample -> gene -> list of tools
        sample_gene_tools = {s: {} for s in sample_map.keys()}

        for sample, path in sample_map.items():
            try:
                with open(path, 'r', newline='') as fh:
                    reader = csv.reader(fh, delimiter='\t')
                    header = next(reader, None)
                    if not header:
                        continue
                    # header expected: Gene, <tool1>, <tool2>, ..., Total_Detections (last)
                    # Identify tool columns (everything except 'Gene' and last 'Total_Detections')
                    gene_idx = 0
                    if len(header) < 2:
                        continue
                    total_idx = len(header) - 1
                    tool_cols = [(i, col) for i, col in enumerate(header) if i != gene_idx and i != total_idx]

                    for row in reader:
                        if not row:
                            continue
                        gene = row[gene_idx]
                        all_genes.add(gene)
                        detected_tools = []
                        for i, colname in tool_cols:
                            try:
                                val = row[i].strip()
                            except Exception:
                                val = ''
                            # Interpret '1' or non-zero/non-empty as detection
                            if val and val != '0':
                                detected_tools.append(colname)
                        sample_gene_tools.setdefault(sample, {})[gene] = detected_tools
            except Exception as e:
                logger.warning(f"Failed to read detection matrix for {sample} ({path}): {e}")

        if not all_genes:
            logger.info(f"No genes found in detection matrices for database '{db}'; skipping")
            continue

        genes_sorted = sorted(all_genes)

        # Prepare ordered sample list preserving discovery order where possible
        ordered_samples = [s for s in sample_names if s in sample_map]
        if not ordered_samples:
            ordered_samples = sorted(sample_map.keys())

        # 1) Write binary combined matrix (backwards-compatible)
        combined_path = out_root / f"{db}_combined_detection_matrix.tsv"
        try:
            with open(combined_path, 'w', newline='') as outfh:
                writer = csv.writer(outfh, delimiter='\t')
                header = ['Gene'] + ordered_samples + ['Total_Samples']
                writer.writerow(header)
                for gene in genes_sorted:
                    row = [gene]
                    total_samples = 0
                    for sample in ordered_samples:
                        tools = sample_gene_tools.get(sample, {}).get(gene, [])
                        present = 1 if tools else 0
                        row.append('1' if present else '0')
                        if present:
                            total_samples += 1
                    row.append(str(total_samples))
                    writer.writerow(row)
            logger.info(f"Combined detection matrix written: {combined_path}")
        except Exception as e:
            logger.warning(f"Failed to write combined detection matrix for {db}: {e}")

        # 2) Write informative tools-per-sample matrix
        tools_combined_path = out_root / f"{db}_combined_detection_matrix_tools.tsv"
        try:
            with open(tools_combined_path, 'w', newline='') as outfh:
                writer = csv.writer(outfh, delimiter='\t')
                header = ['Gene'] + ordered_samples + ['Total_Samples']
                writer.writerow(header)
                for gene in genes_sorted:
                    row = [gene]
                    total_samples = 0
                    for sample in ordered_samples:
                        tools = sample_gene_tools.get(sample, {}).get(gene, [])
                        if tools:
                            cell = '|'.join(sorted(tools))
                            total_samples += 1
                        else:
                            cell = ''
                        row.append(cell)
                    row.append(str(total_samples))
                    writer.writerow(row)
            logger.info(f"Informative tools-per-sample combined matrix written: {tools_combined_path}")
            created_files.extend([str(combined_path), str(tools_combined_path)])
        except Exception as e:
            logger.warning(f"Failed to write informative tools combined matrix for {db}: {e}")

    # Combine evidence-status matrices separately from binary detection calls.
    db_evidence_files = {}
    for sample in sample_names:
        sample_dir = out_root / sample
        if not sample_dir.is_dir():
            continue
        for path in sample_dir.glob('*_evidence_matrix.tsv'):
            database = path.name.replace('_evidence_matrix.tsv', '')
            db_evidence_files.setdefault(database, {})[sample] = path

    for database, sample_map in db_evidence_files.items():
        all_genes = set()
        sample_statuses = {sample: {} for sample in sample_map}
        sample_evidence = {sample: {} for sample in sample_map}
        for sample, path in sample_map.items():
            try:
                with open(path, 'r', newline='') as handle:
                    reader = csv.DictReader(handle, delimiter='\t')
                    for row in reader:
                        gene = row.get('Gene', '')
                        status = row.get('Best_Evidence_Status', 'NOT_DETECTED')
                        if gene:
                            all_genes.add(gene)
                            sample_statuses[sample][gene] = status
                            evidence_count = row.get(
                                'Evidence_Detections',
                                row.get('Evidence_Tools', ''),
                            )
                            try:
                                evidence_present = int(evidence_count) > 0
                            except (TypeError, ValueError):
                                evidence_present = status in {
                                    'EXACT_ALLELE_DETECTED',
                                    'CANDIDATE_ALLELE_DETECTED',
                                    'ALLELE_LIKE',
                                    'FAMILY_DETECTED',
                                    'MIXED_OR_MOSAIC',
                                    'PROFILE_DETECTED',
                                    'DETECTED',
                                    'LEGACY_RELAXED_DETECTED',
                                }
                            sample_evidence[sample][gene] = evidence_present
            except Exception as e:
                logger.warning(
                    f"Failed to read evidence matrix for {sample} ({path}): {e}"
                )

        ordered_samples = [sample for sample in sample_names if sample in sample_map]
        combined_evidence = out_root / f"{database}_combined_evidence_matrix.tsv"
        try:
            with open(combined_evidence, 'w', newline='') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow([
                    'Gene', *ordered_samples,
                    'Evidence_Samples', 'Candidate_Allele_Samples',
                    'Exact_Allele_Samples', 'Profile_Samples',
                    'Strict_Samples',
                ])
                for gene in sorted(all_genes):
                    statuses = [
                        sample_statuses.get(sample, {}).get(
                            gene, 'NOT_DETECTED'
                        )
                        for sample in ordered_samples
                    ]
                    evidence_samples = sum(
                        sample_evidence.get(sample, {}).get(gene, False)
                        for sample in ordered_samples
                    )
                    exact_samples = sum(
                        status == 'EXACT_ALLELE_DETECTED'
                        for status in statuses
                    )
                    candidate_samples = sum(
                        status == 'CANDIDATE_ALLELE_DETECTED'
                        for status in statuses
                    )
                    profile_samples = sum(
                        status == 'PROFILE_DETECTED'
                        for status in statuses
                    )
                    writer.writerow([
                        gene, *statuses,
                        evidence_samples, candidate_samples, exact_samples,
                        profile_samples, exact_samples + profile_samples,
                    ])
            created_files.append(str(combined_evidence))
            logger.info(
                f"Combined qualified evidence matrix written: {combined_evidence}"
            )
        except Exception as e:
            logger.warning(
                f"Failed to write combined evidence matrix for {database}: {e}"
            )

    # Write a small README/manifest into the output root describing the created combined files
    try:
        readme_path = out_root / 'combined_detection_matrices_README.txt'
        with open(readme_path, 'w') as rh:
            rh.write('Combined detection matrices generated by GeneFíor/GeneFíor-Combine\n')
            rh.write('\nFiles created:\n')
            if created_files:
                for p in created_files:
                    rh.write(f"  - {p}\n")
            else:
                rh.write('  (no combined matrices were created)\n')
            rh.write('\nFormats:\n')
            rh.write('  - <database>_combined_detection_matrix.tsv: binary presence/absence matrix. Columns: Gene, <sample1>, <sample2>, ..., Total_Samples. Cell=1 indicates the gene was detected in that sample by at least one tool.\n')
            rh.write("  - <database>_combined_detection_matrix_tools.tsv: informative matrix. Each sample cell lists the tool(s) that detected the gene (pipe-separated), or empty if none. Columns: Gene, <sample1>, <sample2>, ..., Total_Samples.\n")
            rh.write("  - <database>_combined_evidence_matrix.tsv: per-sample threshold-passing evidence status, including ambiguous, mosaic, candidate, and exact calls.\n")
            rh.write('\nIf you prefer Excel (.xlsx) outputs, consider running the combine tool and converting the TSVs to XLSX with your preferred tool.\n')
        logger.info(f"Wrote combined matrices README: {readme_path}")
    except Exception as e:
        logger.debug(f"Failed to write combined matrices README: {e}")
def handle_all_input_files(options, logger):
    # Process input files based on sequence type and tool requirements
    logger.info("=" * 70)
    logger.info("GeneFíor I/O processing...")
    logger.info("=" * 70)

    # Initialise cleanup tracking
    options.temp_files_to_cleanup = []

    # Handle Paired-FASTQ input
    if options.sequence_type == 'Paired-FASTQ':
        logger.info("Input type: Paired-FASTQ")
        options.input_fasta = None

        # Validate paired FASTQ files
        r1_path, r2_path = validate_paired_fastq(options, logger)
        options.input_fastq_original = (r1_path, r2_path)

        if getattr(options, 'fastp_trim', False):
            r1_processing, r2_processing = run_fastp_for_paired_reads(
                options, r1_path, r2_path, logger
            )
            options.temp_files_to_cleanup.extend([r1_processing, r2_processing])
            alignment_temp_dir = os.path.dirname(r1_processing)
        else:
            r1_processing, r2_processing = r1_path, r2_path
            alignment_temp_dir = options.temp_directory
            logger.info(
                "fastp trimming was not requested; non-trimmed reads were used "
                "for this paired FASTQ analysis"
            )

        # Prepare FASTQ files for alignment tools (add suffixes if needed)
        # Do this BEFORE FASTA conversion so both use the same modified files
        r1_prepared, r2_prepared, needs_cleanup = prepare_fastq_for_alignment(
            r1_processing, r2_processing,
            alignment_temp_dir,
            logger,
            getattr(options, 'force_modify_fastq', False)
        )

        # Track modified FASTQ files for cleanup
        if needs_cleanup:
            options.temp_files_to_cleanup.extend([r1_prepared, r2_prepared])

        # Store both original and prepared paths
        options.input_fastq = (r1_prepared, r2_prepared)

        # Check if BLAST-based tools require FASTA conversion
        requires_fasta = requires_fasta_conversion(options.tools)

        if requires_fasta:
            logger.info("BLAST-based tools detected, converting FASTQ to FASTA...")
            # Temporarily set options.input to point to prepared files for FASTA conversion
            original_input = options.input
            options.input = f"{r1_prepared},{r2_prepared}"

            FASTQ_to_FASTA(options, logger)

            # Restore original input
            options.input = original_input

            # Copy converted FASTA to temp directory if specified
            if options.temp_directory and hasattr(options, 'input_fasta') and options.input_fasta:
                logger.info("Copying converted FASTA to temp directory...")
                original_fasta = options.input_fasta
                options.input_fasta = copy_to_temp_directory(
                    options.input_fasta,
                    options.temp_directory,
                    logger
                )

                # Track for cleanup only if we actually copied it
                if options.input_fasta != original_fasta:
                    options.temp_files_to_cleanup.append(options.input_fasta)
                    logger.info(f"FASTA will be read from temp directory: {options.input_fasta}")

        else:
            # No FASTA conversion needed
            logger.info("No BLAST-based tools, keeping FASTQ format")
            options.input_fasta = None

        # # Prepare FASTQ files for alignment tools (add suffixes if needed)
        # r1_prepared, r2_prepared, needs_cleanup = prepare_fastq_for_alignment(
        #     r1_path, r2_path,
        #     options.temp_directory,
        #     logger
        # )
        #
        # # Track modified FASTQ files for cleanup
        # if needs_cleanup:
        #     options.temp_files_to_cleanup.extend([r1_prepared, r2_prepared])
        #
        # # Store both original and prepared paths
        # options.input_fastq = (r1_prepared, r2_prepared)/
        # options.input_fastq_original = (r1_path, r2_path)
        #
        # if not hasattr(options, 'input_fasta'):
        #     options.input_fasta = None
        #
        # # Keep FASTQ paths for reference
        # if not hasattr(options, 'input_fastq'):
        #     options.input_fastq = (r1_path, r2_path)
        # else:
        #     # No FASTA conversion needed
        #     logger.info("No BLAST-based tools, keeping FASTQ format")
        #     options.input_fastq = (r1_path, r2_path)
        #     options.input_fasta = None

    # Handle single FASTA or FASTQ input
    else:
        logger.info(f"Input type: {options.sequence_type}")
        if getattr(options, 'fastp_trim', False):
            logger.warning(
                "fastp trimming was requested, but --fastp-trim only applies to "
                "Paired-FASTQ input; the input will be used unchanged"
            )

        # Check input file exists and is a file
        # Support '-' as stdin. If '-' is provided and data is piped (not TTY), materialise into temp file.
        if options.input == '-':
            try:
                if not sys.stdin.isatty():
                    tmpf = tempfile.NamedTemporaryFile(delete=False, dir=options.temp_directory or None, prefix='genefior_stdin_', suffix='.fasta')
                    data = sys.stdin.buffer.read()
                    tmpf.write(data)
                    tmpf.flush()
                    tmpf.close()
                    logger.info(f"Materialised piped stdin to temporary FASTA: {tmpf.name}")
                    options.input = tmpf.name
                    options.temp_files_to_cleanup.append(tmpf.name)
                else:
                    logger.error("Input '-' was provided but stdin is a TTY (no piped data). Provide a file or pipe FASTA data into stdin.")
                    print("Error: Input '-' provided but no piped data detected", file=sys.stderr)
                    sys.exit(1)
            except Exception as e:
                logger.error(f"Failed to materialise stdin input '-': {e}")
                print(f"Error: Failed to read stdin: {e}", file=sys.stderr)
                sys.exit(1)

        if not os.path.isfile(options.input):
            logger.error(f"Input file '{options.input}' not found")
            print(f"Error: Input file '{options.input}' not found", file=sys.stderr)
            sys.exit(1)

        file_size_mb = os.path.getsize(options.input) / (1024 * 1024)
        logger.info(f"Input file: {options.input} ({file_size_mb:.2f} MB)")

        # Handle based on sequence type
        if options.sequence_type == 'Single-FASTA':
            #logger.info("Input type: Single-FASTA")
            # Copy FASTA to temp directory if specified
            if options.temp_directory:
                logger.info("Copying FASTA to temp directory...")
                original_input = options.input
                options.input_fasta = copy_to_temp_directory(
                    options.input,
                    options.temp_directory,
                    logger
                )

                # Track for cleanup only if it was actually copied
                if options.input_fasta != original_input:
                    options.temp_files_to_cleanup.append(options.input_fasta)
                    logger.info(f"FASTA will be read from temp directory: {options.input_fasta}")
            else:
                options.input_fasta = options.input

            options.input_fastq = None

        elif options.sequence_type == 'Genes-FASTA':
            # Genes-FASTA: treat like Single-FASTA (full-length gene sequences)
            if options.temp_directory:
                logger.info("Copying Genes-FASTA to temp directory...")
                original_input = options.input
                options.input_fasta = copy_to_temp_directory(
                    options.input,
                    options.temp_directory,
                    logger
                )
                if options.input_fasta != original_input:
                    options.temp_files_to_cleanup.append(options.input_fasta)
                    logger.info(f"Genes-FASTA will be read from temp directory: {options.input_fasta}")
            else:
                options.input_fasta = options.input
            options.input_fastq = None

        # elif options.sequence_type == 'Paired-FASTQ':
        #     # Check if BLAST-based tools require FASTA conversion
        #     requires_fasta = requires_fasta_conversion(options.tools)
        #
        #     if requires_fasta:
        #         logger.warning("FASTQ input with BLAST-based tools requires conversion")
        #         logger.warning("Consider using Paired-FASTQ type or pre-convert to FASTA")
        #         # Set for potential conversion
        #         options.input_fastq = options.input
        #         options.input_fasta = None
        #     else:
        #         options.input_fastq = options.input
        #         options.input_fasta = None

        else:
            logger.error(f"Unknown sequence type: {options.sequence_type}")
            print(f"Error: Unknown sequence type: {options.sequence_type}", file=sys.stderr)
            sys.exit(1)

    # Log final configuration
    #logger.info("=" * 70)
    #logger.info("IO configuration:")
    #logger.info(f"  FASTA input: {options.input_fasta if options.input_fasta else 'None'}")
    #logger.info(f"  FASTQ input: {options.input_fastq if options.input_fastq else 'None'}")
    logger.info(f"Output directory: {options.output}")

    if options.temp_directory:
        logger.info(f"Temp directory: {options.temp_directory}")
        logger.info(f"Files tracked for cleanup: {len(options.temp_files_to_cleanup)}")
        if options.temp_files_to_cleanup:
            for f in options.temp_files_to_cleanup:
                logger.info(f"  - {f}")

    logger.info("=" * 70)


def cleanup_all_temp_files(options, logger):
    if getattr(options, 'no_cleanup', False):
        logger.info("Temporary-file cleanup disabled by --no_cleanup")
        return
    if hasattr(options, 'temp_files_to_cleanup') and options.temp_files_to_cleanup:
        #logger.info("=" * 70)
        logger.info("Cleaning up temporary files...")
        logger.info("=" * 70)
        cleanup_temp_files(
            getattr(options, 'temp_directory', None),
            options.temp_files_to_cleanup,
            logger
        )


# Note: discovery helpers are available as module-level functions:
#   discover_samples_from_input_dir(...) and discover_samples_from_subdirs(...)
# Callers should import them directly from GeneFior.utils
