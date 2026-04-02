import os, sys
import gzip
import shutil
import json
import shlex
import time
import tempfile
import subprocess

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
    return any(tool in ('blastn', 'blastx', 'diamond', 'all') for tool in (tools or []))

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

        # Prepare FASTQ files for alignment tools (add suffixes if needed)
        # Do this BEFORE FASTA conversion so both use the same modified files
        r1_prepared, r2_prepared, needs_cleanup = prepare_fastq_for_alignment(
            r1_path, r2_path,
            options.temp_directory,
            logger,
            getattr(options, 'force_modify_fastq', False)
        )

        # Track modified FASTQ files for cleanup
        if needs_cleanup:
            options.temp_files_to_cleanup.extend([r1_prepared, r2_prepared])

        # Store both original and prepared paths
        options.input_fastq = (r1_prepared, r2_prepared)
        options.input_fastq_original = (r1_path, r2_path)

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
    if hasattr(options, 'temp_files_to_cleanup') and options.temp_files_to_cleanup:
        #logger.info("=" * 70)
        logger.info("Cleaning up temporary files...")
        logger.info("=" * 70)
        cleanup_temp_files(
            getattr(options, 'temp_directory', None),
            options.temp_files_to_cleanup,
            logger
        )
