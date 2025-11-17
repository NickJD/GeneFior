import os, sys
import gzip
import shutil
import subprocess


def prepare_fastq_for_alignment(r1_path, r2_path, temp_dir, logger):
    # Check if FASTQ read IDs need suffixes and create temporary modified files if needed.

    import tempfile

    # Check if suffixes are needed
    needs_suffix = check_read_id_uniqueness(r1_path, r2_path, logger, sample_size=10000)

    if not needs_suffix:
        logger.info("FASTQ read IDs are already unique, no modification needed")
        return r1_path, r2_path, False

    logger.warning(
        "FASTQ read IDs are not unique after first space. Creating modified FASTQ files with /1 and /2 suffixes...")

    # Create temp directory if needed
    if temp_dir:
        work_dir = temp_dir
    else:
        work_dir = tempfile.mkdtemp(prefix='amrfior_fastq_')

    os.makedirs(work_dir, exist_ok=True)

    # Create output paths
    r1_modified = os.path.join(work_dir,
                               f"{os.path.basename(r1_path).replace('.fastq', '').replace('.fq', '')}_R1_modified.fastq.gz")
    r2_modified = os.path.join(work_dir,
                               f"{os.path.basename(r2_path).replace('.fastq', '').replace('.fq', '')}_R2_modified.fastq.gz")

    def modify_fastq_ids(input_path, output_path, suffix):
        # dd suffix to FASTQ read IDs
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
    logger.info(f"Modifying R1: {r1_path} -> {r1_modified}")
    modify_fastq_ids(r1_path, r1_modified, '/1')

    logger.info(f"Modifying R2: {r2_path} -> {r2_modified}")
    modify_fastq_ids(r2_path, r2_modified, '/2')

    return r1_modified, r2_modified, True

def check_read_id_uniqueness(r1_path, r2_path, logger, sample_size=10000):
    # Check if read IDs remain unique after truncation at first space.
    from collections import defaultdict

    seen_ids = defaultdict(int)

    for fastq_path in [r1_path, r2_path]:
        cmd = ['seqtk', 'seq', '-A', fastq_path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        count = 0
        for line in proc.stdout:
            if line.startswith(b'>'):
                # Get base ID (portion before first space)
                read_id = line[1:].decode('utf-8', errors='ignore').strip()
                base_id = read_id.split()[0] if ' ' in read_id else read_id
                seen_ids[base_id] += 1
                count += 1
                if count >= sample_size:
                    break

        proc.kill()
        proc.wait()

    # If any ID appears more than once, we need suffixes
    needs_suffix = any(count > 1 for count in seen_ids.values())

    return needs_suffix

def requires_fasta_conversion(tools):
    return any(tool in ('blastn', 'blastx', 'diamond', 'all') for tool in (tools or []))

def FASTQ_to_FASTA(options, logger):
    # If Paired-FASTQ, convert R1/R2 FASTQ -> FASTA and set options.input to the combined FASTA
    logger.info("FASTQ_to_FASTA: using seqtk to convert paired FASTQ -> combined FASTA")

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
    if os.path.isfile(combined_fasta) and os.path.getsize(combined_fasta) > 0:
        logger.info(f"Found existing combined FASTA at `{combined_fasta}`; using it (skipping conversion)")
        options.input_fasta = combined_fasta
        options.input_fastq = (r1_path, r2_path)
        return

    # ensure seqtk is available
    if shutil.which('seqtk') is None:
        logger.error("`seqtk` not found in PATH. Install seqtk or provide a FASTA input.")
        sys.exit(1)

    def seqtk_fastq_to_fasta_stream(fastq_path, out_handle, suffix=None):
        cmd = ['seqtk', 'seq', '-A', fastq_path]
        logger.info(f"Running: {' '.join(cmd)}")
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except Exception as e:
            logger.error(f"Failed to start seqtk: {e}")
            sys.exit(1)

        if proc.stdout is None:
            logger.error("Failed to start seqtk process")
            sys.exit(1)
        try:
            if suffix is None:
                # No modification needed, stream directly
                for chunk in iter(lambda: proc.stdout.read(8192), b''):
                    if not chunk:
                        break
                    out_handle.write(chunk)
            else:
                # Need to modify headers, process line by line
                buffer = b''
                for chunk in iter(lambda: proc.stdout.read(8192), b''):
                    if not chunk:
                        break

                    buffer += chunk
                    lines = buffer.split(b'\n')
                    buffer = lines[-1]  # Keep incomplete line in buffer

                    for line in lines[:-1]:
                        if line.startswith(b'>'):
                            # Modify header: insert suffix after base ID (before first space)
                            read_id = line[1:].decode('utf-8', errors='ignore').rstrip()
                            if ' ' in read_id:
                                parts = read_id.split(' ', 1)
                                modified_id = f">{parts[0]}{suffix} {parts[1]}\n"
                            else:
                                modified_id = f">{read_id}{suffix}\n"
                            out_handle.write(modified_id.encode('utf-8'))
                        else:
                            # Sequence line, write as-is
                            out_handle.write(line + b'\n')

                # Write any remaining buffer
                if buffer:
                    if buffer.startswith(b'>'):
                        read_id = buffer[1:].decode('utf-8', errors='ignore').rstrip()
                        if ' ' in read_id:
                            parts = read_id.split(' ', 1)
                            modified_id = f">{parts[0]}{suffix} {parts[1]}\n"
                        else:
                            modified_id = f">{read_id}{suffix}\n"
                        out_handle.write(modified_id.encode('utf-8'))
                    else:
                        out_handle.write(buffer + b'\n')
        finally:
            try:
                proc.stdout.close()
            except Exception:
                pass
            _, stderr = proc.communicate()
            if proc.returncode != 0:
                stderr_text = stderr.decode(errors='ignore') if stderr else ''
                logger.error(f"seqtk failed: {stderr_text}")
                sys.exit(1)

    # Check if we need suffixes
    needs_suffix = check_read_id_uniqueness(r1_path, r2_path, logger)

    # Define suffixes to use
    r1_suffix = '/1' if needs_suffix else None
    r2_suffix = '/2' if needs_suffix else None

    if needs_suffix:
        logger.warning(
            f"Read IDs are not unique after first space. Appending {r1_suffix} and {r2_suffix} to read names.")

    # Convert both FASTQ files and append into a single FASTA
    with gzip.open(combined_fasta, 'wb') as out:
        seqtk_fastq_to_fasta_stream(r1_path, out, suffix=r1_suffix)
        seqtk_fastq_to_fasta_stream(r2_path, out, suffix=r2_suffix)

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
    for file_path in files_to_remove:
        if file_path and os.path.exists(file_path) and temp_dir in file_path:
            try:
                os.remove(file_path)
                logger.info(f"Removed temporary file: {file_path}")
                cleaned_count += 1
            except Exception as e:
                logger.warning(f"Failed to remove temporary file {file_path}: {e}")

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
    logger.info("AMRfíor I/O processing...")
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
            logger
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
