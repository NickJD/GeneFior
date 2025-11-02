import gzip
import os
import sys
import logging

logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

def FASTQ_to_FASTA(options):
    # If Paired-FASTQ, convert R1/R2 FASTQ -> FASTA and set options.input to the combined FASTA
    logger.info("FASTQ_to_FASTA: using seqtk to convert paired FASTQ -> combined FASTA")
    import subprocess
    import shutil

    def find_pair(input_spec):
        # Doesn't currently check if both are teh same R1/R2 files
        if ',' in input_spec:
            r1, r2 = map(str.strip, input_spec.split(',', 1))
            return r1, r2
        base = input_spec
        candidates = [base, base + '_R1.fastq', base + '_R1.fq', base + '_1.fastq', base + '_1.fq']
        r1_path = None
        for c in candidates:
            if os.path.exists(c):
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
    if not os.path.exists(r1_path) or not os.path.exists(r2_path):
        logger.error(f"Paired FASTQ files not found: {r1_path}, {r2_path}")
        sys.exit(1)

    conv_dir = os.path.join(options.output, 'paired_fastq_fasta')
    os.makedirs(conv_dir, exist_ok=True)
    combined_fasta = os.path.join(conv_dir, 'fastq_to_fasta_combined.fasta.gz')
    if os.path.exists(combined_fasta):
        logger.info(f"Found existing combined FASTA at `{combined_fasta}`; using it (skipping conversion)")
        options.fasta_input = combined_fasta
        options.fastq_input = (r1_path, r2_path)
        return

    # ensure seqtk is available
    if shutil.which('seqtk') is None:
        logger.error("`seqtk` not found in PATH. Install seqtk or provide a FASTA input.")
        sys.exit(1)

    def seqtk_fastq_to_fasta_stream(fastq_path, out_handle):
        cmd = ['seqtk', 'seq', '-A', fastq_path]
        logger.info(f"Running: {' '.join(cmd)}")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if proc.stdout is None:
            logger.error("Failed to start seqtk process")
            sys.exit(1)
        try:
            for chunk in iter(lambda: proc.stdout.read(8192), b''):
                if not chunk:
                    break
                out_handle.write(chunk)
        finally:
            proc.stdout.close()
            _, stderr = proc.communicate()
            if proc.returncode != 0:
                logger.error(f"seqtk failed: {stderr.decode(errors='ignore')}")
                sys.exit(1)

    # convert both FASTQ files and append into a single gzipped FASTA
    with gzip.open(combined_fasta, 'wb') as out:
        seqtk_fastq_to_fasta_stream(r1_path, out)
        seqtk_fastq_to_fasta_stream(r2_path, out)

    logger.info(f"Combined FASTA created at {combined_fasta}")
    options.fasta_input = combined_fasta
    options.fastq_input = (r1_path, r2_path)