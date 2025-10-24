import gzip
import os
import sys
import logging

logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

def FASTQ_to_FASTA(options):
    # If Paired-FASTQ, convert R1/R2 FASTQ -> FASTA and set options.input to the combined FASTA
    logger.info("FASTQ_to_FASTA: starting paired FASTQ -> FASTA conversion")

    def open_maybe_gz(path):
        if path.endswith('.gz'):
            return gzip.open(path, 'rt')
        return open(path, 'r')

    def convert_fastq_to_fasta(fastq_path, fasta_path):
        logger.info(f"Converting FASTQ {fastq_path} -> FASTA {fasta_path}")


    if ',' in options.input:
        r1_path, r2_path = map(str.strip, options.input.split(',', 1))
    else:
        base = options.input
        candidates = [base, base + '_R1.fastq', base + '_R1.fq', base + '_1.fastq', base + '_1.fq']
        r1_path = None
        for c in candidates:
            if os.path.exists(c):
                r1_path = c
                break
        if not r1_path:
            logger.error("Could not locate R1 FASTQ. Provide `R1.fastq,R2.fastq` as `-i`.")
            sys.exit(1)
        # derive R2 from R1 with common patterns
        if '_R1.' in r1_path:
            r2_path = r1_path.replace('_R1.', '_R2.')
        elif '_1.' in r1_path:
            r2_path = r1_path.replace('_1.', '_2.')
        else:
            r2_path = r1_path.replace('_R1', '_R2')

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

    def _convert_to_combined(fastq_path, combined_out):
        logger.info(f"Converting FASTQ {fastq_path} -> combined FASTA {combined_fasta}")
        if fastq_path.endswith('.gz'):
            inf = gzip.open(fastq_path, 'rb')
        else:
            inf = open(fastq_path, 'rb')
        try:
            read = inf.readline
            write = combined_out.write
            while True:
                hdr = read()
                if not hdr:
                    break
                seq = read().rstrip(b'\r\n')
                read()  # plus line
                read()  # quality line
                if not seq:
                    break
                h = hdr.strip()
                if h.startswith(b'@'):
                    h = h[1:].split(b' ')[0]
                rec = b'>' + h + b'\n' + seq + b'\n'
                write(rec)
        finally:
            inf.close()
        logger.info(f"Finished converting {fastq_path}")

    with gzip.open(combined_fasta, 'wb') as combined_out:
        _convert_to_combined(r1_path, combined_out)
        _convert_to_combined(r2_path, combined_out)

    logger.info(f"Combined FASTA created at {combined_fasta}")

    options.fasta_input = combined_fasta
    options.fastq_input = (r1_path, r2_path)