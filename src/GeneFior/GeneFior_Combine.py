import argparse
import sys
import os
import logging
from datetime import datetime
from pathlib import Path

try:
    from .utils import combine_detection_matrices
    from .constants import GENEFIOR_VERSION
except Exception:
    from utils import combine_detection_matrices
    from constants import GENEFIOR_VERSION


def setup_logger(output_dir: Path, verbose: bool = False):
    logger = logging.getLogger('GeneFior-Combine')
    if not getattr(logger, 'handlers', None):
        logger.setLevel(logging.DEBUG if verbose else logging.INFO)
        # Stream handler
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.DEBUG if verbose else logging.INFO)
        sh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(sh)
        # File handler
        fh_path = output_dir / f"combine_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        try:
            fh_path.parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(str(fh_path), mode='a')
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(fh)
        except Exception:
            logger.warning(f"Could not create file handler at {fh_path}")
    logger.propagate = False
    return logger


def discover_sample_dirs(root: Path):
    # Return sorted list of subdirectory names under root to consider as samples.
    names = []
    try:
        for entry in sorted(os.listdir(str(root))):
            p = root / entry
            if p.is_dir() and not entry.startswith('.'):
                names.append(entry)
    except Exception:
        pass
    return names


def main():
    parser = argparse.ArgumentParser(
        description='GeneFíor ' + GENEFIOR_VERSION + ' - GeneFíor-Combine - Combine per-sample detection matrices into per-database combined matrices',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  GeneFior-Combine -i /path/to/output_dir
  GeneFior-Combine -i /path/to/output_dir --samples-file samples.txt
''')

    parser.add_argument('-i', '--input', required=True, help='Root output directory containing per-sample result folders')
    parser.add_argument('--samples-file', required=False, help='Optional file with one sample folder name per line to include (defaults to all subdirectories)')
    parser.add_argument('--output', required=False, help='Output root directory where combined files will be written (defaults to --input)')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose/debug logging')

    options = parser.parse_args()

    input_dir = Path(options.input)
    if not input_dir.exists() or not input_dir.is_dir():
        print(f"Error: Input directory not found or not a directory: {input_dir}", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(options.output) if options.output else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logger(output_dir, verbose=bool(options.verbose))

    # Determine sample names
    sample_names = []
    if options.samples_file:
        sf = Path(options.samples_file)
        if not sf.exists() or not sf.is_file():
            logger.error(f"Samples file not found: {sf}")
            sys.exit(1)
        try:
            with open(sf, 'r') as fh:
                for line in fh:
                    ln = line.strip()
                    if ln and not ln.startswith('#'):
                        sample_names.append(ln)
        except Exception as e:
            logger.error(f"Failed to read samples file: {e}")
            sys.exit(1)
    else:
        sample_names = discover_sample_dirs(input_dir)

    if not sample_names:
        logger.error("No sample directories discovered. Nothing to combine.")
        sys.exit(1)

    logger.info(f"Found {len(sample_names)} sample(s) to include in combine: {', '.join(sample_names)}")

    try:
        combine_detection_matrices(str(output_dir), sample_names, logger)
        logger.info("Combine completed successfully")
        sys.exit(0)
    except Exception as e:
        logger.error(f"Combine failed: {e}")
        sys.exit(2)


if __name__ == '__main__':
    main()

