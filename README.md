# GeneFíor (pronounced Gene "feer", sounds like beer)
This toolkit utilises a combined approach that uses BLAST, BWA, Bowtie2, DIAMOND, and Minimap2 to search DNA and protein sequences against DNA and AA sequence databases - Databases including CARD/RGI and ResFinder are installed with the bioconda package.

Please see GeneFior.pdf for a draft publication.

## Requirements: 
    - python >=3.10
    - samtools >=1.19.2
    - blast >=2.17.0
    - diamond >=2.1.13
    - bowtie2 >=2.5.4
    - bwa >=0.8.09
    - minimap2 >=2.30
    - mmseqs2 >=18-8cc5c
    - hmmer >=3.4
    - seqtk >=1.4
    - pigz >=2.8

### Installation:
GeneFíor is available via bioconda (https://anaconda.org/channels/bioconda/packages/genefior/overview). \
To install, use the following command:
```commandline
conda install -c bioconda genefior
``` 
If there are problems installing pigz, ensure conda-forge is added to your channels and try again: 

GeneFíor is also available via pip, but conda installation is recommended to ensure all dependencies are correctly installed.
```commandline
pip install genefior
```

## Menu for GeneFíor (GeneFíor or GeneFíor):
BLASTn and BLASTx are disabled by default due to their slow speed, but can be enabled if desired.
```commandline
GeneFíor v0.9.0 GeneFíor - The Multi-Tool Gene Detection Toolkit.

options:
  -h, --help            show this help message and exit

Required selection:
  -i INPUT, --input INPUT
                        Input FASTA/FASTAQ file(s) with sequences to analyse - Separate FASTQ R1 and R2 with a comma for Paired-FASTQ or single file path for Single-FASTA - .gz files accepted
  --input-dir INPUT_DIR
                        Path to a directory containing multiple sample FASTA files (for Single-FASTA) or multiple paired FASTQ files (for Paired-FASTQ).
  --input-subdirs INPUT_SUBDIRS
                        Path to a directory where each subdirectory contains one sample (a FASTA or paired FASTQ files).
  -st {Single-FASTA,Paired-FASTQ,Genes-FASTA}, --sequence-type {Single-FASTA,Paired-FASTQ,Genes-FASTA}
                        Specify the input Sequence Type: Single-FASTA, Paired-FASTQ (R1+R2) or Genes-FASTA. When Genes-FASTA is selected the pipeline treats the input as full-length gene FASTA(s) (DNA or protein)
                        and will skip read-mappers (bowtie2, bwa, minimap2).
  --db-path USER_DB_PATH
                        Path to the directory containing user-provided databases in correct format (see build_databases.sh) (can supply multiple paths separated by commas)
  -o OUTPUT, --output OUTPUT
                        Output directory for results

Output selection:
  --report-fasta {all,detected,detected-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected" will report the reads that passed detection thresholds for
                        each detected gene."detected-all" will report all reads for each detected gene. (default: None)

Tool selection:
  --tools {blastn,blastx,blastp,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,blastp,diamond,bowtie2,bwa,minimap2,all} ...]
                        Specify which tools to run - "all" will run all tools (default: all except blastx/p as they can be slow!!)

Query threshold Parameters:
  --q-min-cov QUERY_MIN_COVERAGE, --query-min-coverage QUERY_MIN_COVERAGE
                        Minimum coverage threshold in percent (HSP for blastx/n) (default: 40.0)
  --q-min-id QUERY_MIN_IDENTITY, --query-min-identity QUERY_MIN_IDENTITY
                        Minimum identity threshold in percent (HSP for blast/diamond) (default: 80.0)

Gene Detection Parameters:
  --d-min-cov DETECTION_MIN_COVERAGE, --detection-min-coverage DETECTION_MIN_COVERAGE
                        Minimum coverage threshold in percent (default: 80.0)
  --d-min-id DETECTION_MIN_IDENTITY, --detection-min-identity DETECTION_MIN_IDENTITY
                        Minimum identity threshold in percent (default: 80.0)
  --d-min-base-depth DETECTION_MIN_BASE_DEPTH, --detection-min-base-depth DETECTION_MIN_BASE_DEPTH
                        Minimum average base depth for detection - calculated against regions of the detected gene with at least one read hit (default: 1.0)
  --d-min-reads DETECTION_MIN_NUM_READS, --detection-min-num-reads DETECTION_MIN_NUM_READS
                        Minimum number of reads required for detection (default: 1)

Mode Selection:
  --dna-only            Run only DNA-based tools
  --protein-only        Run only protein-based tools
  --sensitivity {default,very-fast,fast,sensitive,more-sensitive,very-sensitive,ultra-sensitive}
                        Preset sensitivity levels. "default" leaves tools unchanged. Other presets map to common DIAMOND and Bowtie2 presets (e.g. "very-sensitive" -> DIAMOND --very-sensitive / Bowtie2 --very-
                        sensitive-local). Use --tool-param to override specific tool parameters if needed.

Tool-Specific Parameters:
  --minimap2-preset {sr,map-ont,map-pb,map-hifi}
                        Minimap2 preset: sr=short reads, map-ont=Oxford Nanopore, map-pb=PacBio, map-hifi=PacBio HiFi (default: sr)
  --blastx-task {blastx,blastx-fast}
                        Run blastx with task blastx-fast (default: blastx-fast)
  --genes-type {dna,aa}
                        (Required when -st Genes-FASTA) Specify whether provided genes FASTA contains DNA (dna) or amino-acid sequences (aa)
  --tool-param TOOL_PARAM
                        Specify extra parameters for a tool as TOOL="args". Can be provided multiple times. Example: --tool-param 'diamond=--more-sensitive -e 1e-5'
  -e EVALUE, --evalue EVALUE
                        E-value threshold to pass to BLASTn/x and DIAMOND (default: 10)

Runtime Parameters:
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 4)
  --chunk-jobs CHUNK_JOBS
                        Number of concurrent BLAST chunk jobs to run when chunking is active. If unset the pipeline auto-derives concurrency from total threads or defaults to 1
  --chunk-threads-per-job CHUNK_THREADS_PER_JOB
                        If set, reserve this many threads per chunk job; otherwise total threads are divided evenly across concurrent chunk jobs
  --preserve-chunks     Keep chunk files and per-chunk outputs after concatenation (useful for debugging)
  --max-fasta-chunk-mb MAX_FASTA_CHUNK_MB
                        Max FASTA chunk size in MiB (default: 200.0). Inputs larger than this will be chunked for per-chunk BLAST runs
  -tmp TEMP_DIRECTORY, --temp-directory TEMP_DIRECTORY
                        Path to temporary to place input FASTA/Q file(s) for faster IO during BLAST - Path will also be used for all temporary files (default: output directory)
  --force-modify-fastq  Force addition of /1 and /2 suffixes to paired FASTQ read IDs even if they appear unique
  --no_cleanup
  --verbose

Miscellaneous Parameters:
  -v, --version         Show program version and exit

Examples:
  # Basic usage with default tools (runs DNA & protein tools)
  GeneFior -i reads.fasta -st Single-FASTA --db-path <path-to-db> -o results

  # Select specific tools and output detected FASTA sequences
  GeneFior -i reads.fasta -st Single-FASTA --db-path <path-to-db> -o results     --tools diamond bowtie2     --report_fasta detected

  # Custom thresholds, paired-fastq input, threads and dna-only mode
  GeneFior -i reads_R1.fastq,reads_R2.fastq -st Paired-FASTQ --db-path <path-to-db>     -o results -t 16 --d-min-cov 90 --d-min-id 85     --dna-only```
```
# AMRfíor has been absorbed into GeneFíor but is still available as a separate command for backwards compatibility with the same functionality and AMR databases.
## Menu for AMRfíor:
CARD and resfinder databases are used by default, but user-provided databases can also be specified.
The NCBI AMR database is also available as an option.
All 3 databases are prepackaged and formatted as part of the bioconda installation of AMRfíor.

## Menu for AMRfíor (AMRfíor or AMRfíor):
BLASTn and BLASTx are disabled by default due to their slow speed, but can be enabled if desired.

```commandline
GeneFíor v0.9.0 - AMRfíor - The Multi-Tool AMR Gene Detection Toolkit.

options:
  -h, --help            show this help message and exit

Required selection:
  -i INPUT, --input INPUT
                        Input FASTA/FASTAQ file(s) with sequences to analyse - Separate FASTQ R1 and R2 with a comma for Paired-FASTQ or single file path for Single-FASTA - .gz files accepted
  --input-dir INPUT_DIR
                        Path to a directory containing multiple sample FASTA files (for Single-FASTA) or multiple paired FASTQ files (for Paired-FASTQ).
  --input-subdirs INPUT_SUBDIRS
                        Path to a directory where each subdirectory contains one sample (a FASTA or paired FASTQ files).
  -st {Single-FASTA,Paired-FASTQ,Genes-FASTA}, --sequence-type {Single-FASTA,Paired-FASTQ,Genes-FASTA}
                        Specify the input Sequence Type: Single-FASTA, Paired-FASTQ (R1+R2) or Genes-FASTA. When Genes-FASTA is selected the pipeline treats the input as full-length gene FASTA(s) (DNA or protein)
                        and will skip read-mappers (bowtie2, bwa, minimap2).
  -o OUTPUT, --output OUTPUT
                        Output directory for results

Output selection:
  --report-fasta {all,detected,detected-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected" will report the reads that passed detection thresholds for
                        each detected gene."detected-all" will report all reads for each detected gene. (default: None)

Tool selection:
  --tools {blastn,blastx,blastp,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,blastp,diamond,bowtie2,bwa,minimap2,all} ...]
                        Specify which tools to run - "all" will run all tools (default: all except blastx/p as they can be slow!!)

Database selection:
  --databases {resfinder,card,ncbi,user-provided} [{resfinder,card,ncbi,user-provided} ...]
                        Specify which AMR gene databases to use (default: resfinder and card) -If "user-provided" is selected, please ensure the path contains the appropriate databases set up as per the
                        documentation and specify the path with --user-db-path.
  --user-db-path USER_DB_PATH
                        Path to the directory containing user-provided databases (required if --databases includes "user-provided")

Query threshold Parameters:
  --q-min-cov QUERY_MIN_COVERAGE, --query-min-coverage QUERY_MIN_COVERAGE
                        Minimum coverage threshold in percent (HSP for blastx/n) (default: 40.0)
  --q-min-id QUERY_MIN_IDENTITY, --query-min-identity QUERY_MIN_IDENTITY
                        Minimum identity threshold in percent (HSP for blast/diamond) (default: 80.0)

Gene Detection Parameters:
  --d-min-cov DETECTION_MIN_COVERAGE, --detection-min-coverage DETECTION_MIN_COVERAGE
                        Minimum coverage threshold in percent (default: 80.0)
  --d-min-id DETECTION_MIN_IDENTITY, --detection-min-identity DETECTION_MIN_IDENTITY
                        Minimum identity threshold in percent (default: 80.0)
  --d-min-base-depth DETECTION_MIN_BASE_DEPTH, --detection-min-base-depth DETECTION_MIN_BASE_DEPTH
                        Minimum average base depth for detection - calculated against regions of the detected gene with at least one read hit (default: 1.0)
  --d-min-reads DETECTION_MIN_NUM_READS, --detection-min-num-reads DETECTION_MIN_NUM_READS
                        Minimum number of reads required for detection (default: 1)

Mode Selection:
  --dna-only            Run only DNA-based tools
  --protein-only        Run only protein-based tools
  --sensitivity {default,very-fast,fast,sensitive,more-sensitive,very-sensitive,ultra-sensitive}
                        Preset sensitivity levels. "default" leaves tools unchanged. Other presets map to common DIAMOND and Bowtie2 presets (e.g. "very-sensitive" -> DIAMOND --very-sensitive / Bowtie2 --very-
                        sensitive-local). Use --tool-param to override specific tool parameters if needed.

Tool-Specific Parameters:
  --minimap2-preset {sr,map-ont,map-pb,map-hifi}
                        Minimap2 preset: sr=short reads, map-ont=Oxford Nanopore, map-pb=PacBio, map-hifi=PacBio HiFi (default: sr)
  --blastx-task {blastx,blastx-fast}
                        Run blastx with task blastx-fast (default: blastx-fast)
  --genes-type {dna,aa}
                        (Required when -st Genes-FASTA) Specify whether provided genes FASTA contains DNA (dna) or amino-acid sequences (aa)
  --tool-param TOOL_PARAM
                        Specify extra parameters for a tool as TOOL="args". Can be provided multiple times. Example: --tool-param 'diamond=--more-sensitive -e 1e-5'
  -e EVALUE, --evalue EVALUE
                        E-value threshold to pass to BLASTn/x and DIAMOND (default: 10)

Runtime Parameters:
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 4)
  --chunk-jobs CHUNK_JOBS
                        Number of concurrent BLAST chunk jobs to run when chunking is active. If unset the pipeline auto-derives concurrency from total threads or defaults to 1
  --chunk-threads-per-job CHUNK_THREADS_PER_JOB
                        If set, reserve this many threads per chunk job; otherwise total threads are divided evenly across concurrent chunk jobs
  --preserve-chunks     Keep chunk files and per-chunk outputs after concatenation (useful for debugging)
  --max-fasta-chunk-mb MAX_FASTA_CHUNK_MB
                        Max FASTA chunk size in MiB (default: 200.0). Inputs larger than this will be chunked for per-chunk BLAST runs
  -tmp TEMP_DIRECTORY, --temp-directory TEMP_DIRECTORY
                        Path to temporary to place input FASTA/Q file(s) for faster IO during BLAST - Path will also be used for all temporary files (default: output directory)
  --force-modify-fastq  Force addition of /1 and /2 suffixes to paired FASTQ read IDs even if they appear unique
  --no_cleanup
  --verbose

Miscellaneous Parameters:
  -v, --version         Show program version and exit

Examples:
  # Basic usage with default tools (runs DNA & protein tools)
  AMRfíor -i reads.fasta -st Single-FASTA -o results

  # Select specific tools and output detected FASTA sequences
  AMRfíor -i reads.fasta -st Single-FASTA -o results     --tools diamond bowtie2     --report_fasta detected

  # Custom thresholds, paired-fastq input, threads and dna-only mode
  AMRfíor -i reads_R1.fastq,reads_R2.fastq -st Paired-FASTQ -o results/     -t 16 --d-min-cov 90 --d-min-id 85     --dna-only    
```
### Sensitivity presets map: 
```
Sensitivity presets provide convenient combinations of parameters for different use cases. The presets map to specific parameter settings for DIAMOND and Bowtie2 as follows:
| Preset           | DIAMOND Sensitivity Option  | Bowtie2 Sensitivity Option      |
|------------------|-----------------------------|---------------------------------|
| default          | No change (DIAMOND default) | No change (Bowtie2default)      |
| very-fast        | --faster                    | --very-fast-local               |
| fast             | --fast                      | --fast-local                    |
| sensitive        | --sensitive                 | --sensitive-local               |
| more-sensitive   | --more-sensitive            | --very-sensitive-local          |
| very-sensitive   | --very-sensitive            | --very-sensitive-local          |
| ultra-sensitive  | --ultra-sensitive           | --very-sensitive-local          | 
```
## Menu for GeneFíor-Recompute (or genefíor-recompute):

### GeneFíor-Recompute is used to recalculate detection statistics from existing sequence search outputs with different thresholds without needing to rerun the entire analysis.

```commandline
GeneFíor v0.9.0 - GeneFíor-Recompute: Recalculate detection statistics from existing sequence search outputs

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory containing Genefíor results (with raw_outputs/ subdirectory)
  -o OUTPUT, --output OUTPUT
                        Output directory for recomputed results
  --tools {blastn,blastx,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,diamond,bowtie2,bwa,minimap2,all} ...]
                        Specify which tools to recompute - "all" will recompute for all detected tools (default: all)

Query threshold Parameters:
  --q-min-cov QUERY_MIN_COVERAGE, --query-min-coverage QUERY_MIN_COVERAGE
                        Minimum coverage threshold in percent (default: 40.0)
  --q-min-id QUERY_MIN_IDENTITY, --query-min-identity QUERY_MIN_IDENTITY
                        Minimum identity threshold in percent (HSP for blast/diamond) (default: 80.0)

Gene Detection Parameters:
  --d-min-cov DETECTION_MIN_COVERAGE, --detection-min-coverage DETECTION_MIN_COVERAGE
                        Minimum coverage threshold in percent (default: 80.0)
  --d-min-id DETECTION_MIN_IDENTITY, --detection-min-identity DETECTION_MIN_IDENTITY
                        Minimum identity threshold in percent (default: 80.0)
  --d-min-base-depth DETECTION_MIN_BASE_DEPTH, --detection-min-base-depth DETECTION_MIN_BASE_DEPTH
                        Minimum average base depth for detection - calculated against regions of the detected gene with at least one read hit (default: 1.0)
  --d-min-reads DETECTION_MIN_NUM_READS, --detection-min-num-reads DETECTION_MIN_NUM_READS
                        Minimum number of reads required for detection (default: 1)

Output Parameterts:
  --report-fasta {None,all,detected,detected-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected" will report the reads that passed detection thresholds for
                        each detected gene."detected-all" will report all reads for each detected gene. (default: None)
  --query-fasta QUERY_FASTA
                        Specify the original query FASTA/FASTQ file used for alignment (required for reporting mapped sequences for BLAST/DIAMOND).

Miscellaneous Parameters:
  -v, --version         Show program version and exit

Examples:
  # Recompute with different thresholds
  GeneFior-recompute -i original_results/ -o recomputed_90_90/ \
    --d-min-cov 90 --d-min-id 90

  # More stringent depth requirement
  GeneFior-recompute -i original_results/ -o high_depth/ \
    --d-min-base-depth 5.0 --d-min-reads 10
```
## Menu for GeneFíor-Gene-Stats (or genefíor-gene-stats):

### GeneFíor-Gene-Stats is used to generate summary statistics and visualisations from Genefíor results.

```commandline
GeneFíor v0.9.0 - GeneFíor-Gene-Stats: Generate detailed coverage visualisations for searched genes

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory containing GeneFíor results
  -o OUTPUT, --output OUTPUT
                        Output directory for visualisation reports
  -g GENES, --genes GENES
                        Comma-separated gene names (FULL NAMES) or path to file with gene names (one per line)
  --all-genes           Include all genes found in raw outputs (default: only genes listed as detected in detection_matrix.tsv)
  --databases {resfinder,card,ncbi} [{resfinder,card,ncbi} ...]
                        Database(s) to interrogate
  --tools {blastn,blastx,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,diamond,bowtie2,bwa,minimap2,all} ...]
                        Tool(s) to interrogate
  --ref-fasta REF_FASTA
                        NOT IMPLEMENTED YET - Reference FASTA file for variant calling (optional)
  --query-fasta QUERY_FASTA
                        NOT IMPLEMENTED YET - Query FASTA file (your input reads) for BLAST base-level analysis (optional)
  --plot-per-tool       Generate individual per-tool coverage PNGs in addition to combined comparison plots (default: off)

Examples:
  # Visualise specific genes (FULL NAMES) from all tools
  GeneFior-gene-stats -i results/ -o vis/ \
    -g "sul1_2_U12338,tet(W)|ARO:3000194" \
    --databases resfinder card \
    --tools diamond bowtie2 bwa

  # Visualise from gene (FULL NAMES) list file with reference
  GeneFior-gene-stats -i results/ -o vis/ \
    -g genes_of_interest.txt \
    --databases resfinder \
    --tools blastn diamond 
```

## GeneFíor-Combine (genefior-combine)

GeneFíor-Combine is a small standalone tool for combining per-sample *_detection_matrix.tsv files produced by a
multi-sample GeneFíor/AMRfíor run into per-database combined matrices. It writes two outputs per database:

- <database>_combined_detection_matrix.tsv (binary presence/absence matrix compatible with previous behaviour)
- <database>_combined_detection_matrix_tools.tsv (informative matrix where each cell lists which tools detected the gene in that sample)

Usage:
```commandline
GeneFíor-Combine -i /path/to/output_root [--samples-file samples.txt] [--output /path/to/write] [--verbose]
```

```commandline
GeneFíor v0.9.0 - GeneFíor-Combine - Combine per-sample detection matrices into per-database combined matrices

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Root output directory containing per-sample result folders
  --samples-file SAMPLES_FILE
                        Optional file with one sample folder name per line to include (defaults to all subdirectories)
  --output OUTPUT       Output root directory where combined files will be written (defaults to --input)
  --verbose             Enable verbose/debug logging

Examples:
  GeneFior-Combine -i /path/to/output_dir
  GeneFior-Combine -i /path/to/output_dir --samples-file samples.txt
  ```

## Database Setup: See /src/Genefior/databases/ for details on setting up user-provided databases.
### GeneFíor includes an automated script in the Databases directory to automate the setup of user-provided databases.

--------------------------------------------------------------------------------

## Testing and CI (repository-provided test harness) - BETA!!!

This repository includes an integration test harness that runs the GeneFíor CLI against a small set of repository-provided test inputs. The harness is located at:

  `test_data/run_tests.sh`

### Test data layout

Place the canonical test inputs under `test_data/` with this layout:

- `test_data/Genes-AA/all_aa.gz` — protein gene FASTA (gzipped)
- `test_data/Genes-DNA/all.gz` — nucleotide gene FASTA (gzipped)
- `test_data/Paired-FASTQ/Test_R1.fastq.gz` and `Test_R2.fastq.gz`
- `test_data/Single-FASTA/*.fasta.gz`

### What the harness does

- Builds minimal temporary databases under `test_data/tmp_user_dbs/resfinder` when the local tools are available (runs `makeblastdb`, `diamond makedb`, `bowtie2-build`, `bwa index` where present). If binaries are missing the harness will still run GeneFíor but may skip index creation.
- Writes outputs into a centralised directory: `test_data/outputs/out_<testname>/`. If an `out_<testname>` already exists it is moved to `out_<testname>.prev` before the run so the harness can compare new vs previous results.
- After each successful run the harness computes SHA-256 checksums for stable files in the new and `.prev` outputs and compares them. By default transient files are excluded from checksums (`pipeline_*` logs, `run.log`, and `checksum_diff.txt`) so only stable data files are compared. When checksums differ a unified diff is written to `test_data/outputs/out_<testname>/checksum_diff.txt` and the test is marked as failed.
- The harness uses a relative path to the GeneFíor CLI (`../src/GeneFior/GeneFior.py` relative to `test_data/`) so it runs from within the repo without absolute paths. You can override the `PYTHON` environment variable to use a specific interpreter.

### Running tests locally

From the repository `test_data` directory run:

```bash
cd /path/to/GeneFior/test_data
./run_tests.sh
```

On success the script prints a per-test summary and exits `0`. If tests fail it exits non-zero and leaves logs and `checksum_diff.txt` in `test_data/outputs` for inspection.

### Compare-only mode

If you only want to compare existing `out_<test>` directories against their `.prev` counterparts (without re-running the pipeline), use:

```bash
./run_tests.sh --compare-only
# or
./run_tests.sh -c
```

This compares the current outputs against the `.prev` directory and reports any checksum differences.

### Prerequisites and CI

- The harness and the GeneFíor CLI call external bioinformatics tools (BLAST+, DIAMOND, Bowtie2, BWA, samtools). To run full integration tests locally install the toolchain (Conda + Bioconda is recommended). Example:

```bash
conda create -n genefior-test -c conda-forge -c bioconda python=3.11 blast diamond bowtie2 bwa samtools pigz -y
conda activate genefior-test
```

- The GitHub Actions CI included with this repository installs the toolchain on the runner so end-to-end tests run in CI.

### Customisation

- If you want a different set of files to be considered "stable outputs" (for example only `raw_outputs/` and `run_parameters.json`) I can update the harness to use a whitelist instead of the current exclusion rules.
- If you want log normalisation (strip timestamps) included in comparisons I can add a normalisation step before checksumming to reduce false positives.

If you'd like, I can also add a short `TESTING.md` to the repo and commit these harness instructions there. Tell me whether you'd like the current comparison rules changed or committed documentation added.
