# GeneFíor (pronounced Gene "feer", sounds like beer)
This toolkit utilises a combined approach that uses BLAST, BWA, Bowtie2, DIAMOND, and Minimap2 to search DNA and protein sequences against DNA and AA sequence databases - Databases including CARD/RGI and ResFinder are preloaded.

## Requirements: 
    - python >=3.10
    - samtools >=1.19.2
    - blast >=2.17.0
    - diamond >=2.1.13
    - bowtie2 >=2.5.4
    - bwa >=0.7.19
    - minimap2 >=2.30
    - seqtk >=1.4
    - pigz >=2.8
### Installation:
GeneFíor is available via bioconda. To install, use the following command:
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
GeneFíor v0.7.0 GeneFíor - The Multi-Tool Gene Detection Toolkit.

options:
  -h, --help            show this help message and exit

Required selection:
  -i INPUT, --input INPUT
                        Input FASTA/FASTAQ file(s) with sequences to analyse - Separate FASTQ R1 and R2 with a comma for Paired-FASTQ or single file path
                        for Single-FASTA - .gz files accepted
  -st {Single-FASTA,Paired-FASTQ}, --sequence-type {Single-FASTA,Paired-FASTQ}
                        Specify the input Sequence Type: Single-FASTA or Paired-FASTQ (R1+R2) - Will convert Paired-FASTQ to single combined FASTA for BLAST
                        and DIAMOND analyses (SLOW)
  --db-path USER_DB_PATH
                        Path to the directory containing user-provided databases in correct format (see build_databases.sh) (can supply multiple paths
                        separated by commas)
  -o OUTPUT, --output OUTPUT
                        Output directory for results

Output selection:
  --report-fasta {all,detected,detected-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected"
                        will report the reads that passed detection thresholds for each detected gene."detected-all" will report all reads for each detected
                        gene. (default: None)

Tool selection:
  --tools {blastn,blastx,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,diamond,bowtie2,bwa,minimap2,all} ...]
                        Specify which tools to run - "all" will run all tools (default: all except blastx/n as it is very slow!!)

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
  --sensitivity {default,conservative,sensitive,very-sensitive}
                        Preset sensitivity levels - default means each tool uses its own default settings and very-sensitive applies DIAMONDs --ultra-
                        sensitive and Bowtie2s --very-sensitive-local presets

Tool-Specific Parameters:
  --minimap2-preset {sr,map-ont,map-pb,map-hifi}
                        Minimap2 preset: sr=short reads, map-ont=Oxford Nanopore, map-pb=PacBio, map-hifi=PacBio HiFi (default: sr)
  --blastx-task {blastx,blastx-fast}
                        Run blastx with task blastx-fast (default: blastx-fast)

Runtime Parameters:
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 4)
  --chunk-jobs CHUNK_JOBS
                        Number of concurrent BLAST chunk jobs to run when chunking is active. If unset the pipeline auto-derives concurrency from total
                        threads or defaults to 1
  --chunk-threads-per-job CHUNK_THREADS_PER_JOB
                        If set, reserve this many threads per chunk job; otherwise total threads are divided evenly across concurrent chunk jobs
  --preserve-chunks     Keep chunk files and per-chunk outputs after concatenation (useful for debugging)
  --max-fasta-chunk-mb MAX_FASTA_CHUNK_MB
                        Max FASTA chunk size in MiB (default: 200.0). Inputs larger than this will be chunked for per-chunk BLAST runs
  -tmp TEMP_DIRECTORY, --temp-directory TEMP_DIRECTORY
                        Path to temporary to place input FASTA/Q file(s) for faster IO during BLAST - Path will also be used for all temporary files
                        (default: output directory)
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
  GeneFior -i reads_R1.fastq,reads_R2.fastq -st Paired-FASTQ --db-path <path-to-db>     -o results -t 16 --d-min-cov 90 --d-min-id 85     --dna-only
```
# AMRFíor has been absorbed into GeneFíor but is still available as a separate command for backwards compatibility with the same functionality and AMR databases.
## Menu for AMRfíor:
CARD and resfinder databases are used by default, but user-provided databases can also be specified.
The NCBI AMR database is also available as an option.
All 3 databases are prepackaged and formatted as part of the bioconda installation of AMRfíor.

## Menu for AMRfíor (AMRfíor or AMRfíor):
BLASTn and BLASTx are disabled by default due to their slow speed, but can be enabled if desired.

```commandline
GeneFíor v0.7.0 AMRfíor - The Multi-Tool AMR Gene Detection Toolkit.

options:
  -h, --help            show this help message and exit

Required selection:
  -i INPUT, --input INPUT
                        Input FASTA/FASTAQ file(s) with sequences to analyse - Separate FASTQ R1 and R2 with a comma for Paired-FASTQ or single file path
                        for Single-FASTA - .gz files accepted
  -st {Single-FASTA,Paired-FASTQ}, --sequence-type {Single-FASTA,Paired-FASTQ}
                        Specify the input Sequence Type: Single-FASTA or Paired-FASTQ (R1+R2) - Will convert Paired-FASTQ to single combined FASTA for BLAST
                        and DIAMOND analyses (SLOW)
  -o OUTPUT, --output OUTPUT
                        Output directory for results

Output selection:
  --report-fasta {all,detected,detected-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected"
                        will report the reads that passed detection thresholds for each detected gene."detected-all" will report all reads for each detected
                        gene. (default: None)

Tool selection:
  --tools {blastn,blastx,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,diamond,bowtie2,bwa,minimap2,all} ...]
                        Specify which tools to run - "all" will run all tools (default: all except blastx/n as it is very slow!!)

Database selection:
  --databases {resfinder,card,ncbi,user-provided} [{resfinder,card,ncbi,user-provided} ...]
                        Specify which AMR gene databases to use (default: resfinder and card) -If "user-provided" is selected, please ensure the path
                        contains the appropriate databases set up as per the documentation and specify the path with --user-db-path.
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
  --sensitivity {default,conservative,sensitive,very-sensitive}
                        Preset sensitivity levels - default means each tool uses its own default settings and very-sensitive applies DIAMONDs --ultra-
                        sensitive and Bowtie2s --very-sensitive-local presets

Tool-Specific Parameters:
  --minimap2-preset {sr,map-ont,map-pb,map-hifi}
                        Minimap2 preset: sr=short reads, map-ont=Oxford Nanopore, map-pb=PacBio, map-hifi=PacBio HiFi (default: sr)
  --blastx-task {blastx,blastx-fast}
                        Run blastx with task blastx-fast (default: blastx-fast)

Runtime Parameters:
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 4)
  --chunk-jobs CHUNK_JOBS
                        Number of concurrent BLAST chunk jobs to run when chunking is active. If unset the pipeline auto-derives concurrency from total
                        threads or defaults to 1
  --chunk-threads-per-job CHUNK_THREADS_PER_JOB
                        If set, reserve this many threads per chunk job; otherwise total threads are divided evenly across concurrent chunk jobs
  --preserve-chunks     Keep chunk files and per-chunk outputs after concatenation (useful for debugging)
  --max-fasta-chunk-mb MAX_FASTA_CHUNK_MB
                        Max FASTA chunk size in MiB (default: 200.0). Inputs larger than this will be chunked for per-chunk BLAST runs
  -tmp TEMP_DIRECTORY, --temp-directory TEMP_DIRECTORY
                        Path to temporary to place input FASTA/Q file(s) for faster IO during BLAST - Path will also be used for all temporary files
                        (default: output directory)
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


## Menu for Genefíor-Recompute (Genefíor-Recompute or genefíor-recompute):

### Genefíor-Recompute is used to recalculate detection statistics from existing sequence search outputs with different thresholds without needing to rerun the entire analysis.

```commandline
eneFíor v0.6.0 - GeneFíor-Recompute: Recalculate detection statistics from existing sequence search outputs

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory containing Genefíor results (with
                        raw_outputs/ subdirectory)
  -o OUTPUT, --output OUTPUT
                        Output directory for recomputed results
  --tools {blastn,blastx,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,diamond,bowtie2,bwa,minimap2,all} ...]
                        Specify which tools to recompute - "all" will
                        recompute for all detected tools (default: all)

Query threshold Parameters:
  --q-min-cov QUERY_MIN_COVERAGE, --query-min-coverage QUERY_MIN_COVERAGE
                        Minimum coverage threshold in percent (default: 40.0)

Gene Detection Parameters:
  --d-min-cov DETECTION_MIN_COVERAGE, --detection-min-coverage DETECTION_MIN_COVERAGE
                        Minimum coverage threshold in percent (default: 80.0)
  --d-min-id DETECTION_MIN_IDENTITY, --detection-min-identity DETECTION_MIN_IDENTITY
                        Minimum identity threshold in percent (default: 80.0)
  --d-min-base-depth DETECTION_MIN_BASE_DEPTH, --detection-min-base-depth DETECTION_MIN_BASE_DEPTH
                        Minimum average base depth for detection - calculated
                        against regions of the detected gene with at least one
                        read hit (default: 1.0)
  --d-min-reads DETECTION_MIN_NUM_READS, --detection-min-num-reads DETECTION_MIN_NUM_READS
                        Minimum number of reads required for detection
                        (default: 1)

Output Parameterts:
  --report-fasta {None,all,detected,detected-all}
                        Specify whether to output sequences that "mapped" to
                        genes."all" should only be used for deep
                        investigation/debugging."detected" will report the
                        reads that passed detection thresholds for each
                        detected gene."detected-all" will report all reads for
                        each detected gene. (default: None)
  --query-fasta QUERY_FASTA
                        Specify the original query FASTA/FASTQ file used for
                        alignment (required for reporting mapped sequences for
                        BLAST/DIAMOND).

Miscellaneous Parameters:
  -v, --version         Show program version and exit

Examples:
  # Recompute with different thresholds
  Genefior-recompute -i original_results/ -o recomputed_90_90/ \
    --d-min-cov 90 --d-min-id 90

  # More stringent depth requirement
  Genefior-recompute -i original_results/ -o high_depth/ \
    --d-min-base-depth 5.0 --d-min-reads 10

```
## Menu for Genefíor-Gene-Stats (Genefíor-Gene-Stats or genefíor-gene-stats):

### Genefíor-Gene-Stats is used to generate summary statistics and visualizations from Genefíor results.

```commandline
Genefíor v0.6.0 - Genefíor-Gene-Stats: Generate detailed coverage visualisations for Gene genes

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory containing Genefíor results
  -o OUTPUT, --output OUTPUT
                        Output directory for visualisation reports
  -g GENES, --genes GENES
                        Comma-separated gene names (FULL NAMES) or path to file with gene names (one per line)
  --databases {resfinder,card,ncbi} [{resfinder,card,ncbi} ...]
                        Database(s) to interrogate
  --tools {blastn,blastx,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,diamond,bowtie2,bwa,minimap2,all} ...]
                        Tool(s) to interrogate
  --ref-fasta REF_FASTA
                        NOT IMPLEMENTED YET - Reference FASTA file for variant calling (optional)
  --query-fasta QUERY_FASTA
                        NOT IMPLEMENTED YET - Query FASTA file (your input reads) for BLAST base-level analysis (optional)

Examples:
  # Visualise specific genes (FULL NAMES) from all tools
  Genefior-gene-stats -i results/ -o vis/ \
    -g "sul1_2_U12338,tet(W)|ARO:3000194" \
    --databases resfinder card \
    --tools diamond bowtie2 bwa

  # Visualise from gene (FULL NAMES) list file with reference
  Genefior-gene-stats -i results/ -o vis/ \
    -g genes_of_interest.txt \
    --databases resfinder \
    --tools blastn diamond 

```

## Database Setup: See /src/Genefior/databases/ for details on setting up user-provided databases.
### Genefíor includes an automated script in the Databases directory to automate the setup of user-provided databases.