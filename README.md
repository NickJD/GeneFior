[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/genefior/README.html)

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
    - fastp >=1.3.5 (optional, for paired-read trimming)

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
GeneFíor v0.10.1 GeneFíor - The Multi-Tool Gene Detection Toolkit.

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
  --report-fasta {all,detected,detected-all,evidence,evidence-all,candidate,candidate-all,exact,exact-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected" will report the reads that passed detection thresholds for
                        each detected gene. Evidence modes are explicit aliases; candidate modes restrict output to candidate allele calls; exact modes restrict output to exact nucleotide allele calls. (default: None)

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
  --d-min-depth DETECTION_MIN_DEPTH, --detection-min-depth DETECTION_MIN_DEPTH
                        Minimum per-base depth required across the configured detection coverage in legacy-relaxed mode (default: 1)
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
  --fastp-trim, --fastp, --trim-reads
                        For Paired-FASTQ input, run fastp automatic adapter and quality trimming before analysis
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

## Menu for AMRfíor (AMRfíor or AMRFíor):
BLASTn and BLASTx are disabled by default due to their slow speed, but can be enabled if desired.

```commandline
GeneFíor v0.10.1 - AMRfíor - The Multi-Tool AMR Gene Detection Toolkit.

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
  --report-fasta {all,detected,detected-all,evidence,evidence-all,candidate,candidate-all,exact,exact-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected" will report the reads that passed detection thresholds for
                        each detected gene. Evidence modes are explicit aliases; candidate modes restrict output to candidate allele calls; exact modes restrict output to exact nucleotide allele calls. (default: None)

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
  --d-min-depth DETECTION_MIN_DEPTH, --detection-min-depth DETECTION_MIN_DEPTH
                        Minimum per-base depth required across the configured detection coverage in legacy-relaxed mode (default: 1)
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
  --fastp-trim, --fastp, --trim-reads
                        For Paired-FASTQ input, run fastp automatic adapter and quality trimming before analysis
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

### Optional paired-read trimming

Add `--fastp-trim` (also accepted as `--fastp` or `--trim-reads`) to a
`Paired-FASTQ` GeneFior or AMRfior run to apply fastp automatic adapter and
quality trimming before every selected search tool. The `fastp` executable
must be available in `PATH`. HTML and JSON QC reports are retained in
`<output>/fastp/`; the trimmed FASTQ files follow the normal temporary-file
cleanup policy. Without this option, the log explicitly records that
non-trimmed reads were used.

### Multiple custom databases with GeneFior

GeneFior accepts multiple independent custom databases through a
comma-separated `--db-path` value:

```bash
GeneFior \
  -i reads_R1.fastq.gz,reads_R2.fastq.gz \
  -st Paired-FASTQ \
  --db-path /data/databases/amr,/data/databases/virulence \
  --tools diamond blastn bowtie2 \
  -o results
```

Each directory must contain its own tool-specific subdirectories, such as
`diamond/`, `blast_dna/`, and `bowtie2/`. With multiple paths, GeneFior uses
each directory name as the database name in logs and output files. Directory
names must therefore be unique. A single path retains the historical
`user-provided-db` output label.

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
## Evidence qualification and allele resolution

GeneFíor no longer treats union coverage alone as proof that one exact database
allele is present. Every tool now receives a qualified evidence call:

- `EXACT_ALLELE_DETECTED`: literal nucleotide allele support with 100% mean identity, full-length coverage, and allele-resolving support.
- `CANDIDATE_ALLELE_DETECTED`: nucleotide allele candidate with full-length support, mean identity at least 98% by default, and allele-resolving support. For read inputs, every base must be covered at least 3x by default; for full-length `Genes-FASTA` DNA input, one full-length query sequence can support the candidate.
- `ALLELE_LIKE`: strong evidence close to the named allele, but candidate or exact identity is unresolved.
- `FAMILY_DETECTED`: the gene family is strongly supported, but reads do not distinguish its alleles.
- `PARTIAL_OR_DIVERGENT`: an incomplete or divergent mapping pattern that does
  not pass the user's evidence thresholds.
- `MIXED_OR_MOSAIC`: discontinuous or threshold-sensitive evidence inconsistent with one coherent allele.
- `PROFILE_DETECTED`: robust HMM profile evidence.
- `MUST_FLAG_REVIEW`: a must-flag profile passed the significance filter but
  failed ordinary evidence thresholds; it is flagged for review, not counted
  as positive evidence.
- `NOT_DETECTED`: insufficient evidence.

Protein searches (`BLASTx`, `BLASTp`, and DIAMOND) cannot by themselves claim
an exact nucleotide allele because synonymous differences are invisible.
Their strongest normal result is therefore `ALLELE_LIKE` or
`FAMILY_DETECTED`. Nucleotide read tools can make a candidate call when every
base is covered at least `--evidence-candidate-depth` times and the mean
identity is at least `--evidence-candidate-identity`. Full-length DNA
`Genes-FASTA` input uses full-length 1x allele coverage instead of read-depth
corroboration. Exact calls are stricter: candidate support plus literal 100%
identity.

Important output columns:

- `Evidence_Present=1` means the gene passes the user's configured detection
  coverage at the qualified evidence depth, identity, and minimum-read
  thresholds.
- `Candidate_Allele_Detected=1` means the named nucleotide allele is a strong
  candidate under the configured candidate-depth and candidate-identity thresholds
  for read inputs, or under full-length 1x coverage and candidate-identity
  thresholds for DNA `Genes-FASTA` input.
- `Exact_Allele_Detected=1` means the named nucleotide allele was resolved.
- `Profile_Detected=1` means a robust HMM profile was detected.
- `Detected` is retained as a backwards-compatible alias of
  `Evidence_Present`; it reports whether the user's detection thresholds pass.
- `Evidence_Status` and `Evidence_Warnings` should be used for interpretation.
- Qualified tool-stat files include `Detection_System`, `Evidence_Status`,
  `Evidence_Warnings`, `Evidence_Present`, `Candidate_Allele_Detected`,
  `Exact_Allele_Detected`, and the family/allele-resolution fields.
- Read export modes `detected`/`evidence` and
  `detected-all`/`evidence-all` export threshold-passing calls. The `candidate`
  and `candidate-all` modes restrict output to candidate allele calls, while
  `exact` and `exact-all` restrict output to exact nucleotide allele calls.

New per-database outputs:

- `<database>_detection_matrix.tsv`: user-threshold evidence binary calls.
- `<database>_evidence_matrix.tsv`: qualified per-tool evidence statuses.
- `<database>_evidence_summary.tsv`: evidence, candidate-allele, and exact-allele
  counts per tool and across all tools.
- `<database>_allele_resolution.tsv`: family-level top database candidate and competing alleles.

These qualified-only files and fields are produced only by
`--detection-system qualified`.

Detection-system compatibility:

- `--detection-system qualified` is the default. It reports threshold-passing
  evidence, distinguishes family/mosaic/allele-like/candidate-allele results,
  and applies the additional exact-allele requirements.
- `--detection-system legacy-relaxed` reproduces the original binary detector:
  gene coverage, mean identity, covered-region depth, and read count are
  compared directly with the user's thresholds. It does not require
  corroborated depth, continuous coverage, competitive uniqueness, or family
  resolution. For strict historical compatibility, its minimum-read check uses
  the original alignment/HSP count rather than deduplicated passing-read
  support. It never labels the result as an exact allele.
  Legacy tool-stat files report the quantitative mapping statistics followed
  by one binary `Detected` column. They do not emit evidence matrices,
  evidence summaries, allele-resolution reports, or qualified evidence fields.
  Gene-Stats likewise reports only the binary detection value for legacy runs.
- `--detection-mode` is accepted as an alias for `--detection-system`.

## Menu for GeneFíor-Recompute (or genefíor-recompute):

### GeneFíor-Recompute is used to recalculate detection statistics from existing sequence search outputs with different thresholds without needing to rerun the entire analysis.

```commandline
GeneFíor v0.10.1 - GeneFíor-Recompute: Recalculate detection statistics from existing sequence search outputs

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
  --d-min-depth DETECTION_MIN_DEPTH, --detection-min-depth DETECTION_MIN_DEPTH
                        Minimum per-base depth required across the configured detection coverage in legacy-relaxed mode (default: source run, otherwise 1)
  --d-min-reads DETECTION_MIN_NUM_READS, --detection-min-num-reads DETECTION_MIN_NUM_READS
                        Minimum number of reads required for detection (default: 1)

Output Parameters:
  --report-fasta {None,all,detected,detected-all,evidence,evidence-all,candidate,candidate-all,exact,exact-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected" will report the reads that passed detection thresholds for
                        each detected gene. Evidence modes are explicit aliases; candidate modes restrict output to candidate allele calls; exact modes restrict output to exact nucleotide allele calls. (default: None)
  --report-read-names {all,detected,detected-all,evidence,evidence-all,candidate,candidate-all,exact,exact-all}
                        Write deduplicated read-name files per database/tool/gene. Candidate modes restrict output to candidate allele calls; exact modes restrict output to exact nucleotide alleles. Read names are stored on disk during recompute rather than retained in memory.
  --query-fasta QUERY_FASTA
                        Specify the original query FASTA/FASTQ file used for alignment (required for FASTA reporting from BLAST/DIAMOND).

Recompute memory behaviour:
- By default, read names and per-hit identity values are not retained.
- Coverage/depth uses compact interval boundaries rather than a Python object per covered base.
- Read-name and FASTA reporting are opt-in and use a temporary disk-backed store.
- Query-level thresholds can be tightened safely. Relaxing them below the source run is rejected because source BLAST/DIAMOND outputs were already filtered; use `--allow-incomplete-relaxation` only when accepting incomplete results.

Miscellaneous Parameters:
  -v, --version         Show program version and exit

Examples:
  # Recompute with different thresholds
  GeneFior-recompute -i original_results/ -o recomputed_90_90/ \
    --d-min-cov 90 --d-min-id 90

  # More stringent depth requirement
  GeneFior-recompute -i original_results/ -o high_depth/ \
    --d-min-depth 5 --d-min-reads 10

  # Export threshold-passing read names and sequences for detected genes
  GeneFior-recompute -i original_results/ -o recomputed_with_reads/ \
    --tools blastx --report-read-names detected --report-fasta detected \
    --query-fasta reads.fasta.gz
```
## Menu for GeneFíor-Gene-Stats (or genefíor-gene-stats):

### GeneFíor-Gene-Stats is used to generate summary statistics and visualisations from Genefíor results.

```commandline
GeneFíor v0.10.1 - GeneFíor-Gene-Stats: Generate detailed coverage visualisations for searched genes

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory containing GeneFíor results
  -o OUTPUT, --output OUTPUT
                        Output directory for visualisation reports
  -g GENES, --genes GENES
                        Comma-separated gene names (FULL NAMES) or path to file with gene names (one per line)
  --all-genes           Include all genes found in raw outputs (default: genes with meaningful evidence in evidence_matrix.tsv, falling back to strict detections for older runs)
  --gene-selection {evidence,candidate,exact,candidate-or-exact,all}
                        Select genes automatically from qualified evidence outputs. Use candidate-or-exact to report every candidate or exact nucleotide allele call without specifying gene names.
  --databases {resfinder,card,ncbi} [{resfinder,card,ncbi} ...]
                        Database(s) to interrogate
  --tools {blastn,blastx,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,diamond,bowtie2,bwa,minimap2,all} ...]
                        Tool(s) to interrogate
  --ref-fasta REF_FASTA
                        Plain-text reference FASTA file for variant calling (optional)
  --query-fasta QUERY_FASTA
                        Original query FASTA/FASTQ file. Required for FASTA read export from BLAST/DIAMOND results.
  --raw-dir RAW_DIR     Directory containing raw GeneFior alignment outputs. May point directly to raw_outputs/ or to a result directory containing raw_outputs/. Useful when -i/--input is a recompute output.
  --plot-per-tool       Generate individual per-tool coverage PNGs in addition to combined comparison plots (default: off)
  --report-read-names   Write deduplicated read-name files for the selected genes, databases, and tools
  --report-fasta        Write deduplicated read FASTA files for the selected genes, databases, and tools
  --read-selection {all,passing}
                        Export all raw mappings or only mappings passing the original query identity/coverage thresholds (default: passing)

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

  # Export reads for one selected gene from selected tools
  GeneFior-gene-stats -i results/ -o gene_reads/ \
    -g "tet(Q)_1_L33696" \
    --databases resfinder \
    --tools blastx bowtie2 \
    --report-read-names --report-fasta --read-selection passing \
    --query-fasta reads.fasta.gz

  # Report all candidate/exact allele genes automatically
  GeneFior-gene-stats -i results/ -o candidate_exact_gene_stats/ \
    --gene-selection candidate-or-exact \
    --databases resfinder card \
    --tools blastn bwa minimap2

  # Generate Gene-Stats from a recompute output while using the original raw files
  GeneFior-gene-stats -i recomputed_results/ -o gene_stats/ \
    --raw-dir original_results/raw_outputs \
    --databases resfinder card \
    --tools blastn blastx diamond bowtie2 bwa minimap2 \
    --query-fasta reads_R1.fastq.gz,reads_R2.fastq.gz
```

## GeneFíor-Combine (genefior-combine)

GeneFíor-Combine is a small standalone tool for combining per-sample *_detection_matrix.tsv files produced by a
multi-sample GeneFíor/AMRfíor run into per-database combined matrices. It writes:

- <database>_combined_detection_matrix.tsv (binary presence/absence matrix compatible with previous behaviour)
- <database>_combined_detection_matrix_tools.tsv (informative matrix where each cell lists which tools detected the gene in that sample)
- <database>_combined_evidence_matrix.tsv (qualified per-sample evidence statuses with evidence, candidate, exact, profile, and strict sample counts)

Usage:
```commandline
GeneFíor-Combine -i /path/to/output_root [--samples-file samples.txt] [--output /path/to/write] [--verbose]
```

```commandline
GeneFíor v0.10.1 - GeneFíor-Combine - Combine per-sample detection matrices into per-database combined matrices

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
