# AMRfíor (pronounced "feer", sounds like beer)
This toolkit utilises a combined approach that uses BLAST, BWA, Bowtie2, DIAMOND, and Minimap2 to search DNA and protein sequences against AMR databases (DNA and AA) such as CARD/RGI and ResFinder.


## Menu:

```commandline
AMRfíor v0.1.1- The Multi-Tool AMR Gene Detection Workflow.

options:
  -h, --help            show this help message and exit

Required selection:
  -i, --input INPUT     Input FASTA/FASTAQ file(s) with sequences to analyse - separate R1 and R2 with a comma for Paired-FASTQ
  -st, --sequence-type {Single-FASTA,Paired-FASTQ}
                        Specify the input Sequence Type: Single-FASTA or Paired-FASTQ (R1+R2) - Willconvert paired-fastq to single fasta for BLAST and DIAMOND analyses
  -o, --output OUTPUT   Output directory for results

Output selection:
  --report_fasta {None,all,detected,detected-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected" will report the reads that passed detection
                        thresholds for each detected gene."detected-all" will report all reads for each detected gene. (default: None)

Tool selection:
  --tools {blastn,blastx,diamond,bowtie2,bwa,minimap2,all} [{blastn,blastx,diamond,bowtie2,bwa,minimap2,all} ...]
                        Specify which tools to run - "all" will run all tools (default: all except blastx as it is very slow)

Query threshold Parameters:
  --q-min-cov, --query-min-coverage QUERY_MIN_COVERAGE
                        Minimum coverage threshold in percent (default: 40.0)

Gene Detection Parameters:
  --d-min-cov, --detection-min-coverage DETECTION_MIN_COVERAGE
                        Minimum coverage threshold in percent (default: 80.0)
  --d-min-id, --detection-min-identity DETECTION_MIN_IDENTITY
                        Minimum identity threshold in percent (default: 80.0)

Mode Selection:
  --dna-only            Run only DNA-based tools
  --protein-only        Run only protein-based tools
  --sensitivity {default,conservative,sensitive,very-sensitive}
                        Preset sensitivity levels - default means each tool uses its own default settings and very-sensitive applies DIAMONDs --ultra-sensitive and Bowtie2s --very-sensitive-local
                        presets

Tool-Specific Parameters:
  --minimap2-preset {sr,map-ont,map-pb,map-hifi}
                        Minimap2 preset: sr=short reads, map-ont=Oxford Nanopore, map-pb=PacBio, map-hifi=PacBio HiFi (default: sr)

Runtime Parameters:
  -t, --threads THREADS
                        Number of threads to use (default: 4)
  --no_cleanup
  -v, --verbose

Examples:
  # Basic usage with default tools (runs DNA & protein tools)
  AMRfior -i reads.fasta -st Single-FASTA -o results/

  # Select specific tools and output detected FASTA sequences
  AMRfior -i reads.fasta -st Single-FASTA -o results/     --tools blastn diamond bowtie2     --report_fasta detected

  # Custom thresholds, paire-fastq input, threads and dna-only mode
  AMRfior -i reads_R1.fastq,reads_R2.fastq -st Paired-FASTQ -o results/     -t 16 --d-min-cov 90 --d-min-id 85     --dna-only

```