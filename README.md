# AMRfíor (pronounced "feer", sounds like beer)
This toolkit utilises a combined approach that uses BLAST, BWA, Bowtie2, DIAMOND, and Minimap2 to search DNA and protein sequences against AMR databases (DNA and AA) such as CARD/RGI and ResFinder.


## Menu:

```commandline
AMRfíor - The AMR Multi-Tool Gene Detection Workflow.

options:
  -h, --help            show this help message and exit

Required selection:
  -i, --input INPUT     Input FASTA file with sequences to analyse
  -o, --output OUTPUT   Output directory for results

Output selection:
  --report_fasta {None,all,detected,detected-all}
                        Specify whether to output sequences that "mapped" to genes."all" should only be used for deep investigation/debugging."detected" will
                        report the reads that passed detection thresholds for each detected gene."detected-all" will report all reads for each detected gene.
                        (default: None)

Tool selection:
  --tools {blastn,blastx,diamond,bowtie2,bwa,minimap2} [{blastn,blastx,diamond,bowtie2,bwa,minimap2} ...]
                        Specify which tools to run (default: all)

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
                        Preset sensitivity levels - default means each tool uses its own default settings and very-sensitive applies DIAMONDs --ultra-
                        sensitive and Bowtie2s --very-sensitive-local presets

Tool-Specific Parameters:
  --minimap2-preset {sr,map-ont,map-pb,map-hifi}
                        Minimap2 preset: sr=short reads, map-ont=Oxford Nanopore, map-pb=PacBio, map-hifi=PacBio HiFi (default: sr)

Runtime Parameters:
  -t, --threads THREADS
                        Number of threads to use (default: 4)
  --no_cleanup
  -v, --verbose

Examples:
  # Basic usage with default tools
  python amr_pipeline.py -i reads.fasta -o results/

  # Select specific tools
  python amr_pipeline.py -i reads.fasta -o results/ \
    --tools blastn diamond bowtie2

  # Custom thresholds and dna-only mode
  python amr_pipeline.py -i nanopore.fasta -o results/ \
    -t 16 --d-min-cov 90 --d-min-id 85 \
     --dna-only
        

```