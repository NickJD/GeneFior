# AMRfíor (pronounced "feer", sounds like beer)
This toolkit utilises a combined approach that uses BLAST, BWA, Bowtie2, DIAMOND, and Minimap2 to search DNA and protein sequences against AMR databases (DNA and AA) such as CARD/RGI and ResFinder (future work will include more databases).



## Menu:

```commandline
    AMR Gene Detection Pipeline - Multi-tool alignment with detection matrices

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     Input FASTA file with sequences to analyse
  -o, --output OUTPUT   Output directory for results

Tool selection:
  --tools {blastn,blastp,diamond,bowtie2,bwa,minimap2} [{blastn,blastp,diamond,bowtie2,bwa,minimap2} ...]
                        Specify which tools to run (default: all)

Gene Detection Parameters:
  --min-cov, --min-coverage MIN_COVERAGE
                        Minimum coverage threshold in percent (default: 80.0)
  --min-id, --min-identity MIN_IDENTITY
                        Minimum identity threshold in percent (default: 80.0)
  --max_target_seqs MAX_TARGET_SEQS
                        Maximum number of "hits" to return per query sequence (default: 100)

Mode Selection:
  --dna-only            Run only DNA-based tools
  --protein-only        Run only protein-based tools
  --sensitivity {default,conservative,sensitive,very-sensitive}
                        Preset sensitivity levels - default means each tool uses its own default settings and very-
                        sensitive applies DIAMONDs --ultra-sensitive and Bowtie2s --very-sensitive-local presets

Tool-Specific Parameters:
  --minimap2-preset {sr,map-ont,map-pb,map-hifi}
                        Minimap2 preset: sr=short reads, map-ont=Oxford Nanopore, map-pb=PacBio, map-hifi=PacBio HiFi
                        (default: sr)

Runtime Parameters:
  -t, --threads THREADS
                        Number of threads to use (default: 4)

Examples:
  # Basic usage with default tools
  python amr_pipeline.py -i reads.fasta -o results/

  # Select specific tools
  python amr_pipeline.py -i reads.fasta -o results/ \
    --tools blastn diamond bowtie2

  # Custom thresholds with long reads
  python amr_pipeline.py -i nanopore.fasta -o results/ \
    -t 16 --min-cov 90 --min-id 85 \
    --minimap-preset map-ont --dna-only

```