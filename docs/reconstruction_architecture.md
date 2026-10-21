# Reconstruction evaluation workflow

`GeneFior-Reconstruct` builds consensus sequences from read alignments and
`GeneFior-Architect` adjudicates those sequences against the sample's assembled
genomes. Architect reports comparisons to both the long-read assembly and the
enhanced-contig assembly when available. A missing or fragmented truth locus is
reported separately from a failed reconstruction. Architect uses dense short
seed fallback when sparse exact seeds miss and searches each reconstruction
across every assembly contig independently of the expected reference locus.
A `SEARCH_UNRESOLVED` result is explicitly a search limitation, not evidence
that the gene is absent.

Install the optional alignment dependency with `pip install '.[architect]'`.

## Example

```sh
GeneFior-Reconstruct \
  -output_dir /data/FAN_AMR/genefior_out/232758E_18 \
  -gene 'catA1_1_V00622,tet(B)_2_AF326777' \
  -reference_fasta src/GeneFior/databases/resfinder/all.gz \
  -recon_dir outputs/FAN_AMR/reconstruction/232758E_18

GeneFior-Architect \
  --comparison-mode assembly \
  --reconstruction-dir outputs/FAN_AMR/reconstruction \
  --truth-root /data/FAN_AMR \
  --reference-fasta src/GeneFior/databases/resfinder/all.gz \
  --amrfior-root /data/FAN_AMR/genefior_out \
  --output-dir outputs/FAN_AMR/architect
```

Architect has two explicit comparison modes. `assembly` (the default) compares
to matched long-read and enhanced-contig assemblies; it requires `--truth-root`
and `--reference-fasta` (the latter identifies the expected gene locus in each
assembly). `database` compares directly to the named gene sequence in the
original nucleotide FASTA supplied with `--database-fasta`, and does not
require assemblies or `--truth-root`. Pass the FASTA file itself (for example,
`src/GeneFior/databases/resfinder/all.gz`), not its containing directory:

```sh
GeneFior-Architect \
  --comparison-mode database \
  --reconstruction-dir outputs/FAN_AMR/reconstruction \
  --database-fasta src/GeneFior/databases/resfinder/all.gz \
  --output-dir outputs/FAN_AMR/architect/database
```

Keep the output directories separate when retaining reports from both modes.
Every comparison row and JSON summary records its mode and target.
Database mode reports both the fraction of the reconstructed query aligned and
the fraction of the database gene covered (`database_reference_coverage_pct`).
Extra flanking sequence can lower query coverage while the full reference gene
is still covered. A database-only comparison cannot validate those flanks, so
they are called out in the explanation and should be checked in assembly mode
when an assembled genome is available.

Reconstruction files should follow the layout
`<reconstruction-dir>/<sample-id>/<gene>/<mapper>/<gene>_consensus.fasta`.
Architect discovers completed files on each run, so rerunning it includes
newly generated sample directories. It writes a tab-separated table and a JSON
report with the matched assembly, identity, aligned fraction, ambiguous-base
rate, status, explanation, and suggested next action.
When `--amrfior-root` is given, it also writes `architect_candidates.tsv` from
completed per-sample evidence matrices. Incomplete sample runs are skipped and
appear in the candidate list after AMRFior writes their evidence matrix.
Candidates are ranked by resolved allele or protein evidence first, followed
by supported allele, profile, mixed/mosaic, and broader family evidence.

Generate the companion figures after Architect has refreshed its TSV report:

```sh
GeneFior-Reconstruction-Visualise \
  --architect-report outputs/FAN_AMR/architect/architect_evaluation.tsv \
  --reconstruction-dir outputs/FAN_AMR/reconstruction \
  --reference-fasta src/GeneFior/databases/resfinder/all.gz \
  --amrfior-root /data/FAN_AMR/genefior_out \
  --truth-root /data/FAN_AMR \
  --output-dir outputs/FAN_AMR/visuals
```

The figures include SNP positions on reference coordinates, AMRFior candidate
allele labels and competing references alongside read/assembly evidence, and
reconstruction spans mapped to the assemblies. Assembly locus counts require
at least 90% identity and 80% reference coverage; these count matching loci,
not cellular copy dosage. Competing database labels are not treated as
coexisting sample alleles.

The reconstruction pileup uses primary alignments with MAPQ ≥10 and base
quality ≥20, excludes duplicate-flagged alignments, and counts overlapping
mates once per fragment. If overlapping mates disagree, the higher-quality
call is retained only when it is at least 10 Phred points stronger; otherwise
the fragment is excluded at that base. Haplotype assignment combines mate
calls under the shared template name. The SNP evidence TSV reports each
mapper's base call, depth, alternate-allele fraction, and whether mappers
disagree at that position. The visualisation
manifest records the input paths, file sizes, modification times, scope, and
output names. Re-run Architect and then the visualiser after AMRFior finishes
or adds a sample; both rediscover current files on each run, and the manifest
provides a snapshot of which inputs were included. The manifest is a run
snapshot, not a content hash or immutable archive.

AMRFior candidate labels and the number of competing reference labels describe
database resolution, not counts of physical alleles. Informative sites and
emitted haplotypes show sequence diversity in the read evidence, while
reference-like assembly hits identify candidate genome locations. Current
reports deliberately say copy number is not inferred. A validated allele
dosage estimate would need read-level phasing/support and repeat-aware locus
resolution, calibrated against mixtures with known allele composition.

Current boundary handling preserves the full reference interval in the
consensus, leaving uncalled reference positions as `N`. It searches for
soft-clip-supported extension in a window of 20% of reference length at each
end by default, rounded up to a codon boundary; terminal `N`s are stripped from
the extension wings only. It does not trim the gene to inferred start/stop
codons: ORF checks report warnings but do not alter the sequence. Use the
assembly-mapped boundary figure to spot when captured flanks extend beyond the
gene or when the aligned consensus starts inside the reference boundary.

## Statuses

| Status | Interpretation |
| --- | --- |
| `PASS` | At least 95% of the reconstructed sequence aligns at 99% or higher identity. |
| `LIKELY_PASS` | Strong match, with a small boundary or identity discrepancy. |
| `REVIEW` | Partial, ambiguous, or meaningfully divergent sequence. |
| `FAILED` | Assembly has a reference-like locus, but the reconstruction has poor identity or coverage to it. |
| `SEARCH_UNRESOLVED` | Sensitive seeded search found no reference-like locus; this does not establish absence. |
| `UNASSESSABLE` | Missing reference, assembly, or insufficiently confident truth-locus alignment. |

These thresholds are initial defaults for this dataset and should be reviewed
against the broader AMRFior cohort before being treated as calibrated cutoffs.

## Regression checks

Run `PYTHONPATH=src pytest tests/test_reconstruction_evaluation.py` for controlled checks of
codon-aware extension sizing, forward/reverse assembly alignments, divergent
chromosome-sized targets, fragment-level overlap and quality filtering,
truncated matches, and distinct repeated loci. These are unit cases; they do not replace
validation on known mixture samples or manual review of ambiguous repeats.

GeneFior-Reconstruct v2.2.1 removes prior per-gene, per-mapper files before rebuilding them.
This prevents optional artifacts from a previous run (for example, haplotypes
or variants) from surviving when the new run no longer supports them.
