# GeneFíor v0.11.0

This release extends the qualified evidence workflow from gene and allele
detection to quantitative whole-genome mapping. It adds protein-level exact
matching, graded whole-genome completeness, multi-contig genome accounting,
multi-genome database support, richer visualisation, and continued reliability
and compatibility improvements from the `0.10.x` releases.

## Highlights

- Added `EXACT_PROTEIN_DETECTED` for exact protein-reference matches that do
  not prove an exact nucleotide allele.
- Added quantitative whole-genome mapping with coverage, depth, identity,
  read-support, and gap metrics.
- Added graded whole-genome statuses:
  - `WHOLE_GENOME_PARTIAL`
  - `WHOLE_GENOME_NEAR_COMPLETE`
  - `WHOLE_GENOME_COMPLETE`
- Added conservative multi-contig grouping for references named with explicit
  separators or recognised contig affixes.
- Added support for databases containing multiple genomes.
- Added whole-genome coverage files and visualisation support.
- Added combined heatmap generation and improved database building and
  discovery.
- Continued the qualified evidence, streaming Recompute, hAMRonization, and
  read-reporting work introduced in `v0.10.0`.

## Qualified Evidence

The default `qualified` detection system now distinguishes biological evidence
levels instead of reducing every result to a binary detected/not-detected
value.

Evidence statuses include:

- `EXACT_ALLELE_DETECTED`
- `EXACT_PROTEIN_DETECTED`
- `CANDIDATE_ALLELE_DETECTED`
- `ALLELE_LIKE`
- `FAMILY_DETECTED`
- `MIXED_OR_MOSAIC`
- `PARTIAL_OR_DIVERGENT`
- `PROFILE_DETECTED`
- `MUST_FLAG_REVIEW`
- `NOT_DETECTED`

Exact protein evidence requires 100% amino-acid identity and full-length
protein coverage. It is reported separately because synonymous DNA changes are
not visible at the protein level and therefore cannot establish an exact DNA
allele.

Nucleotide exact-allele calls remain stricter: they require literal 100%
identity, full-length coverage, candidate-tier support, corroborating depth,
and sufficient unique support. Candidate calls use the configured candidate
identity and depth thresholds.

The historical detector remains available with:

```text
--detection-system legacy-relaxed
```

## Whole-Genome Mapping

Whole-genome references are selected with:

```text
--db-whole-genome
```

These runs are mapping-oriented. The default tool set is `minimap2`, `bowtie2`,
and `bwa`; `blastn` can be included when explicitly requested. Protein search
and gene/allele resolution tools are not used for nucleotide whole-genome
references.

Whole-genome mapping reports quantitative evidence rather than a binary result.
The reports include:

- 1x, 2x, 3x, 5x, and 10x coverage;
- mean depth, median depth, and depth coefficient of variation;
- mean identity;
- mapped, passing, best, unique-best, and ambiguous-best read support;
- total gaps, internal gaps, and longest gap lengths;
- configured detection depth and coverage at that depth;
- final evidence status and warnings.

The genome-level statuses use these defaults:

- `WHOLE_GENOME_PARTIAL`: at least 20% of the combined genome is covered at
  1x, but near-complete quality or coverage requirements are not met.
- `WHOLE_GENOME_NEAR_COMPLETE`: at least 80% of the combined genome is covered
  at the detection depth, 3x by default, and configured identity and
  read-support thresholds pass.
- `WHOLE_GENOME_COMPLETE`: 100% of the combined genome is covered at 1x, at
  least 95% is covered at the detection depth, and no internal coverage gaps
  exceed the configured tolerance.

`WHOLE_GENOME_MAPPED` remains available for compatibility with older
per-reference results. New genome-level summaries use the graded statuses.

## Multi-Contig And Multi-Genome Databases

Whole-genome references are grouped conservatively into genomes. Supported
examples include:

- `genome_A|contig_1`;
- `genome_A::scaffold_2`;
- `genome_A_contig_1`;
- `genome_A_scaffold_2`;
- `genome_A_chr_1`.

Unrecognised reference names remain separate single-contig genomes rather than
being merged on a potentially unsafe shared prefix. This prevents unrelated
genomes from being treated as one genome.

Whole-genome runs now write:

- `<database>_contig_mapping_summary.tsv`, retaining one row per reference
  sequence;
- `<database>_genome_mapping_summary.tsv`, with one length-weighted row per
  grouped genome.

Missing contigs reduce combined genome coverage. Therefore, a well-covered
single contig cannot by itself produce a complete call for a multi-contig
genome. The contig report makes the underlying evidence and gaps auditable.

## Coverage And Mapping Improvements

- Added per-base sparse bedGraph generation from sorted BAM files, using
  `bedtools genomecov` where available and a `samtools depth` fallback.
- Added coverage and gap metrics to mapping statistics.
- Preserved contig identity in quantitative reports through `Genome_ID` and
  `Contig_ID`.
- Restricted whole-genome default execution to mapping-oriented tools.
- Recorded whole-genome mode in run parameters and generated reports.
- Continued the earlier corrections for merged HSPs, overlapping intervals,
  secondary/supplementary BAM alignments, unavailable MAPQ 255, missing `NM`
  tags, and consistent query-consuming coordinates.

## Visualisation And Database Tools

- Added whole-genome coverage plotting for per-sample mapping outputs.
- Added combined heatmap generation for multi-sample results.
- Improved database build handling for large Bowtie2 indexes, including
  `.bt2l` indexes.
- Improved database auto-discovery, naming, and tool availability reporting.
- Added arbitrary/custom database-name support.

## Recompute And Gene-Stats

- Recompute continues to support compact interval statistics for large BLAST
  and DIAMOND outputs.
- Added disk-backed read-name and sequence export modes.
- Added paired FASTA/FASTQ recovery and paired-query support.
- Recompute inherits source-run settings when overrides are not supplied.
- Added source manifest fields for locating original raw outputs.
- Added safeguards against incomplete relaxed recomputation.
- Parser and samtools failures now return failing status instead of plausible
  partial output.
- Gene-Stats supports automatic gene selection, raw-output directories,
  qualified evidence fields, candidate/exact export modes, and per-tool
  coverage plots with a low-depth panel.

## hAMRonization And Combination

- Added hAMRonization export with required database-version validation.
- Added configurable minimum call tiers for harmonized output.
- Added protein-exact evidence ranking without treating it as exact nucleotide
  evidence.
- Combined detection outputs now retain qualified evidence, candidate, exact,
  profile, and protein-exact sample counts.
- Fixed the GeneFior-Combine output-directory handling.

## Input Handling And Reliability

- Added optional paired-read preprocessing with `fastp`.
- Added isolated batch-sample temporary directories.
- Batch and single-sample workflows return non-zero status when all selected
  tools fail.
- Improved DIAMOND tool naming and BLAST/DIAMOND tabular output consistency.
- Added an `--always-flag-genes` workflow for review-priority targets.

## Outputs

Qualified gene-centric runs produce:

- `<database>_detection_matrix.tsv`;
- `<database>_evidence_matrix.tsv`;
- `<database>_evidence_summary.tsv`;
- `<database>_allele_resolution.tsv`;
- expanded qualified tool-stat tables.

Whole-genome runs additionally produce:

- `<database>_contig_mapping_summary.tsv`;
- `<database>_genome_mapping_summary.tsv`;
- per-tool sparse coverage bedGraph files where BAM mapping output is
  available.

Legacy-relaxed runs continue to produce binary detection outputs without the
qualified evidence matrix, evidence summary, or allele-resolution report.

## Compatibility Notes

- Qualified detection remains the default and is intentionally more informative
  than the historical binary detector.
- Use `--detection-system legacy-relaxed` when a historical binary result is
  required.
- Exact nucleotide allele counts can be lower than broad mapping detections.
- `EXACT_PROTEIN_DETECTED` must not be interpreted as an exact DNA allele.
- Whole-genome completeness is calculated across all grouped contigs, not from
  the best contig alone.
- If a database mixes gene references and whole-genome references, separate
  runs are recommended because `--db-whole-genome` applies to the selected
  database set as a whole.
- Recompute cannot recover alignments removed by stricter source-run query
  filters.
- GeneFior-Reconstruct remains outside the core whole-genome evidence workflow.

## Verification

The current repository verification includes regression coverage for:

- qualified exact, candidate, protein-exact, and whole-genome statuses;
- multi-contig aggregation and multi-genome separation;
- coverage, depth, identity, read-support, and gap calculations;
- HMMER target/profile significance;
- duplicate HSP and secondary BAM accounting;
- Recompute streaming and metadata inheritance;
- paired read recovery and fastp preprocessing;
- custom database selection and output schemas.

The v0.11.0 whole-genome regression tests pass alongside the existing
qualified evidence and hAMRonization tests.
