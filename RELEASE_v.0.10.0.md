# Qualified AMR Evidence and Workflow 

This release substantially strengthens GeneFíor and AMRFíor detection,
recomputation, read reporting, and gene-statistics generation. The central change is a clear
separation between simple detection (now named legacy) and qualified `biological'
evidence.

## Highlights

- Added a qualified evidence system that distinguishes exact alleles,
  candidate alleles, allele-like evidence, unresolved families, mixed or mosaic
  mappings, partial or divergent hits, profile detections, and review-only HMM
  calls.
- Preserved the historical detector through
  `--detection-system legacy-relaxed`.
- Reworked GeneFior-Recompute to process very large alignment outputs without
  retaining all hit identities or read names in memory.
- Added disk-backed read-name and read-sequence export to Recompute and
  Gene-Stats, including candidate-allele export modes.
- Corrected several coverage, identity, HMMER, BAM, and multi-HSP accounting
  errors that could exaggerate or misrepresent gene support.

## Detection And Allele Resolution

- Added the default `qualified` detection system.
- Added evidence statuses:
  - `EXACT_ALLELE_DETECTED`
  - `CANDIDATE_ALLELE_DETECTED`
  - `ALLELE_LIKE`
  - `FAMILY_DETECTED`
  - `MIXED_OR_MOSAIC`
  - `PARTIAL_OR_DIVERGENT`
  - `PROFILE_DETECTED`
  - `MUST_FLAG_REVIEW`
  - `NOT_DETECTED`
- Candidate nucleotide-allele calls require full-length support and 98% mean
  identity by default through `--evidence-candidate-identity`. Read inputs must
  cover every base at `--evidence-candidate-depth` (default min 3x); full-length DNA
  `Genes-FASTA` inputs use full-length 1x allele coverage instead.
- Exact nucleotide-allele calls now require literal 100% mean identity, 100%
  gene coverage, candidate-tier support, and sufficient competitive support.
- Protein searches can support a family or allele-like call but cannot claim
  a candidate or exact nucleotide allele.
- Added family inference, competing-allele reporting, unique-best support,
  ambiguous-best support, internal-gap checks, depth-CV warnings, and
  corroborated-coverage checks.
- Evidence remains tied to the user's coverage, identity, depth, and
  minimum-read thresholds. Partial and review states are visible diagnostics,
  not positive evidence.

## Legacy Compatibility

- Added `--detection-system legacy-relaxed` and the
  `--detection-mode` alias to GeneFior, AMRFior, and Recompute.
- Legacy mode applies the original direct coverage, identity, covered-depth,
  and sequence-count thresholds.
- Legacy outputs now remain deliberately simple:
  - tool-stat files end with one binary `Detected` field;
  - no evidence matrix, evidence summary, or allele-resolution report is
    generated;
  - Gene-Stats reports binary detection rather than qualified evidence labels.
- Qualified mode retains the complete evidence and allele-resolution fields,
  including `Candidate_Allele_Detected`.

## Coverage And Read-Support Corrections

- Multiple HSPs from the same read to the same gene are merged before depth,
  coverage, identity, and read-support accounting.
- Overlapping HSPs no longer inflate depth or sequence counts.
- Subject intervals are merged while preserving union coverage.
- Read support is counted once per read/gene combination.
- Secondary and supplementary BAM alignments are excluded.
- BAM alignments without an `NM` tag no longer silently receive 100% identity.
- MAPQ 255 is treated as unavailable mapping quality rather than unique
  support.
- BAM and BLAST query-coverage calculations now use consistent query-consuming
  coordinates.

## HMMER Corrections

- HMMER significance is now applied to each target/profile pair.
- A weak domain can no longer borrow significance from another target that hit
  the same profile.
- Only significant, quality-passing domains contribute to profile coverage and
  read support.
- Must-flag profiles that fail ordinary evidence thresholds are surfaced as
  review items rather than positive detections.

## Recompute

- Recompute now uses compact interval-based statistics for large BLAST and
  DIAMOND tables.
- Read names are retained only when explicitly requested.
- Optional read output uses a temporary SQLite-backed store.
- Added `--report-read-names` and expanded FASTA modes for all, evidence,
  detected, candidate, and exact calls.
- Added paired FASTA/FASTQ sequence recovery for exported reads.
- Recompute inherits source-run sequence type, gene type, query thresholds,
  detection thresholds, evidence settings, E-value, and bit-score settings
  when the user does not override them.
- Recompute manifests now record `source_input` and `source_raw_outputs` so
  Gene-Stats can locate the original large alignment files without copying them.
- DIAMOND blastx/blastp interpretation is inherited from source metadata where
  available.
- Relaxing query-level filters below those used by the source run is rejected
  by default because discarded hits cannot be recovered.
- Added `--allow-incomplete-relaxation` for explicit acceptance of incomplete
  relaxed results.
- Parser and samtools failures now terminate recomputation instead of yielding
  plausible-looking partial output.

## Gene-Stats

- Added `--gene-selection` so Gene-Stats can automatically report evidence,
  candidate, exact, candidate-or-exact, or all genes without requiring a
  hand-written gene list.
- Added `--raw-dir` so Gene-Stats can analyse recompute outputs while reading
  raw BLAST/DIAMOND/BAM files from the original run directory.
- Added read-name and FASTA export for selected genes and tools.
- Added support for paired query files.
- Added BLASTp and DIAMOND-blastp coordinate handling in amino-acid units.
- DIAMOND-blastx coordinates continue to be converted to nucleotide space for
  comparison with nucleotide mappers.
- Added arbitrary/custom database-name support and database auto-discovery.
- Excluded secondary and supplementary BAM alignments.
- Deduplicated overlapping HSP depth from the same read/gene pair.
- Parser and samtools failures now return a failing exit status.
- Coverage plots now include a dedicated low-depth panel so small non-zero
  coverage is visible instead of being flattened by high-depth peaks.
- Qualified runs display evidence status and warnings; legacy runs display
  only binary detection.

## Main Workflow And Input Handling

- Added optional paired-read preprocessing with
  `--fastp-trim`, `--fastp`, or `--trim-reads`.
- fastp performs automatic paired-end adapter and quality trimming and writes
  JSON and HTML reports.
- Runs without trimming explicitly log that untrimmed reads were used.
- Batch samples now receive isolated temporary directories, preventing
  basename collisions and cross-sample contamination.
- Batch and single-sample workflows now fail with a non-zero status when all
  selected tools fail.
- Parser errors are no longer converted into successful partial runs.
- Added multiple custom GeneFior databases through comma-separated database
  paths.
- Improved database discovery and tool availability reporting.
- DIAMOND logging now identifies `DIAMOND-BLASTX` and `DIAMOND-BLASTP`
  correctly.
- BLAST/DIAMOND result tables are consistently written as tab-separated files
  without mixed-format comment headers.

## New Outputs

Qualified runs produce:

- `<database>_detection_matrix.tsv`
- `<database>_evidence_matrix.tsv`
- `<database>_evidence_summary.tsv`
- `<database>_allele_resolution.tsv`
- expanded qualified tool-stat tables

Qualified summaries include candidate-allele counts separately from exact
alleles. Candidate calls are deliberately not exact calls.

Legacy runs produce:

- `<database>_detection_matrix.tsv`
- tool-stat tables with quantitative metrics and one binary `Detected` column

## Reliability And Tests

- Added regression tests for:
  - exact, candidate, and allele-like classification;
  - evidence versus partial diagnostic calls;
  - HMMER target/profile significance;
  - duplicate HSP accounting;
  - secondary BAM alignment exclusion;
  - DIAMOND blastp/blastx coordinate interpretation;
  - batch temporary-directory isolation;
  - source-run metadata inheritance;
  - paired read recovery;
  - fastp preprocessing;
  - custom database selection;
  - legacy and qualified output schemas.
- All four primary command-line entry points now render `--help` successfully.
- Current verification result: 166 tests passing.

## Compatibility Notes

- Qualified detection is the default and is intentionally more informative
  than the historical binary detector.
- Use `--detection-system legacy-relaxed` when a simple historical
  detected/not-detected result is required.
- Exact allele counts will generally be lower than broad mapping-based
  detections because exact calls now require literal sequence identity and
  candidate-tier support. Use `Candidate_Allele_Detected=1` for the middle tier
  between legacy detection and exact allele resolution.
- Recompute cannot recover alignments removed by stricter source-run query
  filters.
- GeneFior-Reconstruct was not part of these workflow changes.
