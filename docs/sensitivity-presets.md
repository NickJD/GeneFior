# Search sensitivity presets in plain language

GeneFíor and AMRfíor offer the same `--sensitivity` option. It selects a set of
search settings for the tools that are running. The preset can change how much
work a search tool does to find possible matches. More search effort can take
longer and may return more weak or ambiguous matches.

The preset does **not** change GeneFíor's identity, query-coverage, gene-coverage,
or evidence thresholds. Those answer a different question: after an alignment
has been found, does it meet the criteria for reporting evidence or calling a
gene? A more sensitive search can find additional alignments, but it does not
automatically make the detection criteria more permissive.

## Choosing a preset

| Preset | In everyday terms | Typical trade-off |
| --- | --- | --- |
| `default` | Use each tool's normal search settings. | Good starting point when you do not have a reason to change search effort. BLASTx retains GeneFíor's existing `blastx-fast` default. |
| `very-fast` | Spend less time looking for less obvious matches. | Fastest search-oriented choice; may miss more distant or less exact matches. |
| `fast` | Favor speed, with a little more search effort than `very-fast` in some tools. | Useful for a quick first pass. Some tools share settings with another preset. |
| `mid-sensitive` | Use an intermediate search setting where the tool provides one. | A middle-ground option; for some tools it is the same as the native default or the next named preset. |
| `sensitive` | Look harder for matches. | Can recover more divergent matches, usually with more run time and potentially more ambiguous hits. |
| `more-sensitive` | Increase search effort further where the tool supports a separate setting. | Often slower; the exact change depends on the tool. |
| `very-sensitive` | Use a high-effort search setting. | Can find more distant matches, with higher run time and more candidates that need careful interpretation. |
| `ultra-sensitive` | Use the highest-effort preset GeneFíor provides. | Usually the slowest choice. Some tools use custom high-effort settings rather than an official preset. |

The names are a convenient shared interface. They do **not** mean every tool
uses the same algorithm or that, for example, `sensitive` produces equivalent
results in BLAST, BWA, and HMMER. Some steps in the list intentionally use the
same settings for a particular tool.

## What the preset changes for each tool

### DIAMOND

DIAMOND receives its corresponding search flag for `very-fast`,
`mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, and
`ultra-sensitive`. The `fast` profile leaves DIAMOND at its native default;
that default sits between DIAMOND's `--fast` and `--mid-sensitive` settings.
`default` also leaves DIAMOND at its native default.

DIAMOND's `--more-sensitive` variation is not simply another uniformly
stronger step than `--sensitive`; it changes how motif soft-masking is used.
Treat these labels as practical presets, not a precise ranking of biological
recovery.

### Bowtie2

Bowtie2 has two alignment modes. GeneFíor uses end-to-end mode by default;
`--bowtie2-mode local` allows the ends of a read to be clipped when they do not
align. The sensitivity preset chooses a search preset within the selected mode.
For example:

```bash
genefior ... --sensitivity very-sensitive --bowtie2-mode local
```

uses Bowtie2's very-sensitive local preset. `mid-sensitive` and `sensitive`
currently share Bowtie2's sensitive preset. `more-sensitive`,
`very-sensitive`, and `ultra-sensitive` use its very-sensitive preset; the
ultra-sensitive profile adds extra search settings on top.

Changing to local mode can produce partial read alignments, so consider whether
that fits the question you are asking. GeneFíor's query-coverage and detection
criteria still apply separately.

### BLASTn, BLASTx, and BLASTp

BLAST uses task names rather than a single shared sensitivity scale. The
`very-fast` and `fast` profiles select its faster tasks. `mid-sensitive` and
higher profiles select the general tasks where available; those upper profiles
therefore do not represent an increasing sequence of BLAST-specific settings.
With `default`, BLASTn and BLASTp use their native defaults, while BLASTx keeps
GeneFíor's established `blastx-fast` behavior.

You can choose a task directly when you need more control:

```bash
genefior ... --blastn-task dc-megablast
genefior ... --blastx-task blastx
genefior ... --blastp-task blastp-short
```

An explicit task overrides the task selected by `--sensitivity` for that BLAST
program. Task availability and behavior depend on the BLAST program; a task
choice is not a change to GeneFíor's detection thresholds.

### BWA-MEM

BWA does not have named sensitivity presets like DIAMOND or Bowtie2. GeneFíor
uses heuristic seed and reseeding settings for the faster and higher-effort
profiles. The default, `mid-sensitive`, and `sensitive` profiles keep BWA's
normal settings. These BWA values are tuning suggestions, not upstream BWA
presets, and should be benchmarked against representative samples and a known
reference set before being used for production comparisons.

### HMMER

HMMER uses search filters to decide which candidate matches receive more
expensive analysis. The faster profiles make those filters stricter; the
sensitive profiles allow more candidates through. The ultra-sensitive profile
uses HMMER's `--max` option, which disables its acceleration filters and can be
substantially slower. These are proposed GeneFíor profiles, not named HMMER
presets.

HMMER reporting thresholds are separate. For example, `--hmmer-threshold-mode`
and E-value options control which reported matches are accepted; they are not
the same as the search filters controlled by `--sensitivity`.

### Minimap2

The shared sensitivity setting does not change Minimap2. Its `--minimap2-preset`
describes the input technology or mapping scenario, such as short reads
(`sr`), Oxford Nanopore (`map-ont`), or PacBio (`map-pb` / `map-hifi`). Choose
that option based on your reads.

## Practical starting points

- Start with `default` for routine runs.
- Try `sensitive` if you expect divergent sequences and can accept a longer run.
- Use `very-fast` or `fast` when you need a quick screen and can tolerate missing
  some less obvious matches.
- Reserve `ultra-sensitive` for cases where finding additional candidates is
  worth substantially more search time.
- Compare presets on a small, representative dataset with known expected genes
  before adopting a new setting for a large study.

Always review the result in its biological context. A hit is evidence for a
possible match; it is not by itself proof that the gene is present or functional.
