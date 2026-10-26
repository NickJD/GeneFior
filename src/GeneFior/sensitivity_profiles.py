"""Tool-aware search sensitivity profiles for GeneFíor and AMRfíor.

The named profiles are practical presets, not claims that sensitivity is
identical across different alignment/search algorithms.

Bowtie2 mode and minimap2 read-technology presets are independent settings.
BWA and HMMER profiles are heuristic configurations that should be benchmarked
on representative data before production use.
"""

SENSITIVITY_LEVELS = (
    "default",
    "very-fast",
    "fast",
    "mid-sensitive",
    "sensitive",
    "more-sensitive",
    "very-sensitive",
    "ultra-sensitive",
)


# DIAMOND documents --fast, --mid-sensitive, --sensitive,
# --more-sensitive, --very-sensitive and --ultra-sensitive.
# No flag means use the DIAMOND native default.
DIAMOND_ARGS = {
    "default": [],
    "very-fast": ["--fast"],
    "fast": [],
    "mid-sensitive": ["--mid-sensitive"],
    "sensitive": ["--sensitive"],
    "more-sensitive": ["--more-sensitive"],
    "very-sensitive": ["--very-sensitive"],
    "ultra-sensitive": ["--ultra-sensitive"],
}


# Each Bowtie2 preset also specifies the appropriate alignment mode:
# end-to-end and local presets must not be mixed.
BOWTIE2_PRESETS = {
    "default": {
        "end-to-end": None,
        "local": None,
    },
    "very-fast": {
        "end-to-end": "--very-fast",
        "local": "--very-fast-local",
    },
    "fast": {
        "end-to-end": "--fast",
        "local": "--fast-local",
    },
    "mid-sensitive": {
        "end-to-end": "--sensitive",
        "local": "--sensitive-local",
    },
    "sensitive": {
        "end-to-end": "--sensitive",
        "local": "--sensitive-local",
    },
    "more-sensitive": {
        "end-to-end": "--very-sensitive",
        "local": "--very-sensitive-local",
    },
    "very-sensitive": {
        "end-to-end": "--very-sensitive",
        "local": "--very-sensitive-local",
    },
    "ultra-sensitive": {
        "end-to-end": "--very-sensitive",
        "local": "--very-sensitive-local",
    },
}

# Bowtie2 has no native "--ultra-sensitive" preset.
# This custom profile increases search effort beyond its very-sensitive preset.
# Benchmark this profile before using it for large production runs.
BOWTIE2_EXTRA_ARGS = {
    "ultra-sensitive": [
        "-N", "1",
        "-L", "18",
        "-D", "25",
        "-R", "4",
    ],
}


# BLAST exposes program-specific tasks rather than one shared sensitivity
# ladder. Higher tiers use the general task, where available; they do not
# pretend to provide progressively stronger BLAST-specific sensitivity levels.
BLAST_TASKS = {
    "default": {},

    "very-fast": {
        "blastn": "megablast",
        "blastx": "blastx-fast",
        "blastp": "blastp-fast",
    },
    "fast": {
        "blastn": "megablast",
        "blastx": "blastx-fast",
        "blastp": "blastp-fast",
    },

    "mid-sensitive": {
        "blastn": "blastn",
        "blastx": "blastx",
        "blastp": "blastp",
    },
    "sensitive": {
        "blastn": "blastn",
        "blastx": "blastx",
        "blastp": "blastp",
    },
    "more-sensitive": {
        "blastn": "blastn",
        "blastx": "blastx",
        "blastp": "blastp",
    },
    "very-sensitive": {
        "blastn": "blastn",
        "blastx": "blastx",
        "blastp": "blastp",
    },
    "ultra-sensitive": {
        "blastn": "blastn",
        "blastx": "blastx",
        "blastp": "blastp",
    },
}


# BWA-MEM has no equivalent named sensitivity presets.
# These are heuristic starting configurations using its seed length (-k)
# and reseeding ratio (-r). The native defaults are approximately -k 19 -r 1.5.
BWA_ARGS = {
    "default": [],

    "very-fast": [
        "-k", "25",
        "-r", "2.0",
    ],
    "fast": [
        "-k", "21",
        "-r", "1.75",
    ],

    # Native BWA-MEM defaults.
    "mid-sensitive": [],
    "sensitive": [],

    "more-sensitive": [
        "-k", "17",
        "-r", "1.25",
    ],
    "very-sensitive": [
        "-k", "15",
        "-r", "1.0",
    ],
    "ultra-sensitive": [
        "-k", "13",
        "-r", "1.0",
        "-c", "1000",
    ],
}


# HMMER acceleration-filter profiles.
# Higher F1/F2/F3 P-value thresholds let more targets proceed to later,
# more expensive search stages. --max disables heuristic filters.
#
# These are proposed profiles, not upstream HMMER named presets.
HMMER_ARGS = {
    "default": [],

    "very-fast": [
        "--F1", "0.005",
        "--F2", "0.0005",
        "--F3", "0.000005",
    ],
    "fast": [
        "--F1", "0.01",
        "--F2", "0.00075",
        "--F3", "0.0000075",
    ],

    # Use HMMER's native filter defaults.
    "mid-sensitive": [],
    "sensitive": [
        "--F1", "0.05",
        "--F2", "0.005",
        "--F3", "0.00005",
    ],
    "more-sensitive": [
        "--F1", "0.10",
        "--F2", "0.01",
        "--F3", "0.0001",
    ],
    "very-sensitive": [
        "--F1", "0.25",
        "--F2", "0.025",
        "--F3", "0.0005",
    ],
    "ultra-sensitive": [
        "--max",
    ],
}


def build_sensitivity_params(
    level="default",
    bowtie2_mode="end-to-end",
    blastn_task_override=None,
    blastx_task_override=None,
    blastp_task_override=None,
):
    """Build the structured sensitivity configuration passed to Workflow."""

    level = str(level or "default").strip().lower()
    bowtie2_mode = str(bowtie2_mode or "end-to-end").strip().lower()

    if level not in SENSITIVITY_LEVELS:
        raise ValueError(
            f"Unknown sensitivity profile '{level}'. "
            f"Choose from: {', '.join(SENSITIVITY_LEVELS)}"
        )

    if bowtie2_mode not in {"end-to-end", "local"}:
        raise ValueError(
            "Bowtie2 mode must be either 'end-to-end' or 'local'."
        )

    # Start with the BLAST task for the selected profile.
    blast_tasks = dict(BLAST_TASKS[level])

    # Explicit task settings take precedence over the global sensitivity profile.
    for tool, override in (
        ("blastn", blastn_task_override),
        ("blastx", blastx_task_override),
        ("blastp", blastp_task_override),
    ):
        if override:
            blast_tasks[tool] = str(override)

    return {
        "diamond": {
            "args": list(DIAMOND_ARGS[level]),
        },
        "bowtie2": {
            "mode": bowtie2_mode,
            "preset": BOWTIE2_PRESETS[level][bowtie2_mode],
            "args": list(BOWTIE2_EXTRA_ARGS.get(level, [])),
        },
        "blastn": {
            "task": blast_tasks.get("blastn"),
        },
        "blastx": {
            "task": blast_tasks.get("blastx"),
        },
        "blastp": {
            "task": blast_tasks.get("blastp"),
        },
        "bwa": {
            "args": list(BWA_ARGS[level]),
        },
        "hmmer": {
            "args": list(HMMER_ARGS[level]),
        },

        # Intentionally no generic minimap2 sensitivity flag.
        # Its -x preset should continue to describe the read technology.
    }