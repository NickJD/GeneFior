import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_sample_data(sample_path: Path, database: str):
    matrix_file = sample_path / f"{database}_detection_matrix.tsv"

    if not matrix_file.exists():
        return None

    try:
        df = pd.read_csv(matrix_file, sep="\t")

        if df.empty or df.shape[1] == 0:
            return None

        # Normalise column names while retaining the original gene values.
        original_columns = list(df.columns)
        normalised = [str(c).strip().lower() for c in df.columns]
        df.columns = normalised

        # Identify gene column.
        if "gene" in df.columns:
            gene_col = "gene"
        else:
            # Fall back to first column, matching the baseline behaviour.
            gene_col = df.columns[0]

        # Identify value/count column.
        possible_count_cols = [
            "count",
            "reads",
            "depth",
            "abundance",
            "status",
            "detected",
            "detection",
            "value",
        ]

        value_col = next(
            (col for col in possible_count_cols if col in df.columns),
            None,
        )

        if value_col is None:
            if df.shape[1] < 2:
                print(
                    f"Warning: {matrix_file} contains no value column; skipping.",
                    file=sys.stderr,
                )
                return None
            value_col = df.columns[1]

        result = df[[gene_col, value_col]].copy()
        result.columns = ["gene", "value"]

        # Remove missing/empty gene identifiers.
        result["gene"] = result["gene"].astype(str).str.strip()
        result = result[
            result["gene"].notna()
            & (result["gene"] != "")
            & (result["gene"].str.lower() != "nan")
        ]

        if result.empty:
            return None

        # Detection matrices may contain numeric counts or binary values.
        # Convert numeric values directly. For textual status fields, common
        # positive labels are interpreted as detected.
        numeric = pd.to_numeric(result["value"], errors="coerce")

        if numeric.notna().any():
            result["value"] = numeric.fillna(0.0)
        else:
            positive_labels = {
                "1",
                "true",
                "yes",
                "y",
                "detected",
                "present",
                "positive",
                "hit",
            }
            result["value"] = (
                result["value"]
                .astype(str)
                .str.strip()
                .str.lower()
                .isin(positive_labels)
                .astype(float)
            )

        # If a gene occurs more than once within a sample, retain the
        # maximum value. This avoids duplicate pivot entries and preserves
        # detection if any hit is present.
        result = (
            result.groupby("gene", as_index=False)["value"]
            .max()
        )

        return result

    except Exception as exc:
        print(
            f"Warning: Could not read {matrix_file}: {exc}",
            file=sys.stderr,
        )
        return None


def find_sample_folders(input_path: Path):
    """Return sample directories in deterministic order."""
    return sorted(
        (p for p in input_path.iterdir() if p.is_dir()),
        key=lambda p: p.name,
    )


# ---------------------------------------------------------------------------
# Matrix construction
# ---------------------------------------------------------------------------

def build_combined_matrix(input_dir: str, database: str):
    """
    Load all sample detection matrices and return a combined numeric matrix.

    Rows are genes and columns are samples.
    """
    input_path = Path(input_dir)

    if not input_path.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_path}")

    if not input_path.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_path}")

    all_data = []
    sample_folders = find_sample_folders(input_path)

    print(f"Searching for samples in {input_path}...")

    for sample_folder in sample_folders:
        df = load_sample_data(sample_folder, database)

        if df is None or df.empty:
            continue

        df["sample_id"] = sample_folder.name
        all_data.append(df)

    if not all_data:
        raise RuntimeError(
            f"No data found for database '{database}' in {input_path}. "
            f"Expected files named '{database}_detection_matrix.tsv' "
            "inside sample directories."
        )

    combined = pd.concat(all_data, ignore_index=True)

    # If a gene/sample combination occurs more than once, retain the maximum.
    heatmap_df = combined.pivot_table(
        index="gene",
        columns="sample_id",
        values="value",
        aggfunc="max",
        fill_value=0,
    )

    # Deterministic sample order.
    heatmap_df = heatmap_df.reindex(
        sorted(heatmap_df.columns),
        axis=1,
    )

    # Remove genes with no detection in any sample.
    heatmap_df = heatmap_df.loc[
        heatmap_df.sum(axis=1) > 0
    ]

    if heatmap_df.empty:
        raise RuntimeError(
            f"No detected genes remain for database '{database}'."
        )

    return heatmap_df


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def prepare_detection_matrix(heatmap_df: pd.DataFrame):
    """
    Convert numeric abundance/count values into a binary detection matrix.

    Any value > 0 is considered detected.
    """
    return (heatmap_df.astype(float) > 0).astype(int)


def sort_detection_matrix(
    detection_df: pd.DataFrame,
    sort_method: str,
):
    """
    Sort genes according to the requested ordering.

    prevalence:
        Most widely detected genes first.

    alphabetical:
        Alphabetical gene order.

    none:
        Preserve the input matrix order.
    """
    if sort_method == "prevalence":
        prevalence = detection_df.sum(axis=1)
        return detection_df.assign(
            __prevalence=prevalence
        ).sort_values(
            ["__prevalence"],
            ascending=False,
            kind="stable",
        ).drop(columns="__prevalence")

    if sort_method == "alphabetical":
        return detection_df.sort_index()

    return detection_df


def make_heatmap(
    detection_df: pd.DataFrame,
    database: str,
    output_dir: Path,
    dpi: int = 300,
    fmt: str = "both",
    cmap: str = "Greys",
    cell_lines: bool = True,
    max_width: float = 40.0,
    max_height: float = 30.0,
):
    """
    Create the publication-style detection heatmap with a prevalence panel.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    n_genes = len(detection_df.index)
    n_samples = len(detection_df.columns)

    # Size scales with matrix dimensions but is capped to avoid enormous
    # figures for very large matrices.
    fig_width = min(max(10.0, 0.18 * n_samples + 6.0), max_width)
    fig_height = min(max(8.0, 0.095 * n_genes + 3.5), max_height)

    fig = plt.figure(
        figsize=(fig_width, fig_height)
    )

    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=[8.8, 1.7],
        wspace=0.06,
        left=0.22,
        right=0.96,
        top=0.95,
        bottom=0.065,
    )

    ax = fig.add_subplot(grid[0])

    # Binary detection: 0 = not detected, 1 = detected.
    image = ax.imshow(
        detection_df.values,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=1,
    )

    # Compact sample labels.
    sample_labels = [
        str(sample).replace(f"_{database}_genefior/", "")
        for sample in detection_df.columns
    ]
    sample_labels = [
        label.replace("_genefior/", "")
        for label in sample_labels
    ]

    ax.set_xticks(np.arange(n_samples))
    ax.set_xticklabels(
        sample_labels,
        rotation=90,
        fontsize=6,
        ha="center",
    )

    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels(
        detection_df.index,
        fontsize=6.2,
    )

    ax.tick_params(axis="both", length=0)

    ax.set_xlabel(
        "Sample",
        fontsize=10,
        labelpad=8,
    )
    ax.set_ylabel(
        "AMR gene / target",
        fontsize=10,
        labelpad=8,
    )

    ax.set_title(
        f"{database.upper()} antimicrobial resistance gene detection",
        loc="left",
        fontsize=15,
        fontweight="bold",
        pad=14,
    )

    # Optional fine cell boundaries.
    if cell_lines:
        ax.set_xticks(
            np.arange(-0.5, n_samples, 1),
            minor=True,
        )
        ax.set_yticks(
            np.arange(-0.5, n_genes, 1),
            minor=True,
        )
        ax.grid(
            which="minor",
            linewidth=0.18,
            alpha=0.35,
        )
        ax.tick_params(
            which="minor",
            bottom=False,
            left=False,
        )

    # Prevalence panel.
    prevalence = detection_df.sum(axis=1)

    ax2 = fig.add_subplot(grid[1], sharey=ax)

    ax2.barh(
        np.arange(n_genes),
        prevalence.values,
        height=0.72,
    )

    ax2.invert_yaxis()
    ax2.set_xlim(0, n_samples)
    ax2.set_xlabel(
        "Samples\ndetected",
        fontsize=9,
    )
    ax2.set_title(
        "Prevalence",
        fontsize=10,
        pad=14,
    )

    ax2.tick_params(
        axis="y",
        left=False,
        labelleft=False,
    )
    ax2.tick_params(
        axis="x",
        labelsize=8,
    )

    ax2.grid(
        axis="x",
        alpha=0.25,
        linewidth=0.6,
    )

    # Summary.
    total_detections = int(detection_df.values.sum())

    fig.text(
        0.22,
        0.018,
        (
            f"{n_genes} AMR targets across {n_samples} samples • "
            f"{total_detections:,} gene–sample detections • "
            "Targets ordered by prevalence."
        ),
        fontsize=8.5,
        ha="left",
    )

    # Keep the binary nature explicit in the colourbar.
    cbar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.027,
        pad=0.025,
        ticks=[0, 1],
    )
    cbar.ax.set_yticklabels(
        ["Not detected", "Detected"],
        fontsize=8,
    )
    cbar.set_label(
        "Detection",
        fontsize=9,
    )

    # Save requested formats.
    if fmt in ("png", "both"):
        png_path = output_dir / f"{database}_detection_heatmap.png"
        fig.savefig(
            png_path,
            bbox_inches="tight",
            dpi=dpi,
        )
        print(f"Success! PNG saved to: {png_path}")

    if fmt in ("pdf", "both"):
        pdf_path = output_dir / f"{database}_detection_heatmap.pdf"
        fig.savefig(
            pdf_path,
            bbox_inches="tight",
        )
        print(f"Success! PDF saved to: {pdf_path}")

    plt.close(fig)


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def generate_heatmap(
    input_dir,
    output_dir,
    database,
    tools="all",
    sort_method="prevalence",
    dpi=300,
    fmt="both",
    cmap="Greys",
    no_cell_lines=False,
):
    """
    Main generation function.

    'tools' is retained for compatibility with the GeneFíor toolkit CLI and
    future integration with tool-specific metadata.
    """
    del tools  # Reserved for future toolkit integration.

    output_path = Path(output_dir)

    heatmap_df = build_combined_matrix(
        input_dir=input_dir,
        database=database,
    )

    detection_df = prepare_detection_matrix(heatmap_df)

    detection_df = sort_detection_matrix(
        detection_df,
        sort_method=sort_method,
    )

    # Save the combined binary detection matrix alongside the figure.
    matrix_path = (
        output_path / f"{database}_combined_detection_matrix.tsv"
    )
    output_path.mkdir(parents=True, exist_ok=True)

    detection_df.to_csv(
        matrix_path,
        sep="\t",
        index_label="Gene",
    )

    print(f"Combined detection matrix saved to: {matrix_path}")

    make_heatmap(
        detection_df=detection_df,
        database=database,
        output_dir=output_path,
        dpi=dpi,
        fmt=fmt,
        cmap=cmap,
        cell_lines=not no_cell_lines,
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "GeneFíor Combine Heatmap: combine per-sample AMR detection "
            "matrices and generate a binary sample × gene heatmap."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-i", "--input",
        required=True, help="Directory containing sample result folders.",)

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Directory in which to save the combined matrix and figures.",
    )

    parser.add_argument(
        "-d",
        "--database",
        required=True,
        help=(
            "Database name. The script searches each sample folder for "
            "<database>_detection_matrix.tsv, e.g. resfinder or card."
        ),
    )

    parser.add_argument(
        "--tools",
        default="all",
        help=(
            "Tools used. Retained for compatibility with GeneFíor; "
            "currently informational/reserved."
        ),
    )

    parser.add_argument(
        "--sort",
        dest="sort_method",
        choices=("prevalence", "alphabetical", "none"),
        default="prevalence",
        help="Gene ordering in the heatmap.",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG resolution.",
    )

    parser.add_argument(
        "--format",
        dest="fmt",
        choices=("png", "pdf", "both"),
        default="both",
        help="Output figure format.",
    )

    parser.add_argument(
        "--cmap",
        default="Greys",
        help=(
            "Matplotlib colour map. For example: Greys, Blues, Purples, "
            "viridis."
        ),
    )

    parser.add_argument(
        "--no-cell-lines",
        action="store_true",
        help="Do not draw fine cell boundaries.",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    try:
        generate_heatmap(
            input_dir=args.input,
            output_dir=args.output,
            database=args.database,
            tools=args.tools,
            sort_method=args.sort_method,
            dpi=args.dpi,
            fmt=args.fmt,
            cmap=args.cmap,
            no_cell_lines=args.no_cell_lines,
        )
    except (FileNotFoundError, NotADirectoryError, RuntimeError) as exc:
        parser.error(str(exc))
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        sys.exit(130)


if __name__ == "__main__":
    main()
