import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from pathlib import Path
from tqdm import tqdm


def load_sample_data(sample_path, database):
    """
    Loads detection data for a single sample.
    """
    matrix_file = sample_path / f"{database}_detection_matrix.tsv"

    if not matrix_file.exists():
        return None

    try:
        df = pd.read_csv(matrix_file, sep='\t')
        df.columns = [c.lower() for c in df.columns]

        if 'gene' not in df.columns:
            if df.shape[1] > 0:
                df = df.rename(columns={df.columns[0]: 'gene'})
            else:
                return None

        possible_count_cols = ['count', 'reads', 'depth', 'abundance', 'status']
        count_col = None
        for col in possible_count_cols:
            if col in df.columns:
                count_col = col
                break
        if count_col is None and df.shape[1] > 1:
            count_col = df.columns[1]

        return df[['gene', count_col]].copy()
    except Exception as e:
        print(f"Warning: Could not read {matrix_file}: {e}")
        return None


def generate_genome_coverage_plot(input_dir, output_dir, database):
    """Generate whole-genome per-base coverage plot across samples.
    Expects per-sample coverage files in each sample folder matching
    patterns like '<database>_*coverage*.tsv' or '*coverage*.bedgraph'.
    The reader is tolerant to common two-column formats (pos, depth).
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    coverage_series = {}
    print(f"Searching for per-base coverage files for DB '{database}' in {input_dir}...")

    patterns = [f"{database}*coverage*.tsv", f"{database}*coverage*.txt", f"{database}*coverage*.bedgraph", "*coverage*.tsv", "*coverage*.bedgraph"]

    for sample_folder in input_path.iterdir():
        if not sample_folder.is_dir():
            continue
        sample_name = sample_folder.name
        found = False
        for pat in patterns:
            for f in sample_folder.glob(pat):
                try:
                    # Try reading with header first
                    df = pd.read_csv(f, sep='\t', comment='#')
                except Exception:
                    try:
                        df = pd.read_csv(f, sep='\t', header=None)
                    except Exception:
                        continue
                # Normalize to position, depth
                cols = [c.lower() for c in df.columns]
                if 'pos' in cols or 'position' in cols or 'start' in cols:
                    # find pos and depth columns
                    pos_col = None
                    depth_col = None
                    for c in df.columns:
                        lc = str(c).lower()
                        if lc in ('pos','position','start') and pos_col is None:
                            pos_col = c
                        if lc in ('depth','coverage','cov','count','reads') and depth_col is None:
                            depth_col = c
                    if pos_col is None or depth_col is None:
                        # try first two columns
                        pos_col = df.columns[0]
                        depth_col = df.columns[1] if len(df.columns) > 1 else None
                else:
                    # fallback to first two columns
                    pos_col = df.columns[0]
                    depth_col = df.columns[1] if len(df.columns) > 1 else None

                if depth_col is None:
                    continue

                try:
                    pos = df[pos_col].astype(int).to_numpy()
                    depth = df[depth_col].astype(float).to_numpy()
                except Exception:
                    # coerce using values
                    arr = df.to_numpy()
                    if arr.shape[1] < 2:
                        continue
                    pos = arr[:,0].astype(int)
                    depth = arr[:,1].astype(float)

                if len(pos) == 0:
                    continue

                # Downsample long genomes for plotting clarity
                max_points = 200000
                if len(pos) > max_points:
                    step = max(1, len(pos) // max_points)
                    pos = pos[::step]
                    depth = depth[::step]

                coverage_series[sample_name] = (pos, depth)
                found = True
                break
            if found:
                break

    if not coverage_series:
        print(f"No per-base coverage files found for database '{database}' in samples under {input_dir}.")
        return

    # Combined plot
    plt.figure(figsize=(16, 6))
    for sample_name, (pos, depth) in coverage_series.items():
        plt.plot(pos, depth, label=sample_name, alpha=0.7, linewidth=0.8)

    plt.xlabel('Genome position')
    plt.ylabel('Coverage (depth)')
    plt.title(f'Whole-genome coverage: {database}')
    plt.legend(fontsize='small', ncol=2)
    plt.tight_layout()
    out_file = Path(output_dir) / f"{database}_whole_genome_coverage.png"
    plt.savefig(out_file, dpi=200)
    print(f"Saved whole-genome coverage plot to {out_file}")


def generate_heatmap(input_dir, output_dir, database, tools, db_whole_genome=False):
    # If the database is declared as whole-genome, generate coverage plots instead
    if db_whole_genome:
        generate_genome_coverage_plot(input_dir, output_dir, database)
        return

    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    all_data = []
    print(f"Searching for samples in {input_dir}...")
    for sample_folder in input_path.iterdir():
        if sample_folder.is_dir():
            df = load_sample_data(sample_folder, database)
            if df is not None and not df.empty:
                df['sample_id'] = sample_folder.name
                all_data.append(df)

    if not all_data:
        print(f"No data found for database '{database}' in {input_dir}.")
        return

    combined_df = pd.concat(all_data, ignore_index=True)
    value_col_name = combined_df.columns[1]

    # Pivot the data
    heatmap_df = combined_df.pivot(index='gene', columns='sample_id', values=value_col_name)
    heatmap_df = heatmap_df.fillna(0)

    # 1. Filter out genes with zero reads across all samples
    heatmap_df = heatmap_df[heatmap_df.sum(axis=1) > 0]

    # 2. LOG TRANSFORMATION
    heatmap_df_log = np.log1p(heatmap_df.astype(float))

    # Dynamic sizing
    num_samples = len(heatmap_df_log.columns)
    num_genes = len(heatmap_df_log.index)
    fig_width = min(40, 0.8 * num_samples)
    fig_height = min(30, 0.5 * num_genes)

    plt.figure(figsize=(fig_width, fig_height))

    sns.heatmap(heatmap_df_log,
                annot=False,
                cmap="Greys",
                linewidths=0.5,
                cbar_kws={'label': 'Log1p(Read Count)'})

    plt.title(f"AMR Gene Detection (Log Scaled)\nDB: {database.upper()}", pad=20)
    plt.xlabel("Samples")
    plt.ylabel("Genes")

    plt.xticks(rotation=90, ha='center', fontsize=8)
    plt.yticks(rotation=0, fontsize=8)

    plot_file = output_path / f"{database}_heatmap_log.png"
    plt.savefig(plot_file, bbox_inches='tight', dpi=300)
    print(f"Success! Log-scaled heatmap saved to: {plot_file}")


def main():
    parser = argparse.ArgumentParser(description='GeneFíor Heatmap: Generate log-scaled read count visualisations.')
    parser.add_argument('-i', '--input', required=True, help='Directory containing sample result folders')
    parser.add_argument('-o', '--output', required=True, help='Directory to save plots')
    parser.add_argument('-d', '--database', required=True, help='Database to visualize')
    parser.add_argument('--tools', default='all', help='Tools used')
    parser.add_argument('--db-whole-genome', action='store_true', dest='db_whole_genome', default=False,
                        help='Indicate that the database is a whole-genome reference and generate whole-genome coverage plots')

    args = parser.parse_args()
    generate_heatmap(args.input, args.output, args.database, args.tools, db_whole_genome=args.db_whole_genome)


if __name__ == "__main__":
    main()
