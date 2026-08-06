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


def generate_heatmap(input_dir, output_dir, database, tools):
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
    # We use log1p (log(1+x)) to handle zeros and make subtle differences
    # in low-count samples visible while preventing high-counts from exploding.
    heatmap_df_log = np.log1p(heatmap_df.astype(float))

    # Dynamic sizing
    num_samples = len(heatmap_df_log.columns)
    num_genes = len(heatmap_df_log.index)
    fig_width = min(40, 0.8 * num_samples)
    fig_height = min(30, 0.5 * num_genes)

    plt.figure(figsize=(fig_width, fig_height))

    # 3. COLOR SCALE
    # 'Greys' provides a true White -> Black scale.
    # 'Blues' or 'Purples' are also great alternatives if you want color.
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

    args = parser.parse_args()
    generate_heatmap(args.input, args.output, args.database, args.tools)


if __name__ == "__main__":
    main()
