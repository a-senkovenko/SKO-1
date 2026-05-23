import math

import pandas as pd
import matplotlib.pyplot as plt


def prepare_pie_data(
    df_counts: pd.DataFrame,
    gene: str,
    min_top_n: int = 5,
    min_pct: float = 10,
    target_pct: float = 90,
) -> tuple[list[str], list[float]]:

    # Aggregate data by alleles to exclude duplicates.
    gene_data = (
        df_counts.loc[df_counts['gene'] == gene]
        .groupby('allele', as_index=False)['allele_copies']
        .sum()
        .sort_values('allele_copies', ascending=False, ignore_index=True)
    )

    total = gene_data['allele_copies'].sum()

    gene_data['pct'] = gene_data['allele_copies'] / total * 100

    # guaranteed top-N alleles
    top_n = min(min_top_n, len(gene_data))
    selected = gene_data.iloc[:top_n].copy()
    current_pct = selected['pct'].sum()

    # additional alleles
    if current_pct < target_pct:
        remaining = gene_data.iloc[top_n:]
        additional = remaining[remaining['pct'] >= min_pct]
        cumulative_pct = current_pct

        extra_rows = []

        for row in additional.itertuples():
            if cumulative_pct >= target_pct:
                break

            extra_rows.append(row.Index)
            cumulative_pct += row.pct

        if extra_rows:
            selected = pd.concat(
                [selected, gene_data.loc[extra_rows]], ignore_index=True
            )

    # other
    other = gene_data.loc[~gene_data['allele'].isin(selected['allele'])]

    if not other.empty:
        other_row = pd.DataFrame(
            {'allele': ['Other'], 'allele_copies': [other['allele_copies'].sum()]}
        )

        selected = pd.concat([selected, other_row], ignore_index=True)

    return (
        selected['allele'].tolist(),
        selected['allele_copies'].tolist(),
    )


def make_pie(
    df_counts: pd.DataFrame,
    genes: str | list[str],
    cohort_name: str = 'cohort',
) -> None:

    def autopct_func(pct: float) -> str:
        return (f'{pct:.1f}%') if pct > 5 else ''

    if isinstance(genes, str):
        plt.figure(figsize=(8, 6))
        labels, counts = prepare_pie_data(df_counts, genes)
        plt.pie(
            counts, labels=None, autopct=autopct_func, pctdistance=0.85, startangle=270
        )
        plt.legend(
            labels, loc='center left', bbox_to_anchor=(1, 0, 0.5, 1), fontsize=10
        )
        plt.title(f'Pie Chart for {genes} gene', fontsize=14)
        plt.axis('equal')

    elif isinstance(genes, list):
        ncols = 3
        nrows = math.ceil(len(genes) / ncols)
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(ncols * 8, max(5, nrows * 6)), facecolor='white'
        )
        axes = axes.flatten()

        for i, g in enumerate(genes):
            ax = axes[i]
            labels, counts = prepare_pie_data(df_counts, g)
            if not counts:
                ax.set_visible(False)
                continue

            ax.pie(
                counts,
                labels=None,
                autopct=autopct_func,
                pctdistance=0.75,
                startangle=270,
                textprops={'fontsize': 16},
                wedgeprops={'width': 0.6},
            )
            ax.text(0, 0, g, ha='center', va='center', fontsize=24, fontweight='bold')
            ax.legend(labels, loc='upper left', bbox_to_anchor=(1, 1), fontsize=14)

        # Hide unused subplots
        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        fig.suptitle(
            f"Diversity of HLA alleles by genes for {cohort_name}\n",
            fontsize=32,
            fontweight='bold',
        )

        plt.tight_layout(pad=2, h_pad=1.2, w_pad=2)
        plt.show()
