import matplotlib.pyplot as plt
import math


def prepare_pie_data(df_counts, gene, min_top_n=5, min_pct=10, target_pct=90):

    # Aggregate data by alleles to exclude duplicates.
    gene_data = (
        df_counts[df_counts['gene'] == gene]
        .groupby('allele')['allele_copies']
        .sum()
        .reset_index()
    )
    gene_data = gene_data.sort_values('allele_copies', ascending=False).reset_index(
        drop=True
    )

    labels = gene_data['allele'].tolist()
    values = gene_data['allele_copies'].tolist()
    total = sum(values)

    percentages = [(v / total) * 100 for v in values]

    min_top_n = min(min_top_n, len(values))
    selected_indices = [i for i in range(min_top_n)]
    selected_labels = labels[:min_top_n]
    selected_values = values[:min_top_n]
    selected_pcts = percentages[:min_top_n]

    current_sum_pct = sum(selected_pcts)

    if current_sum_pct < target_pct:
        i = min_top_n
        while (
            i < len(values)
            and percentages[i] >= min_pct
            and current_sum_pct < target_pct
        ):
            selected_indices.append(i)
            selected_values.append(values[i])
            selected_labels.append(labels[i])
            selected_pcts.append(percentages[i])
            current_sum_pct = sum(selected_pcts)
            i += 1

    other_indices = [j for j in range(len(values)) if j not in selected_indices]

    if other_indices:
        other_value = sum(values[j] for j in other_indices)
        other_pct = (other_value / total) * 100
        selected_values.append(other_value)
        selected_labels.append('Other')
        selected_pcts.append(other_pct)

    return selected_labels, selected_values


def make_pie(df_counts, genes, cohort_name='cohort'):

    def autopct_func(pct):
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

        fig.subplots_adjust(wspace=0.8, top=0.92)
        plt.tight_layout(pad=2, h_pad=1.2, w_pad=2)
        plt.show()
