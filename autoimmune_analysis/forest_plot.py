import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import matplotlib.transforms as transforms
from matplotlib.ticker import ScalarFormatter


# ============================================================
# SETTINGS
# ============================================================

SIGNIFICANCE_ALPHA = 0.05

DIAMOND_HEIGHT = 0.007

GENE_GAP = 0.1
ALLELE_GAP = 0.005
ROW_GAP = 0.05

POPULATION_COLORS = {
    'H': '#f1640cf5',
    'C': '#474df1',
    'AA': '#2ca02c',
}

DEFAULT_CMAP = plt.cm.tab10


# ============================================================
# PREPARE DATA
# ============================================================


def prepare_data(df):

    df = df.copy()
    df = df[(df['OR'] > 0) & (df['CI_lower'] > 0) & (df['CI_upper'] > 0)]

    df['significant'] = ~((df['CI_lower'] <= 1) & (df['CI_upper'] >= 1))

    df['log_or'] = np.log(df['OR'])

    df['SE'] = (np.log(df['CI_upper']) - np.log(df['CI_lower'])) / (2 * 1.96)

    return df


# ============================================================
# SUMMARY EFFECT
# ============================================================


def compute_summary_effect(subdf):

    log_or = subdf['log_or'].values
    se = subdf['SE'].values

    weights = 1 / (se**2)

    pooled_log_or = np.sum(weights * log_or) / np.sum(weights)

    pooled_se = np.sqrt(1 / np.sum(weights))

    lower = pooled_log_or - 1.96 * pooled_se
    upper = pooled_log_or + 1.96 * pooled_se

    return {
        'OR': np.exp(pooled_log_or),
        'CI_lower': np.exp(lower),
        'CI_upper': np.exp(upper),
        'significant': not (np.exp(lower) <= 1 <= np.exp(upper)),
    }


# ============================================================
# BUILD PLOTTING TABLE
# ============================================================


def build_plot_rows(df):

    rows = []
    grouped = df.groupby(['gene', 'allele'], sort=False)

    for (gene, allele), subdf in grouped:
        for _, r in subdf.iterrows():
            rows.append(
                {
                    'type': 'population',
                    'gene': gene,
                    'allele': allele,
                    'population': r['population'],
                    'OR': r['OR'],
                    'CI_lower': r['CI_lower'],
                    'CI_upper': r['CI_upper'],
                    'significant': r['significant'],
                }
            )

        summary = compute_summary_effect(subdf)

        rows.append(
            {
                'type': 'summary',
                'gene': gene,
                'allele': allele,
                'population': 'Summary',
                **summary,
            }
        )

    return pd.DataFrame(rows)


# ============================================================
# DIAMOND
# ============================================================


def draw_diamond(ax, center_x, lower, upper, y, height=0.08, color='black'):

    diamond = Polygon(
        [
            (lower, y),
            (center_x, y + height),
            (upper, y),
            (center_x, y - height),
        ],
        closed=True,
        facecolor=color,
        edgecolor=color,
        lw=0.8,
        alpha=0.9,
        zorder=4,
    )

    ax.add_patch(diamond)


# ============================================================
# CI CLIPPING
# ============================================================


def clip_ci(lower, upper, xmin, xmax):

    clipped_lower = max(lower, xmin)
    clipped_upper = min(upper, xmax)

    left_trunc = lower < xmin
    right_trunc = upper > xmax

    return (
        clipped_lower,
        clipped_upper,
        left_trunc,
        right_trunc,
    )


# ============================================================
# MAIN PLOT
# ============================================================


def plot_forest(
    df,
    diag='Diag',
    figsize=(9, 7),
    xmin=0.05,
    xmax=50,
):
    """
    Create and display a forest plot of odds ratios with 95% confidence intervals.

    Args:
        df (pd.DataFrame): Input data containing columns "gene", "allele", "population",
            "OR", "CI_lower", and "CI_upper".
        diag (str, optional): Diagnosis label used in the plot title. Defaults to 'Diag'.
        figsize (tuple, optional): Figure size passed to matplotlib. Defaults to (9, 7).
        xmin (float, optional): Minimum x-axis limit for the odds ratio display. Defaults to 0.05.
        xmax (float, optional): Maximum x-axis limit for the odds ratio display. Defaults to 50.

    Returns:
        None: The plot is shown with plt.show().
    """

    df = prepare_data(df)
    plot_df = build_plot_rows(df)

    # COLORS
    populations = plot_df.loc[plot_df['type'] == 'population', 'population'].unique()

    color_map = {}
    for i, pop in enumerate(populations):
        if pop in POPULATION_COLORS:
            color_map[pop] = POPULATION_COLORS[pop]
        else:
            color_map[pop] = DEFAULT_CMAP(i)

    # Y POSITIONS
    y_positions = []
    y = 1

    prev_gene = None
    prev_allele = None

    for _, row in plot_df.iterrows():
        gene = row['gene']
        allele = row['allele']

        if gene != prev_gene:
            y += GENE_GAP

        elif allele != prev_allele:
            y += ALLELE_GAP

        y += ROW_GAP

        y_positions.append(y)

        prev_gene = gene
        prev_allele = allele

    plot_df['y'] = y_positions

    # FIGURE
    fig, ax = plt.subplots(figsize=figsize)

    # REFERENCE LINE
    ax.axvline(
        1,
        color='black',
        linestyle='--',
        lw=1,
        alpha=0.7,
        zorder=1,
    )

    # DRAW EFFECTS
    for _, row in plot_df.iterrows():
        y = row['y']
        OR = row['OR']
        lower = row['CI_lower']
        upper = row['CI_upper']

        lower_clip, upper_clip, left_trunc, right_trunc = clip_ci(
            lower, upper, xmin, xmax
        )

        if row['type'] == 'summary':
            # diamond
            draw_diamond(
                ax=ax,
                center_x=OR,
                lower=lower_clip,
                upper=upper_clip,
                y=y,
                height=DIAMOND_HEIGHT,
                color='black',
            )

        else:
            color = color_map[row['population']]

            # CI line
            ax.plot(
                [lower_clip, upper_clip],
                [y, y],
                color=color,
                lw=1.8,
                solid_capstyle='round',
                zorder=2,
            )

            # arrows
            if right_trunc:
                ax.plot(
                    xmax,
                    y,
                    marker='>',
                    color=color,
                    markersize=5,
                    clip_on=False,
                )

            if left_trunc:
                ax.plot(
                    xmin,
                    y,
                    marker='<',
                    color=color,
                    markersize=5,
                    clip_on=False,
                )

            # filled vs hollow
            facecolor = color if row['significant'] else 'white'

            ax.scatter(
                OR,
                y,
                s=38,
                marker='o',
                edgecolor=color,
                facecolor=facecolor,
                linewidth=1.4,
                zorder=3,
            )

    # AXES
    ax.set_xscale('log')
    ax.set_xlim(xmin, xmax)

    xticks = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, xmax]
    ax.set_xticks(xticks)

    ax.get_xaxis().set_major_formatter(ScalarFormatter())

    ax.tick_params(
        axis='x',
        labelsize=10,
    )

    ax.set_yticks([])
    ax.invert_yaxis()

    ax.set_xlabel(
        'Odds Ratio (95% CI)',
        fontsize=12,
        fontweight='bold',
        labelpad=10,
    )

    ax.set_title(
        f"HLA Allele Associations with IFN Status for {diag} \n",
        fontsize=18,
        fontweight='bold',
        pad=16,
    )

    # LEFT LABELS
    trans = transforms.blended_transform_factory(
        ax.transAxes,
        ax.transData,
    )

    seen_gene = set()
    seen_allele = set()

    allele_y_map = (
        plot_df.groupby(['gene', 'allele'])['y']
        .mean()
        .to_dict()
    )

    for _, row in plot_df.iterrows():
        y = row['y']

        gene = row['gene']
        allele = row['allele']

        # gene label
        if gene not in seen_gene:
            ax.text(
                -0.08,
                y - 0.05,
                gene,
                transform=trans,
                ha='right',
                va='center',
                fontsize=14,
                fontweight='bold',
            )

            seen_gene.add(gene)

        # allele label
        if (gene, allele) not in seen_allele:
            allele_y = allele_y_map[(gene, allele)]

            ax.text(
                -0.08,
                allele_y,
                allele,
                transform=trans,
                ha='right',
                va='center',
                fontsize=12,
            )

            seen_allele.add((gene, allele))

    # RIGHT LABELS
    ax.text(
        1.05,
        0.99,
        'OR (95% CI)',
        transform=ax.transAxes,
        ha='left',
        va='bottom',
        fontsize=12,
        fontweight='bold',
    )

    for _, row in plot_df.iterrows():
        y = row['y']

        txt = f'{row["OR"]:.2f} ({row["CI_lower"]:.2f}–{row["CI_upper"]:.2f})'

        ax.text(
            1.05,
            y,
            txt,
            transform=trans,
            ha='left',
            va='center',
            fontsize=9.5,
        )

        # IFN LABELS
        ax.text(
            0.10,
            -0.085,
            'IFN Low (OR < 1)',
            transform=ax.transAxes,
            color='#1f77b4',
            fontsize=10,
            ha='center',
        )

    ax.text(
        0.90,
        -0.085,
        'IFN High (OR > 1)',
        transform=ax.transAxes,
        color='#d62728',
        fontsize=10,
        ha='center',
    )

    # LEGEND
    legend_elements = []

    for pop, color in color_map.items():
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker='o',
                color=color,
                label=pop,
                markerfacecolor=color,
                markersize=6,
                lw=1.8,
            )
        )

    legend_elements.append(
        Line2D(
            [0],
            [0],
            marker='o',
            color='black',
            label='Non-significant',
            markerfacecolor='white',
            markersize=6,
            lw=0,
        )
    )

    ax.legend(
        handles=legend_elements,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.18),
        frameon=False,
        ncol=len(legend_elements),
        fontsize=9,
    )

    # GRID
    ax.grid(
        axis='x',
        linestyle=':',
        alpha=0.35,
    )

    # REMOVE SPINES
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # LAYOUT
    plt.subplots_adjust(
        left=0.27,
        right=0.80,
        top=0.90,
        bottom=0.16,
    )

    # NOTE
    if xmax < df['CI_upper'].max():
        fig.text(
            0.07,
            0.01,
            f"Confidence intervals truncated to [{xmin}, {xmax}] for visualization.",
            fontsize=8.5,
            color='gray',
        )

    plt.show()
