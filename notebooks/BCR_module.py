import os
import re


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def extract_top_segment(hit_string):
    if pd.isna(hit_string):
        return np.nan

    hit_string = str(hit_string).strip()
    if not hit_string:
        return np.nan

    first_hit = hit_string.split(',')[0].strip()
    first_hit = re.sub(r'\(.*?\)', '', first_hit).strip()
    first_hit = first_hit.split('*')[0].strip()

    return first_hit if first_hit else np.nan

def build_ifn_comparison_df(meta_df, usage_df, top_segments_dict, segment_type):
    frames = []

    usage_sub = usage_df[usage_df['segment_type'] == segment_type].copy()

    for diagnosis in group_order:
        selected_segments = top_segments_dict[(segment_type, diagnosis)]

        meta_sub = meta_df[
            (meta_df['group'] == diagnosis) &
            (meta_df['IFN status'].isin(ifn_order))
        ][['SRR', 'group', 'IFN status']].drop_duplicates().copy()

        if meta_sub.empty or len(selected_segments) == 0:
            continue

        base = pd.MultiIndex.from_product(
            [meta_sub['SRR'].unique(), selected_segments],
            names=['SRR', 'segment']
        ).to_frame(index=False)

        base = base.merge(meta_sub, on='SRR', how='left')

        obs = usage_sub[
            (usage_sub['group'] == diagnosis) &
            (usage_sub['segment'].isin(selected_segments))
        ][['SRR', 'segment', 'segment_freq']].copy()

        merged = base.merge(obs, on=['SRR', 'segment'], how='left')
        merged['segment_freq'] = merged['segment_freq'].fillna(0.0)
        merged['segment_type'] = segment_type

        frames.append(merged)

    return pd.concat(frames, ignore_index=True)


def plot_ifn_boxplots_with_stats(compare_df, stats_df, segment_type, top_segments_dict):
    fig, axes = plt.subplots(1, 4, figsize=(22, 10), sharex=False)

    for j, diagnosis in enumerate(group_order):
        ax = axes[j]

        selected = top_segments_dict[(segment_type, diagnosis)]

        sub = compare_df[
            (compare_df['segment_type'] == segment_type) &
            (compare_df['group'] == diagnosis) &
            (compare_df['segment'].isin(selected))
        ].copy()

        stat_sub = stats_df[
            (stats_df['segment_type'] == segment_type) &
            (stats_df['group'] == diagnosis) &
            (stats_df['segment'].isin(selected))
        ].copy()

        if sub.empty:
            ax.set_visible(False)
            continue

        order = (
            stat_sub.set_index('segment')
            .loc[selected, 'p_adj_fdr']
            .sort_values(na_position='last')
            .index
            .tolist()
        )

        sns.boxplot(
            data=sub,
            y='segment',
            x='segment_freq',
            hue='IFN status',
            order=order,
            hue_order=ifn_order,
            palette=ifn_palette,
            showfliers=False,
            ax=ax
        )

        sns.stripplot(
            data=sub,
            y='segment',
            x='segment_freq',
            hue='IFN status',
            order=order,
            hue_order=ifn_order,
            palette=ifn_palette,
            dodge=True,
            alpha=0.35,
            size=2.5,
            ax=ax
        )

        ax.set_title(f'{segment_type} | {diagnosis}')
        ax.set_ylabel('Segment')
        ax.set_xlabel('Within-sample frequency')

        handles, labels = ax.get_legend_handles_labels()
        uniq = {}
        for h, l in zip(handles, labels):
            if l in ifn_order and l not in uniq:
                uniq[l] = h

        ax.legend(
            [uniq[l] for l in ifn_order if l in uniq],
            [l for l in ifn_order if l in uniq],
            title='IFN status',
            loc='lower right'
        )

        x_max = sub['segment_freq'].max()
        x_pad = x_max * 0.08 if x_max > 0 else 0.02

        for y_pos, seg in enumerate(order):
            row = stat_sub[stat_sub['segment'] == seg]
            if row.empty:
                continue
            label = row['signif'].iloc[0]
            if label == 'NA':
                continue

            local_max = sub.loc[sub['segment'] == seg, 'segment_freq'].max()
            ax.text(local_max + x_pad, y_pos, label, va='center', ha='left', fontsize=12)

        ax.set_xlim(0, x_max + x_pad * 3)

    fig.suptitle(f'{segment_type}: IFN High vs Low for top {top_n} segments', y=1.02)
    plt.tight_layout()
    plt.show()


def build_ifn_compare_df_with_gse(meta_df, usage_df, top_segments_dict, segment_type):
    frames = []

    usage_sub = usage_df[usage_df['segment_type'] == segment_type].copy()

    for diagnosis in group_order:
        selected_segments = top_segments_dict[(segment_type, diagnosis)]

        meta_sub = meta_df[
            (meta_df['group'] == diagnosis) &
            (meta_df['IFN status'].isin(ifn_order))
        ][['SRR', 'GSE', 'group', 'IFN status']].drop_duplicates().copy()

        if meta_sub.empty or len(selected_segments) == 0:
            continue

        base = pd.MultiIndex.from_product(
            [meta_sub['SRR'].unique(), selected_segments],
            names=['SRR', 'segment']
        ).to_frame(index=False)

        base = base.merge(meta_sub, on='SRR', how='left')

        obs = usage_sub[
            (usage_sub['group'] == diagnosis) &
            (usage_sub['segment'].isin(selected_segments))
        ][['SRR', 'GSE', 'segment', 'segment_freq']].copy()

        merged = base.merge(obs, on=['SRR', 'GSE', 'segment'], how='left')
        merged['segment_freq'] = merged['segment_freq'].fillna(0.0)
        merged['segment_type'] = segment_type

        frames.append(merged)

    return pd.concat(frames, ignore_index=True)

def plot_segment_effect_by_gse(gse_stats_df, segment_type, diagnosis, segment):
    sub = gse_stats_df[
        (gse_stats_df['segment_type'] == segment_type) &
        (gse_stats_df['group'] == diagnosis) &
        (gse_stats_df['segment'] == segment)
    ].copy()

    if sub.empty:
        print('Нет данных для выбранного сегмента')
        return

    sub = sub.sort_values('delta_median')

    plt.figure(figsize=(9, max(4, 0.45 * len(sub))))
    sns.barplot(
        data=sub,
        y='GSE',
        x='delta_median',
        color="#336bb9"
    )
    plt.axvline(0, color='black', linestyle='--', linewidth=1)
    plt.xlabel('Median frequency difference: High - Low')
    plt.ylabel('GSE')
    plt.title(f'{segment_type} | {diagnosis} | {segment}')

    for i, (_, row) in enumerate(sub.iterrows()):
        label = f"nL={row['n_low']}, nH={row['n_high']}"
        plt.text(row['delta_median'], i, f'  {label}', va='center')

    plt.tight_layout()
    plt.show()