"""Utilities for loading, processing, and visualizing BCR clonotype data"""

from __future__ import annotations

import glob
import os
import re
import tarfile
from collections import defaultdict
from typing import Any

import gdown
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def download_clonotype_archive(
    file_id: str,
    archive_path: str = '../data/clonotypes.tar.gz',
) -> str:
    """
    Download a clonotype archive from Google Drive if it is not available locally.

    Parameters
    ----------
    file_id : str
        Google Drive file ID.
    archive_path : str, optional
        Local path where the downloaded archive will be stored.

    Returns
    -------
    str
        Path to the downloaded or existing archive.
    """
    os.makedirs(os.path.dirname(archive_path), exist_ok=True)

    if not os.path.exists(archive_path):
        print(f'Downloading archive to: {archive_path}')
        gdown.download(id=file_id, output=archive_path, quiet=False)
    else:
        print(f'Archive already exists: {archive_path}')

    return archive_path


def extract_clonotype_archive(
    archive_path: str,
    extract_dir: str = '../data/clonotypes',
) -> str:
    """
    Extract a clonotype tar.gz archive into the target directory.

    Extraction is performed only if the target directory is empty.

    Parameters
    ----------
    archive_path : str
        Path to the clonotype archive.
    extract_dir : str, optional
        Directory where archive contents will be extracted.

    Returns
    -------
    str
        Path to the extraction directory.
    """
    os.makedirs(extract_dir, exist_ok=True)

    if not os.listdir(extract_dir):
        print(f'Extracting archive to: {extract_dir}')
        with tarfile.open(archive_path, 'r:gz') as tar:
            tar.extractall(path=extract_dir)
    else:
        print(f'Extraction directory is not empty: {extract_dir}')

    return extract_dir


def collect_clonotype_files(
    samples_df: pd.DataFrame,
    extract_dir: str,
) -> tuple[list[str], list[str], dict[str, list[str]]]:
    """
    Match clonotype tables in the extracted directory to SRR IDs from samples_df.

    Supported filename formats
    ---------------------------
    - SRR123456.clones_IGH.tsv
    - SRR123456_trimmed.clones_IGH.tsv

    If multiple matching files are found for the same SRR, the function
    preferentially selects the '_trimmed.clones_IGH.tsv' file when exactly
    one such file is present.

    Parameters
    ----------
    samples_df : pandas.DataFrame
        Sample metadata table containing the 'SRR' column.
    extract_dir : str
        Root directory containing extracted clonotype files.

    Returns
    -------
    tuple[list[str], list[str], dict[str, list[str]]]
        A tuple containing:
        - file_list: matched clonotype file paths
        - missing_samples: SRR IDs for which no file was found
        - duplicate_samples: SRR IDs with multiple candidate files
    """
    all_files = glob.glob(
        os.path.join(extract_dir, '**', '*.clones_IGH.tsv'),
        recursive=True,
    )

    print(f'Total clonotype files found: {len(all_files)}')

    srr_to_files: dict[str, list[str]] = defaultdict(list)

    for file_path in all_files:
        filename = os.path.basename(file_path)
        match = re.match(r'^(SRR\d+)(?:_trimmed)?\.clones_IGH\.tsv$', filename)

        if match:
            srr = match.group(1)
            srr_to_files[srr].append(file_path)

    file_list: list[str] = []
    missing_samples: list[str] = []
    duplicate_samples: dict[str, list[str]] = {}

    for srr in samples_df['SRR'].astype(str).str.strip():
        matches = srr_to_files.get(srr, [])

        if len(matches) == 1:
            file_list.append(matches[0])

        elif len(matches) > 1:
            trimmed_matches = [
                path
                for path in matches
                if os.path.basename(path).endswith('_trimmed.clones_IGH.tsv')
            ]

            if len(trimmed_matches) == 1:
                file_list.append(trimmed_matches[0])
            else:
                duplicate_samples[srr] = matches
                file_list.append(matches[0])
        else:
            missing_samples.append(srr)

    print(f'Found files: {len(file_list)}')
    print(f'Missing samples: {len(missing_samples)}')
    print(f'Samples with multiple candidate files: {len(duplicate_samples)}')

    return file_list, missing_samples, duplicate_samples


def extract_top_segment(hit_string: Any) -> str | float:
    """
    Extract the top V/D/J segment name from a MiXCR hit string.

    The function keeps only the first hit, removes alignment scores in
    parentheses, and strips allele annotation after '*'.

    Example
    -------
    'IGHV3-23*01(320.7),IGHV3-30*01(303.0)' -> 'IGHV3-23'

    Parameters
    ----------
    hit_string : Any
        Raw MiXCR hit string.

    Returns
    -------
    str | float
        Segment name or numpy.nan if parsing fails.
    """
    if pd.isna(hit_string):
        return float('nan')

    hit_string = str(hit_string).strip()
    if not hit_string:
        return float('nan')

    first_hit = hit_string.split(',')[0].strip()
    first_hit = re.sub(r'\(.*?\)', '', first_hit).strip()
    first_hit = first_hit.split('*')[0].strip()

    return first_hit if first_hit else float('nan')


def build_ifn_comparison_df(
    meta_df: pd.DataFrame,
    usage_df: pd.DataFrame,
    top_segments_dict: dict[tuple[str, str], list[str]],
    segment_type: str,
    group_order: list[str],
    ifn_order: list[str],
) -> pd.DataFrame:
    """
    Build a zero-filled sample-segment table for IFN group comparisons.

    For each diagnostic group and selected top segment, the output includes
    all sample-segment combinations. Missing segments are assigned frequency 0.

    Parameters
    ----------
    meta_df : pandas.DataFrame
        Sample metadata containing 'SRR', 'group', and 'IFN status'.
    usage_df : pandas.DataFrame
        Per-sample segment usage table.
    top_segments_dict : dict[tuple[str, str], list[str]]
        Mapping from (segment_type, group) to selected top segments.
    segment_type : str
        Segment class to analyze ('IGHV', 'IGHD', or 'IGHJ').
    group_order : list[str]
        Ordered list of diagnostic groups to process.
    ifn_order : list[str]
        Ordered list of IFN states to retain.

    Returns
    -------
    pandas.DataFrame
        Zero-filled comparison table for the requested segment type.
    """
    frames: list[pd.DataFrame] = []
    usage_sub = usage_df[usage_df['segment_type'] == segment_type].copy()

    for diagnosis in group_order:
        selected_segments = top_segments_dict.get((segment_type, diagnosis), [])
        if not selected_segments:
            continue

        meta_sub = meta_df[
            (meta_df['group'] == diagnosis)
            & (meta_df['IFN status'].isin(ifn_order))
        ][['SRR', 'group', 'IFN status']].drop_duplicates()

        if meta_sub.empty:
            continue

        base = pd.MultiIndex.from_product(
            [meta_sub['SRR'].unique(), selected_segments],
            names=['SRR', 'segment'],
        ).to_frame(index=False)

        base = base.merge(meta_sub, on='SRR', how='left')

        obs = usage_sub[
            (usage_sub['group'] == diagnosis)
            & (usage_sub['segment'].isin(selected_segments))
        ][['SRR', 'segment', 'segment_freq']].copy()

        merged = base.merge(obs, on=['SRR', 'segment'], how='left')
        merged['segment_freq'] = merged['segment_freq'].fillna(0.0)
        merged['segment_type'] = segment_type

        frames.append(merged)

    if not frames:
        return pd.DataFrame(
            columns=[
                'SRR',
                'segment',
                'group',
                'IFN status',
                'segment_freq',
                'segment_type',
            ]
        )

    return pd.concat(frames, ignore_index=True)


def build_ifn_compare_df_with_gse(
    meta_df: pd.DataFrame,
    usage_df: pd.DataFrame,
    top_segments_dict: dict[tuple[str, str], list[str]],
    segment_type: str,
    group_order: list[str],
    ifn_order: list[str],
) -> pd.DataFrame:
    """
    Build a zero-filled sample-segment table with GSE labels preserved.

    This table is designed for assessing whether IFN-associated segment
    differences are consistent across independent datasets.

    Parameters
    ----------
    meta_df : pandas.DataFrame
        Sample metadata containing 'SRR', 'GSE', 'group', and 'IFN status'.
    usage_df : pandas.DataFrame
        Per-sample segment usage table.
    top_segments_dict : dict[tuple[str, str], list[str]]
        Mapping from (segment_type, group) to selected top segments.
    segment_type : str
        Segment class to analyze ('IGHV', 'IGHD', or 'IGHJ').
    group_order : list[str]
        Ordered list of diagnostic groups to process.
    ifn_order : list[str]
        Ordered list of IFN states to retain.

    Returns
    -------
    pandas.DataFrame
        Zero-filled comparison table with GSE labels.
    """
    frames: list[pd.DataFrame] = []
    usage_sub = usage_df[usage_df['segment_type'] == segment_type].copy()

    for diagnosis in group_order:
        selected_segments = top_segments_dict.get((segment_type, diagnosis), [])
        if not selected_segments:
            continue

        meta_sub = meta_df[
            (meta_df['group'] == diagnosis)
            & (meta_df['IFN status'].isin(ifn_order))
        ][['SRR', 'GSE', 'group', 'IFN status']].drop_duplicates()

        if meta_sub.empty:
            continue

        base = pd.MultiIndex.from_product(
            [meta_sub['SRR'].unique(), selected_segments],
            names=['SRR', 'segment'],
        ).to_frame(index=False)

        base = base.merge(meta_sub, on='SRR', how='left')

        obs = usage_sub[
            (usage_sub['group'] == diagnosis)
            & (usage_sub['segment'].isin(selected_segments))
        ][['SRR', 'GSE', 'segment', 'segment_freq']].copy()

        merged = base.merge(obs, on=['SRR', 'GSE', 'segment'], how='left')
        merged['segment_freq'] = merged['segment_freq'].fillna(0.0)
        merged['segment_type'] = segment_type

        frames.append(merged)

    if not frames:
        return pd.DataFrame(
            columns=[
                'SRR',
                'GSE',
                'segment',
                'group',
                'IFN status',
                'segment_freq',
                'segment_type',
            ]
        )

    return pd.concat(frames, ignore_index=True)


def plot_ifn_boxplots_with_stats(
    compare_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    segment_type: str,
    top_segments_dict: dict[tuple[str, str], list[str]],
    group_order: list[str],
    ifn_order: list[str],
    ifn_palette: dict[str, str],
    top_n: int,
) -> None:
    """
    Plot IFN-Low vs IFN-High boxplots for selected segments across groups.

    Parameters
    ----------
    compare_df : pandas.DataFrame
        Zero-filled sample-segment table for plotting.
    stats_df : pandas.DataFrame
        Statistical results table containing FDR-adjusted p-values and
        significance labels.
    segment_type : str
        Segment class to plot ('IGHV', 'IGHD', or 'IGHJ').
    top_segments_dict : dict[tuple[str, str], list[str]]
        Mapping from (segment_type, group) to selected top segments.
    group_order : list[str]
        Ordered list of diagnostic groups to display.
    ifn_order : list[str]
        Ordered list of IFN states.
    ifn_palette : dict[str, str]
        Mapping from IFN states to colors.
    top_n : int
        Number of top segments shown in the title.

    Returns
    -------
    None
    """
    n_groups = len(group_order)
    fig, axes = plt.subplots(1, n_groups, figsize=(7 * n_groups, 8), sharex=False)

    if n_groups == 1:
        axes = [axes]

    for j, diagnosis in enumerate(group_order):
        ax = axes[j]
        selected = top_segments_dict.get((segment_type, diagnosis), [])

        sub = compare_df[
            (compare_df['segment_type'] == segment_type)
            & (compare_df['group'] == diagnosis)
            & (compare_df['segment'].isin(selected))
        ].copy()

        stat_sub = stats_df[
            (stats_df['segment_type'] == segment_type)
            & (stats_df['group'] == diagnosis)
            & (stats_df['segment'].isin(selected))
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
            ax=ax,
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
            ax=ax,
        )

        ax.set_title(f'{segment_type} | {diagnosis}')
        ax.set_ylabel('Segment')
        ax.set_xlabel('Within-sample frequency')

        handles, labels = ax.get_legend_handles_labels()
        unique_handles: dict[str, Any] = {}
        for handle, label in zip(handles, labels):
            if label in ifn_order and label not in unique_handles:
                unique_handles[label] = handle

        ax.legend(
            [unique_handles[label] for label in ifn_order if label in unique_handles],
            [label for label in ifn_order if label in unique_handles],
            title='IFN status',
            loc='lower right',
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
            ax.text(
                local_max + x_pad,
                y_pos,
                label,
                va='center',
                ha='left',
                fontsize=12,
            )

        ax.set_xlim(0, x_max + x_pad * 3)

    fig.suptitle(
        f'{segment_type}: IFN High vs Low for top {top_n} segments',
        y=1.02,
    )
    plt.tight_layout()
    plt.show()


def plot_segment_effect_by_gse(
    gse_stats_df: pd.DataFrame,
    segment_type: str,
    diagnosis: str,
    segment: str,
) -> None:
    """
    Plot dataset-specific effect sizes for one segment.

    Parameters
    ----------
    gse_stats_df : pandas.DataFrame
        Table with per-dataset comparison statistics.
    segment_type : str
        Segment class ('IGHV', 'IGHD', or 'IGHJ').
    diagnosis : str
        Diagnostic group to display.
    segment : str
        Segment name to display.

    Returns
    -------
    None
    """
    sub = gse_stats_df[
        (gse_stats_df['segment_type'] == segment_type)
        & (gse_stats_df['group'] == diagnosis)
        & (gse_stats_df['segment'] == segment)
    ].copy()

    if sub.empty:
        print('No data available for the selected segment.')
        return

    sub = sub.sort_values('delta_median')

    plt.figure(figsize=(9, max(4, 0.45 * len(sub))))
    sns.barplot(
        data=sub,
        y='GSE',
        x='delta_median',
        color='#336BB9',
    )
    plt.axvline(0, color='black', linestyle='--', linewidth=1)
    plt.xlabel('Median frequency difference: High - Low')
    plt.ylabel('GSE')
    plt.title(f'{segment_type} | {diagnosis} | {segment}')

    for i, (_, row) in enumerate(sub.iterrows()):
        label = f'nL={row["n_low"]}, nH={row["n_high"]}'
        plt.text(row['delta_median'], i, f'  {label}', va='center')

    plt.tight_layout()
    plt.show()