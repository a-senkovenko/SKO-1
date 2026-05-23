import pandas as pd
import numpy as np

from skbio.stats.composition import clr, multi_replace

import os
import sys
import re
import glob

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

def get_files_dict(
    results_dir: str,
    genes: list = ["TRA", "TRB", "TRD", "TRG"],
    metadata: pd.DataFrame = None
):
    """
    A function to generate a dictionary of filename lists.

    - results_dir: a path to the directory with the results files
    - genes: a list of gene names to load the results for
    - metadata: a pandas DataFrame with 'SRR' column to mark which samples to load

    Reads filenames recursively from a results_dir. Outputs a dictionary with lists of filenames for each gene
    
    Returns:
        dict: {gene: [filepath, ...], ...}
        Example: {"TRB": ["path/SRR123.clones_TRB.tsv", ...], "TRA": [...]}
    """
    files = {}
    valid_srr = set(metadata["SRR"].astype(str))

    for gene in genes:
        filepath = f"SRR[0-9]*_mixcr/SRR[0-9]*.clones_{gene}.tsv"
        pattern = os.path.join(results_dir, filepath)

        all_files = glob.glob(pattern, recursive=True)

        filtered = []

        for f in all_files:
            if "_trimmed_mixcr/" in f:
                continue

            filename = os.path.basename(f) 
            srr = filename.split(".clones_")[0]

            if srr in valid_srr:
                filtered.append(f)

        files[gene] = filtered

    return files

SEGMENTS = {
    "V": "allVHitsWithScore",
    "D": "allDHitsWithScore",
    "J": "allJHitsWithScore"
}

def get_gene_count_matrix(
        gene: str,
        segment: str,
        results_dir: str,
        metadata: pd.DataFrame = None
) -> pd.DataFrame:
    """
    A function to generate a gene count matrix for V, D or J segment genes across samples.

    - gene: a string specifying for which gene to generate the count matrix. Must be one of the 'TRB', 'TRA', 'TRD', 'TRG'
    - segment: a string specifying for which segment to get counts for. Must be one of the 'V', 'D' or 'J'
    - results_dir: a path to the directory with the results files to load data from
    - metadata: a pandas DataFrame with 'SRR' column to mark which samples to load

    Returns a pd.DataFrame with genes as rows and samples as columns. 
    Values are raw counts from readCount columns of MiXCR clonotype tables.
    """
    if gene not in ("TRA", "TRB", "TRD", "TRG"):
        raise ValueError(f"gene must be one of TRA, TRB, TRD, TRG, got '{gene}'")
    if segment not in SEGMENTS:
        raise ValueError(f"segment must be one of {list(SEGMENTS.keys())}, got '{segment}'")
    
    files = get_files_dict(results_dir, [gene], metadata)
    all_counts = []

    for file in files[gene]:
        report = pd.read_csv(file, sep="\t")

        report = report[
            ~report.aaSeqCDR3.str.contains(r"[_*]", na=False)
            ]
        
        srr = os.path.basename(file).split(".")[0]
        
        report["gene"] = report[SEGMENTS[segment]].str.split("*").str[0]

        gene_counts = (
            report.groupby("gene")["readCount"]
            .sum()
            .reset_index()
            .assign(srr=srr)
        )

        all_counts.append(gene_counts)

    if not all_counts:
        raise ValueError(f"No data found for gene={gene} in {results_dir}")
    
    count_matrix = (
        pd.concat(all_counts)
        .pivot(index="gene", columns="srr", values="readCount")
        .fillna(0)
    )   

    return count_matrix

def clr_transform(freq_matrix: pd.DataFrame):
    """
    A function that perform CLR-transformation on a frequency matrics.

    Takes one argument: the frequency matrix - a pandas DataFrame with genes as columns and samples as rows.

    Returns a CLR-transformed frequency matrix, a pandas DataFrame
    """
    clr_freqs = freq_matrix.copy()
    SRRs = freq_matrix.index
    gene_cols = freq_matrix.columns

    nonzero_min = freq_matrix[freq_matrix > 0].min().min()
    delta = nonzero_min / 2

    clr_freqs = multi_replace(clr_freqs, delta = delta)

    clr_freqs = clr(clr_freqs)
    clr_freqs = pd.DataFrame(clr_freqs, index=SRRs, columns=gene_cols)

    return clr_freqs

def plot_clustermap(df: pd.DataFrame, gene_cols: list, title: str = "") -> None:
    """
    Plots a clustermap for CLR-transformed gene frequencies.

    - df: DataFrame with SRR, diagnosis, predicted_ifn_status and gene columns
    - gene_cols: list of gene column names
    - title: optional plot title
    """
    diagnosis_palette = {'H': 'greenyellow', 'MS': 'darkorange', 'SLE': 'm', 'CLE': 'thistle'}
    ifn_palette = {"High": "#d62728", "Low": "#1f77b4"}

    col_order = df.sort_values(["diagnosis", "predicted_ifn_status"])["SRR"].values

    col_colors = pd.DataFrame({
        "Diagnosis": df.set_index("SRR").loc[col_order, "diagnosis"].map(diagnosis_palette),
        "IFN status": df.set_index("SRR").loc[col_order, "predicted_ifn_status"].map(ifn_palette)
    }, index=col_order)

    data_matrix = df.set_index("SRR")[gene_cols].T[col_order]

    g = sns.clustermap(
        data=data_matrix,
        col_colors=col_colors,
        col_cluster=False,
        row_cluster=True,
        cmap="viridis",
        figsize=(15, 10),
        cbar_kws={"label": "CLR-transformed frequency"}
    )

    diag_patches = [mpatches.Patch(color=c, label=l) for l, c in diagnosis_palette.items()]
    ifn_patches = [mpatches.Patch(color=c, label=l) for l, c in ifn_palette.items()]

    leg1 = g.ax_heatmap.legend(
        handles=diag_patches, title="Diagnosis",
        loc="upper left", bbox_to_anchor=(1.15, 1.1),
        frameon=True
    )
    g.ax_heatmap.add_artist(leg1)

    g.ax_heatmap.legend(
        handles=ifn_patches, title="IFN status",
        loc="upper left", bbox_to_anchor=(1.15, 0.75),
        frameon=True
    )

    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")

    if title:
        plt.suptitle(title, fontweight="bold", y=1.02)

    plt.show()