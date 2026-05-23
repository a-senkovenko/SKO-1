import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, zscore, f_oneway, tukey_hsd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def prepare_count_matrix(filepath, mapping_file='../data/genes_ensembl_length.csv',):
    counts = pd.read_csv(filepath, sep=',', index_col=0)
    print(f'Loaded counts: {counts.shape[0]} genes, {counts.shape[1]} samples')
    mapping_df = pd.read_csv(mapping_file, sep=',')
    lengths = mapping_df.set_index('Geneid')['Length']
    genes_to_keep = lengths.index.tolist()
    mask = counts.index.isin(genes_to_keep)
    filtered_counts = counts.loc[mask]
    print(f"Filtered counts matrix: {filtered_counts.shape[0]} genes, {filtered_counts.shape[1]} samples")
    tpm = pd.DataFrame(index=filtered_counts.index, columns=filtered_counts.columns)
    for sample in filtered_counts.columns:
        counts_s = filtered_counts[sample]
        rpk = (counts_s * 1e9) / (lengths.loc[counts_s.index] * counts_s.sum())
        tpm[sample] = (rpk / rpk.sum()) * 1e6
    log2_tpm = np.log2(tpm + 1)
    print("Log2 TPM table prepared")
    return log2_tpm, filtered_counts

def parse_gmt(gmt_file):
    signatures = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            signature_name = parts[0]
            description = parts[1]
            genes = parts[2:]
            signatures[signature_name] = genes
    return signatures

def calculate_signature_scores(log2_tpm, signatures):
    zscore_matrix = log2_tpm.apply(zscore, axis=1)
    signature_scores = {}
    for sig_name, genes in signatures.items():
        available = [g for g in genes if g in zscore_matrix.index]
        if len(available) > 0:
            signature_scores[sig_name] = zscore_matrix.loc[available].mean(axis=0)
            print(f"'{sig_name}': {len(available)}/{len(genes)} genes")
    if not signature_scores:
        return pd.DataFrame(), pd.Series()
    scores_df = pd.DataFrame(signature_scores).T
    return scores_df

def plot_signatures_heatmap(scores_df, dataset_name, diagnoses_list=None, figsize=(12, 8)):
    plt.figure(figsize=figsize)
    if diagnoses_list is not None and len(diagnoses_list) == scores_df.shape[1]:
        sort_idx = np.argsort(diagnoses_list)
        scores_sorted = scores_df.iloc[:, sort_idx]
        diagnoses_sorted = [diagnoses_list[i] for i in sort_idx]
        ax = sns.heatmap(scores_sorted, 
                         cmap='RdBu_r', 
                         center=0,
                         xticklabels=False,
                         yticklabels=True,
                         cbar_kws={'label': 'Signature Score'},
                         vmin=-2, vmax=2)
        color_map = {'SLE': 'red', 'H': 'green', 'MS': 'blue', 'CLE': 'purple'}
        for i, diag in enumerate(diagnoses_sorted):
            ax.plot(i + 0.5, -0.5, 's', markersize=6, color=color_map.get(diag, 'gray'), 
                    transform=ax.transData, clip_on=False)  
    plt.title(f'{dataset_name}: Signature Scores\n{scores_df.shape[0]} signatures, {scores_df.shape[1]} samples',
              y=1.05)
    plt.ylabel('Signatures')
    plt.xlabel('Samples')
    plt.tight_layout()
    plt.show()
    return plt

def run_deseq2(counts_df, metadata_df):
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design_factors='Predicted_IFN_status',
        ref_level='Low',
        fit_type='parametric'
    )
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=['Predicted_IFN_status', 'High', 'Low'])
    stat_res.summary()
    res = stat_res.results_df
    return res

def plot_volcano(res, padj_threshold=0.05, lfc_threshold=1):
    res['-log10_padj'] = -np.log10(res['padj'].clip(lower=1e-300))
    sig_up = (res['padj'] < padj_threshold) & (res['log2FoldChange'] > lfc_threshold)
    sig_down = (res['padj'] < padj_threshold) & (res['log2FoldChange'] < -lfc_threshold)
    not_sig = ~(sig_up | sig_down)
    plt.figure(figsize=(10, 8))
    plt.scatter(res[not_sig]['log2FoldChange'], res[not_sig]['-log10_padj'], 
                c='lightgray', alpha=0.5, s=20)
    plt.scatter(res[sig_up]['log2FoldChange'], res[sig_up]['-log10_padj'], 
                c='#E63946', alpha=0.7, s=30, label=f'Up in High ({sig_up.sum()})')
    plt.scatter(res[sig_down]['log2FoldChange'], res[sig_down]['-log10_padj'], 
                c='#457B9D', alpha=0.7, s=30, label=f'Down in High ({sig_down.sum()})')
    top20_genes = res.nsmallest(20, 'padj').index.tolist()
    for gene in top20_genes:
        x = res.loc[gene, 'log2FoldChange']
        y = res.loc[gene, '-log10_padj']
        plt.annotate(gene, (x, y), fontsize=8, alpha=0.8, 
                     xytext=(5, 5), textcoords='offset points',
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
    plt.axhline(-np.log10(0.05), linestyle='--', color='gray', label='padj = 0.05')
    plt.axvline(-1, linestyle='--', color='gray')
    plt.axvline(1, linestyle='--', color='gray')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10(padj)')
    plt.title('DESeq2: IFN-High vs IFN-Low (Healthy)')
    plt.legend()
    plt.tight_layout()
    plt.show()

def get_significance_stars(pval):
    if pval < 0.0001:
        return '****'
    elif pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    else:
        return 'ns'
    
def add_significance_bars(ax, x1, x2, y_pos, pval, bar_height=0.1):
    star = get_significance_stars(pval)
    if star != 'ns':
        ax.plot([x1, x1, x2, x2],
                [y_pos - bar_height, y_pos, y_pos, y_pos - bar_height],
                color='black', linewidth=0.8)
        ax.text((x1 + x2) / 2, y_pos + 0.05, star,
                ha='center', va='bottom', fontsize=8, fontweight='bold')

def plot_boxplot_genes(log2_tpm, metadata, gene_list,
                       condition_col='diagnosis'):
    group_col = 'Predicted_IFN_status'
    available_genes = [g for g in gene_list if g in log2_tpm.index]
    expr_data = []
    for gene in available_genes:
        for _, row in metadata.iterrows():
            srr = row['SRR']
            if srr in log2_tpm.columns:
                expr_data.append({
                    'SRR': srr,
                    condition_col: row[condition_col],
                    group_col: row[group_col],
                    'gene': gene,
                    'expression': log2_tpm.loc[gene, srr]
                })
    expr_df = pd.DataFrame(expr_data)
    n_cols = 3
    n_genes = len(available_genes)
    n_rows = (n_genes + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 10))
    axes = axes.flatten() if n_rows > 1 else [axes]
    for i, gene in enumerate(available_genes):
        ax = axes[i]
        subset = expr_df[expr_df['gene'] == gene]
        sns.boxplot(data=subset, x=condition_col, y='expression', hue=group_col,
                    palette={'High': '#d62728', 'Low': '#1f77b4'}, ax=ax)
        sns.stripplot(data=subset, x=condition_col, y='expression', hue=group_col,
                      palette={'High': '#d62728', 'Low': '#1f77b4'}, 
                      dodge=True, alpha=0.2, size=2, ax=ax)
        y_max = subset['expression'].max()
        conditions = subset[condition_col].unique()
        for j, diag in enumerate(conditions):
            high_vals = subset[(subset[condition_col] == diag) & 
                                (subset[group_col] == 'High')]['expression'].dropna()
            low_vals = subset[(subset[condition_col] == diag) &
                                (subset[group_col] == 'Low')]['expression'].dropna()
            if len(high_vals) > 0 and len(low_vals) > 0:
                stat, pval = mannwhitneyu(high_vals, low_vals)
                x1 = j - 0.2
                x2 = j + 0.2
                y_pos = y_max + 0.3
                add_significance_bars(ax, x1, x2, y_pos, pval)
        ax.set_title(gene, fontsize=11)
        ax.set_xlabel('')
        ax.set_ylabel('log2(TPM+1)')
        ax.tick_params(axis='x', rotation=45)
        ax.legend_.remove()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc='upper right', title='IFN Status')
    for i in range(len(available_genes), len(axes)):
        axes[i].set_visible(False)
    plt.suptitle(f'IFN Signature Genes Expression by {condition_col} and IFN Status', 
                 fontsize=12, y=1.02)
    plt.tight_layout()
    plt.show()
    
def calculate_ifn_score(metadata, log2_tpm, ifn_genes):
    metadata['IFN_score'] = np.nan
    for gse in metadata['GSE'].unique():
        gse_samples = metadata[metadata['GSE'] == gse]['SRR'].tolist()
        gse_samples = [s for s in gse_samples if s in log2_tpm.columns]
        available_genes = [g for g in ifn_genes if g in log2_tpm.index]
        gse_data = log2_tpm.loc[available_genes, gse_samples].copy()
        gse_zscored = gse_data.T.apply(zscore)
        ifn_scores = gse_zscored.mean(axis=1)
        for sample in gse_samples:
            metadata.loc[metadata['SRR'] == sample, 'IFN_score'] = ifn_scores[sample]
    return metadata

def plot_boxplot_score(metadata, condition_col='diagnosis', test='mannwhitney'):
    plt.figure(figsize=(7, 4))
    ax = plt.gca()
    sns.boxplot(data=metadata, x=condition_col, y='IFN_score', hue='Predicted_IFN_status', 
                palette={'High': '#d62728', 'Low': '#1f77b4'}, ax=ax)
    sns.stripplot(data=metadata, x=condition_col, y='IFN_score', hue='Predicted_IFN_status',
                  palette={'High': '#d62728', 'Low': '#1f77b4'}, dodge=True, alpha=0.3, size=3,
                  legend=False, ax=ax)
    diagnoses = metadata[condition_col].unique()
    y_max = metadata['IFN_score'].max()
    y_min = metadata['IFN_score'].min()
    y_range = y_max - y_min
    for i, diag in enumerate(diagnoses):
        high_vals = metadata[(metadata[condition_col] == diag) & 
                            (metadata['Predicted_IFN_status'] == 'High')]['IFN_score'].dropna()
        low_vals = metadata[(metadata[condition_col] == diag) & 
                            (metadata['Predicted_IFN_status'] == 'Low')]['IFN_score'].dropna()
        stat, pval = mannwhitneyu(high_vals, low_vals)
        x1 = i - 0.2
        x2 = i + 0.2
        y_pos = y_max + (y_range * 0.05)
        add_significance_bars(ax, x1, x2, y_pos, pval)
    plt.axhline(0, linestyle='--', color='black', alpha=0.5)
    plt.xlabel('Diagnosis', fontsize=12)
    plt.ylabel('IFN Score', fontsize=12)
    plt.title(f'IFN Score by {condition_col} and IFN Status', fontsize=14)
    plt.legend(title='IFN Status', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45)
    plt.ylim(y_min - y_range*0.05, y_max + y_range*0.4)
    plt.tight_layout()
    plt.show()
