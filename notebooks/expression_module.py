import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, zscale
import gzip
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def prepare_count_matrix(filepath, mapping_file='../data/ensembl_gene_names.txt.gz',
                         filtering_file='../data/genes_filtering.tsv'):
    counts = pd.read_csv(filepath, sep='\t', compression='gzip')
    lengths = counts[['Geneid_clean', 'Length']].copy()
    counts_matrix = counts[(counts.sum(axis=1) > 0)]
    mapping_df = pd.read_csv(mapping_file, sep='\t', compression='gzip')
    id_to_symbol = dict(zip(mapping_df['Gene stable ID'], mapping_df['Gene name']))
    original_ids_clean = [id.split('.')[0] for id in counts_matrix.index]
    new_gene_names = [id_to_symbol.get(id, id) for id in original_ids_clean]
    counts_matrix.index = new_gene_names
    if counts_matrix.index.duplicated().sum():
        counts_matrix = counts_matrix.groupby(level=0).sum()
    genes_to_keep = pd.read_csv(filtering_file, sep="\t", header=None)[0].tolist()
    filtered_counts = counts_matrix[counts_matrix.index.isin(genes_to_keep)]
    symbol_to_id = {v: k for k, v in id_to_symbol.items()}
    filtered_ensembl_ids = []
    for gene in filtered_counts.index:
        ensembl_id = symbol_to_id.get(gene, None)
        if ensembl_id is None:
            filtered_ensembl_ids.append(gene)
        else:
            filtered_ensembl_ids.append(ensembl_id)
    filtered_lengths = pd.Series(index=filtered_counts.index, dtype=float)
    for i, (gene_symbol, ensembl_id) in enumerate(zip(filtered_counts.index, filtered_ensembl_ids)):
        if ensembl_id in lengths.index:
            filtered_lengths[gene_symbol] = lengths[ensembl_id]
        else:
            median_length = lengths.median()
            filtered_lengths[gene_symbol] = median_length
    print(f"Filtered counts matrix: {filtered_counts.shape[0]} genes, {filtered_counts.shape[1]} samples")
    tpm = pd.DataFrame(index=filtered_counts.index, columns=filtered_counts.columns)
    for sample in filtered_counts.columns:
        counts_s = filtered_counts[sample]
        rpk = (counts_s * 1e9) / (filtered_lengths.loc[counts_s.index] * counts_s.sum())
        tpm[sample] = (rpk / rpk.sum()) * 1e6
    log2_tpm = np.log2(tpm + 1)
    print("Log2 TPM table prepared")
    return log2_tpm

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
            print(f"'{sig_name}': {len(available)}/{len(genes)} генов")
    
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

def run_deseq2(counts_df, metadata_df,
               condition_col='predicted_ifn_status',
               ref_condition='Low'):
    metadata_df[condition_col] = metadata_df[condition_col].astype('category')
    categories = metadata_df[condition_col].cat.categories.tolist()
    categories.remove(ref_condition)
    new_categories = [ref_condition] + categories
    metadata_df[condition_col] = metadata_df[condition_col].cat.reorder_categories(new_categories)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design_factors=condition_col,
        fit_type='parametric'
    )
    dds.deseq2()
    stat_res = DeseqStats(dds)
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
    return plt

def plot_boxplot_stats(log2_tpm, metadata, gene_list,
                       condition_col='diagnosis'):
    available_genes = [g for g in gene_list if g in log2_tpm.index]
    expr_data = []
    for gene in available_genes:
        for _, row in metadata.iterrows():
            srr = row['SRR']
            if srr in log2_tpm.columns:
                expr_data.append({
                    'SRR': srr,
                    'diagnosis': row['diagnosis'],
                    'IFN_status': row['predicted_ifn_status'],
                    'gene': gene,
                    'expression': log2_tpm.loc[gene, srr]
                })
    expr_df = pd.DataFrame(expr_data)
    fig, axes = plt.subplots(2, 4, figsize=(20, 5))
    axes = axes.flatten()
    for i, gene in enumerate(available_genes):
        ax = axes[i]
        subset = expr_df[expr_df['gene'] == gene]
        sns.boxplot(data=subset, x='diagnosis', y='expression', hue='IFN_status',
                    palette={'High': '#d62728', 'Low': '#1f77b4'}, ax=ax)
        sns.stripplot(data=subset, x='diagnosis', y='expression', hue='IFN_status',
                      palette={'High': '#d62728', 'Low': '#1f77b4'}, 
                      dodge=True, alpha=0.2, size=2, ax=ax)
        y_max = subset['expression'].max()
        diagnoses_gene = subset['diagnosis'].unique()
        for j, diag in enumerate(diagnoses_gene):
            high_vals = subset[(subset['diagnosis'] == diag) & (subset['IFN_status'] == 'High')]['expression']
            low_vals = subset[(subset['diagnosis'] == diag) & (subset['IFN_status'] == 'Low')]['expression']
            if len(high_vals) > 0 and len(low_vals) > 0:
                stat, pval = mannwhitneyu(high_vals, low_vals)
                if pval < 0.0001:
                    star = '****'
                elif pval < 0.001:
                    star = '***'
                elif pval < 0.01:
                    star = '**'
                elif pval < 0.05:
                    star = '*'
                else:
                    star = 'ns'
                if star != 'ns':
                    x1 = j - 0.2
                    x2 = j + 0.2
                    y_pos = y_max + 0.3
                    ax.plot([x1, x1, x2, x2], [y_pos - 0.1, y_pos, y_pos, y_pos - 0.1], 
                            color='black', linewidth=0.8)
                    ax.text((x1 + x2) / 2, y_pos + 0.05, star, ha='center', va='bottom', 
                            fontsize=8, fontweight='bold')
        ax.set_title(gene, fontsize=11)
        ax.set_xlabel('')
        ax.set_ylabel('log2(TPM+1)')
        ax.tick_params(axis='x', rotation=45)
        ax.legend_.remove()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc='upper right', title='IFN Status')
    for i in range(len(available_genes), len(axes)):
        axes[i].set_visible(False)
    plt.suptitle('IFN Signature Genes Expression by Diagnosis and IFN Status', 
                 fontsize=12, y=1.02)
    plt.tight_layout()
    plt.show()
    return fig, axes
    
