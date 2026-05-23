import sys

import pandas as pd
from loguru import logger

from autoimmune_analysis import HLA_stat_tests
from autoimmune_analysis import HLA_pie_charts


logger.remove()
custom_format = (
    "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
    "<level>{level: <8}</level> \n"
    "<level>{message}</level> \n"
)
logger.add(sys.stderr, format=custom_format, colorize=True)


def check_genes_type(genes):
    """
    Checks the format of the input genes
    and returns a list of gene names.
    """
    if isinstance(genes, str):
        genes_list = [g.strip() for g in genes.split(',')]
    elif isinstance(genes, list):
        genes_list = genes
    else:
        raise ValueError(
            "Genes should be either a comma-separated string or a list of strings."
        )
    return genes_list


def drop_low_observation_columns(df: pd.DataFrame, min_obs: int) -> pd.DataFrame:
    """
    Filters out columns from the DataFrame
    that have less than min_obs non-NA observations.
    """
    genes_to_delete = []
    for c in df.columns:
        if df[c].count() < min_obs:
            genes_to_delete.append(c)
    df = df.drop(columns=genes_to_delete)
    if genes_to_delete:
        logger.info(
            f"Columns {genes_to_delete} were deleted, "
            f"because less than {min_obs} observations were in them."
        )
    return df


def make_allele_table(df: pd.DataFrame, genes: list) -> pd.DataFrame:
    """
    Transforms the original DataFrame into a long format DataFrame with
    one row per one allele within each diagnosis and IFN status group.
    """

    def make_tmp_df(col):
        tmp = df[['SRR', 'diagnosis', 'IFN_status', col]].copy()
        tmp.columns = ['sample', 'diagnosis', 'IFN_status', 'allele']
        tmp['gene'] = gene
        tmp['copy'] = 1 if col.endswith('1') else 2
        return tmp

    rows = []

    for gene in genes:
        col1, col2 = f'{gene}1', f'{gene}2'
        if col1 not in df.columns:
            continue

        tmp1 = make_tmp_df(col1)
        tmp2 = make_tmp_df(col2)

        rows.extend([tmp1, tmp2])

    long_df = pd.concat(rows, ignore_index=True)
    long_df = long_df.dropna(subset=['allele'])

    long_df = long_df.sort_values(
        ['sample', 'gene'], ascending=[True, True], ignore_index=True
    )

    long_df['allele'] = long_df['allele'].str.split(':').str[:2].str.join(':')

    return long_df


class HLA:
    """
    A class for analyzing aracasHLA data, including counting allele frequencies
    and performing statistical tests.

    The names of alleles are converted to a two-field format.
    """

    def __init__(self, df: pd.DataFrame, filtration: bool = False, min_obs: int = 0):
        """
        Initializes the HLA class with the given DataFrame and parameters.

        Args:
            df (pd.DataFrame): The input DataFrame containing the results
                of genotyping by arcasHLA in tsv-format.
            filtration (bool, optional): Whether to filter the DataFrame
                based on the minimum number of non-NA observations per gene.
                Defaults to False.
            min_obs (int, optional): The minimum number of observations
                required for a column to be included.
                Defaults to 0.
        """
        if filtration:
            self.df = drop_low_observation_columns(df, min_obs)
        else:
            self.df = df

        genes = (
            'A, B, C, DMA, DMB, DOA, DOB, '
            'DPA1, DPB1, DQA1, DQB1, DRA, DRB1, DRB3, DRB5, '
            'E, F, G, H, J, K, L'
        )

        self.genes_list = [g.strip() for g in genes.split(',')]

        self.long_df = make_allele_table(self.df, self.genes_list)

        self.allele_counts = None
        self.fisher_test_results = None
        self.fisher_test_results_2 = None
        self.chi2_test_results = None

    def count_alleles(self, genes: str | list[str] | None = None) -> pd.DataFrame:
        """
        Counts the number of carriers and allele copies for each allele,
        as well as their frequencies within each diagnosis and IFN status group.

        Args:
            genes (str | list, optional): A comma-separated string or
                a list of gene names to include in the count.
                If None, all genes in the DataFrame will be included.
                Defaults to None.

        Returns:
            pd.DataFrame: A DataFrame containing the allele counts
                and frequencies for each allele within
                each diagnosis and IFN status group.
        """
        if genes is None:
            long_df = self.long_df.copy()
        else:
            genes = check_genes_type(genes)
            long_df = self.long_df[self.long_df['gene'].isin(genes)].copy()

        colnames_1 = ['diagnosis', 'gene', 'allele', 'IFN_status']

        carriers = (
            long_df.groupby(colnames_1)['sample'].nunique().reset_index(name='carriers')
        )

        def make_col_by_size(df, colnames, name):
            return df.groupby(colnames).size().reset_index(name=name)

        allele_counts = make_col_by_size(long_df, colnames_1, 'allele_copies')

        # ZYGOSITY
        colnames_2 = ['sample', 'diagnosis', 'gene', 'allele', 'IFN_status']

        dosage = make_col_by_size(long_df, colnames_2, 'dosage')

        homozygotes = make_col_by_size(dosage[dosage['dosage'] == 2], colnames_1, 'hom')
        heterozygotes = make_col_by_size(
            dosage[dosage['dosage'] == 1], colnames_1, 'het'
        )

        # MERGE
        def merge_cols(df, col, on, how):
            return df.merge(col, on=on, how=how)

        final_counts = merge_cols(carriers, allele_counts, colnames_1, 'outer')
        final_counts = merge_cols(final_counts, homozygotes, colnames_1, 'left')
        final_counts = merge_cols(final_counts, heterozygotes, colnames_1, 'left')
        final_counts = final_counts.fillna(0)
        final_counts[['het', 'hom']] = final_counts[['het', 'hom']].astype(int)

        # FREQUENCES
        n_samples = (
            long_df.groupby(['diagnosis', 'IFN_status'])['sample']
            .nunique()
            .reset_index(name='n_samples')
        )

        final_counts = final_counts.merge(n_samples, on=['diagnosis', 'IFN_status'])

        final_counts['carrier_freq'] = (
            final_counts['carriers'] / final_counts['n_samples']
        ).round(decimals=4)

        final_counts['allele_freq'] = (
            final_counts['allele_copies'] / (2 * final_counts['n_samples'])
        ).round(decimals=4)

        # SORT
        final_counts = final_counts.sort_values(
            ['diagnosis', 'gene', 'carriers', 'allele_copies'],
            ascending=[True, True, False, False],
        ).reset_index(drop=True)

        self.allele_counts = final_counts

        return final_counts

    def fisher_test(self, min_carriers: int = 3) -> pd.DataFrame:
        """
        Applies the Fisher's exact test to the allele counts
        for each diagnosis group.

        Args:
            min_carriers (int): Minimum number of carriers
                required for an allele to be included in the test.
        """
        if self.allele_counts is None:
            raise RuntimeError(
                "Run count_alleles() before fisher_test()."
            )
        
        self.fisher_test_results = HLA_stat_tests.get_fisher_test_results(
            self.allele_counts, min_carriers
        )

        return self.fisher_test_results

    def chi2_test(self, min_carriers: int = 3) -> pd.DataFrame:
        """
        Applies the chi-square test to the allele counts
        for each diagnosis group.

        Args:
            min_carriers (int): Minimum number of carriers
                required for an allele to be included in the test.
        """
        if self.allele_counts is None:
            raise RuntimeError(
                "Run count_alleles() before chi2_test()."
            )
        
        self.chi2_test_results = HLA_stat_tests.get_chi_test_results(self.allele_counts, min_carriers)

        return self.chi2_test_results

    def make_pie(
        self, genes: int | list[str] | None = None, cohort_name: str | None = None
    ) -> None:
        """
        Plots pie charts for specified HLA genes

        Args:
            genes (int or list): the gene or genes to display on the chart
                1 (int) - classic genes of the first class (A, B, C)
                2 (int) - classic genes of the second class (DP, DR, DQ)
                nonc(str) - nonclassic genes
                custom value (list) - a list of the one or several genes.
                Default - None.
            cohort_name (str): the name of the analyzed cohort
                to display in the plot title.

        Returns:
            None: The plot is shown with plt.show().

        """
        if self.allele_counts is None:
            raise RuntimeError(
                "Run count_alleles() before make_pie()."
            )
        
        classes_of_allele_genes = {
            1: ['A', 'B', 'C'],
            2: ['DPA1', 'DPB1', 'DRA', 'DRB1', 'DRB3', 'DRB5', 'DQA1', 'DQB1'],
            'nonc': ['DMA', 'DMB', 'DOA', 'DOB', 'E', 'F', 'G', 'H', 'J', 'K', 'L'],
        }

        if not genes:
            genes = self.allele_counts['gene'].unique().tolist()
        elif isinstance(genes, list):
            genes = genes
        elif genes not in classes_of_allele_genes.keys():
            raise ValueError(
                "`genes` must be int (1 or 2), str (nonc) or list."
            )
        else:
            genes = classes_of_allele_genes[genes]

        genes_from_counts = set(self.allele_counts['gene'].tolist())
        genes = [x for x in genes if x in genes_from_counts]

        pie_charts.make_pie(self.allele_counts, genes, cohort_name)
