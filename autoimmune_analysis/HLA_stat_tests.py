import sys
sys.path.append('..')
from typing import Any, Dict, Tuple

import pandas as pd
import numpy as np

from loguru import logger
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import Table2x2


logger.remove()
custom_format = (
    "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
    "<level>{level: <8}</level> \n"
    "<level>{message}</level> \n"
)
logger.add(sys.stderr, format=custom_format, colorize=True)


def make_tables(
    df: pd.DataFrame,
    diagnosis: str,
    test: str = 'fisher',
    min_carriers: int = 3,
) -> Dict[str, np.ndarray]:
    """Creating a contingency table for statistical tests"""

    logger.info(f"PROCESSING {diagnosis}-SAMPLES FOR {test.upper()} TEST")

    low = df.loc[df['IFN_status'] == 'Low', 'n_samples']
    high = df.loc[df['IFN_status'] == 'High', 'n_samples']

    total_low = low.iat[0] if not low.empty else 0
    total_high = high.iat[0] if not high.empty else 0

    tables = {}

    pivot = df.pivot_table(
        index=['allele'],
        columns='IFN_status',
        values=['carriers', 'allele_copies'],
        aggfunc='first',
        sort=False,
    )

    n = 2
    alleles_to_skip = []
    for allele, row in pivot.iterrows():
        carriers_high = row[('carriers', 'High')]
        carriers_low = row[('carriers', 'Low')]

        carriers_high = 0 if pd.isna(carriers_high) else carriers_high
        carriers_low = 0 if pd.isna(carriers_low) else carriers_low

        if carriers_high + carriers_low < min_carriers:
            alleles_to_skip.append(allele)
            continue

        high_positive = row[('allele_copies', 'High')]
        low_positive = row[('allele_copies', 'Low')]

        high_positive = 0 if pd.isna(high_positive) else high_positive
        low_positive = 0 if pd.isna(low_positive) else low_positive

        high_negative = total_high * n - high_positive
        low_negative = total_low * n - low_positive

        table = np.array([[high_positive, low_positive], [high_negative, low_negative]])

        table = table.astype(int)

        tables[allele] = table

    logger.info(
        f"{len(alleles_to_skip)} alleles were skipped "
        f"due to low carrier count (< {min_carriers}):\n"
        f"{', '.join(alleles_to_skip)}"
    )

    return tables


def apply_haldane_anscombe_correction(table: np.ndarray) -> np.ndarray:
    if np.any(table == 0):
        return table + 0.5
    return table

def count_odds_ratio(corrected_table: np.ndarray) -> float:
    OR = (
        corrected_table[0, 0] * corrected_table[1, 1]
    ) / (
        corrected_table[0, 1] * corrected_table[1, 0]
    )
    return OR


def make_row(
    diagnosis: str,
    allele: str,
    stat: float,
    pval: float,
    contingency_table: Table2x2,
    test: str,
) -> Dict[str, Any]:
    
    lower, upper = contingency_table.oddsratio_confint()
    return {
        'diagnosis': diagnosis,
        'gene': allele.split('*')[0],
        'allele': allele,
        'test': test,
        'pval': pval,
        'OR': stat,
        'CI_lower': lower,
        'CI_upper': upper,
    }


def get_adjusted_pvals(df_stat: pd.DataFrame) -> None:

    # global correction
    df_stat['p_adj'] = multipletests(
        df_stat['pval'],
        method='fdr_bh'
    )[1]

    # correction within each gene
    adjusted_series = []

    grouped = df_stat.groupby(
        ['diagnosis', 'gene'],
        sort=False
    )

    for _, group in grouped:

        if len(group) == 1:
            corrected = group['pval']

        else:
            corrected = pd.Series(
                multipletests(group['pval'], method='fdr_bh')[1],
                index=group.index
            )

        adjusted_series.append(corrected)

    df_stat['p_adj_by_gene'] = pd.concat(adjusted_series).sort_index()


def get_fisher_test_results(
    allele_counts: pd.DataFrame,
    min_carriers: int = 3,
) -> pd.DataFrame:

    stats = []

    for diagnosis, sub in allele_counts.groupby('diagnosis'):
        tables = make_tables(sub, diagnosis, test='fisher', min_carriers=min_carriers)

        rows = []
        for allele, table in tables.items():
            OR, pval = fisher_exact(table)
            corrected_table = apply_haldane_anscombe_correction(table)
            if np.any(table == 0):
                OR = count_odds_ratio(corrected_table)

            contingency_table = Table2x2(corrected_table)

            rows.append(
                make_row(diagnosis, allele, OR, pval, contingency_table, test='fisher')
            )

        df_stat = pd.DataFrame(rows)
        get_adjusted_pvals(df_stat)
        stats.append(df_stat)

    return pd.concat(stats, ignore_index=True)


def choose_test(table: np.ndarray) -> Tuple[str, float, float, Table2x2]:

    row_sums = table.sum(axis=1)
    col_sums = table.sum(axis=0)

    is_fisher = (row_sums == 0).any() or (col_sums == 0).any()

    if not is_fisher:
        chi2, pval, _, expected = chi2_contingency(table)
        if (expected < 5).any():
            is_fisher = True

    corrected_table = apply_haldane_anscombe_correction(table)

    if is_fisher:
        method = 'fisher'
        OR, pval = fisher_exact(table)
        if np.any(table == 0):
            OR = count_odds_ratio(corrected_table)
        stat = OR

    else:
        method = 'chi2'
        stat = chi2

    contingency_table = Table2x2(corrected_table)

    return method, stat, pval, contingency_table


def get_chi_test_results(
    allele_counts: pd.DataFrame,
    min_carriers: int = 3,
) -> pd.DataFrame:

    stats = []

    for diagnosis, sub in allele_counts.groupby('diagnosis'):
        tables = make_tables(sub, diagnosis, test='chi2', min_carriers=min_carriers)

        rows = []

        alleles_counted_by_fisher = []
        for allele, table in tables.items():
            method, stat, pval, contingency_table = choose_test(table)

            if method == 'fisher':
                alleles_counted_by_fisher.append(allele)

            rows.append(
                make_row(diagnosis, allele, stat, pval, contingency_table, test=method)
            )

        df_stat = pd.DataFrame(rows)

        get_adjusted_pvals(df_stat)

        stats.append(df_stat)
        logger.info(
            f"The Fisher-test instead of chi2-test was applied to {len(alleles_counted_by_fisher)} alleles "
            f"due to the low expected counts (< 5) or zero expected frequencies:\n"
            f"{', '.join(alleles_counted_by_fisher)}"
        )

    return pd.concat(stats, ignore_index=True)
