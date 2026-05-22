import pandas as pd
import numpy as np
import glob
import re

from dataclasses import dataclass

import argparse
import os

@dataclass
class Report:
    """
    Trust4/MiXCR Report DataClass

    Arguments:
    - report: pandas DataFrame with trust4 or MiXCR report table

    Attributes:
    - report: the same DataFrame
    - counts: pd.Series - the '#counts' column of the of a trust4 report (or 'readCount' column if using MiXCR)
    - total_reads: sum of all counts
    - n_obs: number of found clones
    """
    report: pd.DataFrame
    counts_col: str
    total_reads: int = None
    n_obs: int = None

    shannon_index = None
    chao1_index = None
    clonality = None
    simpson_index = None

    def __post_init__(self):
        self.total_reads = self.report[self.counts_col].sum()
        self.n_obs = len(self.report[self.counts_col])

    def calculate_shannon_index(self):
        """
        Calculates Shannon diversity index

        Returns a float
        """
        rel_abund = self.report[self.counts_col] / self.total_reads
        self.shannon_index = -np.sum(rel_abund * np.log(rel_abund))

        return self.shannon_index

    def calculate_chao1_index(self):
        """
        Calculates chao1 index for total clonotype richness

        Returns a float
        """
        f1 = (self.report[self.counts_col] == 1).sum()
        f2 = (self.report[self.counts_col] == 2).sum()

        if f2 > 0:
            self.chao1_index = self.n_obs + (f1 ** 2) / (2 * f2)
        else:
            self.chao1_index = self.n_obs + (f1 * (f1 - 1)) / 2

        return self.chao1_index

    def calculate_clonality(self):
        """
        Calculates clonality of a samples based on the given shannon index

        Returns a float
        """
        if self.shannon_index is None:
            self.shannon_index = self.calculate_shannon_index()

        if self.n_obs > 1:
            pielou = self.shannon_index / np.log(self.n_obs)
            self.clonality = 1 - pielou
        else:
            pielou = 0
            self.clonality = 1

        return self.clonality
    
    def calculate_simpson(self):
        """
        Calculates simpson index of a sample.

        Returns a float
        """
        counts = np.asarray(self.report[self.counts_col], dtype=float)

        counts = counts[counts > 0]
        p = counts / counts.sum()

        self.simpson_index = 1 - np.sum(p ** 2)
        return self.simpson_index
        
    def calculate_metrics(self):
        """
        Calculates all diversity metrics given above.

        Return a dictionaty with informative metrics.
        """
        shannon = self.calculate_shannon_index()

        self.div_metrics = {
            "Total_reads": self.total_reads,
            "N_obs": self.n_obs,
            "Shannon": self.calculate_shannon_index(),
            "Chao1": self.calculate_chao1_index(),
            "Clonality": self.calculate_clonality(),
            "Simpson": self.calculate_simpson()
        }

        return self.div_metrics
        
