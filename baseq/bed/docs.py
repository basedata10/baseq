import subprocess, re, os
from baseq.utils.runcommand import run_it, run_generator
import pandas as pd

"""
baseq dev bed ./bed
"""

class BEDFILE:
    def __init__(self, path):
        self.bed = pd.read_table(path, usecols=range(3), names=['chr', 'start', 'end'], comment='#')
        self.stats()
        self.sampling(100)

    def stats(self):
        lengths = []
        for index, row in self.bed.iterrows():
            length = row['end'] - row['start']
            lengths.append(length)
        self.length = sum(lengths)
        self.counts = len(lengths)
        print("[info] Intervels {} Length {}.".format(self.counts, self.length))

    def sampling(self, numbers=100):
        df_s = self.bed.sample(n=numbers)
        return df_s.values.tolist()

    def merge(self):
        pass