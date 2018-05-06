import subprocess, re, os
from baseq.utils.runcommand import run_it, run_generator
import pandas as pd
import random
"""
baseq dev bed ./bed
"""

import click, os, sys
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

class BEDFILE:
    def __init__(self, path):
        self.bed = pd.read_table(path, usecols=range(3), names=['chr', 'start', 'end'], comment='@', converters={'chr':str})
        self.stats()

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

    def sample_split_files(self, lines=100, files=10):
        paths = []
        for x in range(files):
            path = "sample.{}.bed".format(x)
            paths.append(path)
            self.bed.sample(n=lines).to_csv(path, index=False, sep="\t", header=False)
        return paths