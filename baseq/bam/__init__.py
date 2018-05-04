import click, re, os
import numpy as np
from jinja2 import Template
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

@cli.command(short_help="check_bam")
@click.option('--path', '-b', default='./', help='File Path')
@click.option('--chr', '-c', default='./', help='chromesome')
@click.option('--start', '-s', default='./', help='start')
@click.option('--end', '-e', default='./', help='end')
def infos(path, chr, start, end):
    bam = BAMFILE(path)
    bam.bam_stats()
    bam.length_distribution()

@cli.command(short_help="split reads not in region...")
@click.option('--path', '-d', default='./', help='File Path')
@click.option('--species', '-s', default='human', help='Species: human, mouse, zebrafish')
def split_reads_not_in_region(path, species):
    print('RNA-Seq')

import subprocess, re, os
import numpy as np
from baseq.utils.runcommand import run_it, run_generator

class BAMFILE:

    def __init__(self, path):
        lines = run_it("samtools view -H {}".format(path))
        tt = [re.split("\t|:", x) for x in lines]
        chrs = [(x[2], int(x[4])) for x in tt if len(x)>2 and x[1]=="SN"]

        self.chrs = {chr[0]:chr[1] for chr in chrs}
        self.index = os.path.abspath(path+".bai")
        self.path = path
        self.reads = 0
        self.match_bases = 0
        if not os.path.exists(self.index):
            self.index = ""

    def match_length(self):
        from baseq.bam.cigar import match_length
        cmd = Template("samtools view {{path}} | awk '{print $6}' | head -n 100000").render(path=self.path)
        lines = run_it(cmd)
        length = [match_length(line) for line in lines if 'M' in line]
        self.match_length = round(np.mean(length), 3)
        print("[info] Mean match length".format(self.match_length))

    def bam_stats(self):
        lines = run_it("samtools flagstat {}".format(self.path))
        infos = [x.split() for x in lines]
        total = int(infos[0][0])
        mapped = int(infos[4][0])
        self.total_reads = total
        self.mapped_reads = mapped
        print("[info] Total {}, Mapped {}, Ratio {}".format(total, mapped, round(mapped/total, 3)))

    def bin_mean_depth(self, chr, start, end):
        if chr not in self.chrs:
            print("[error] {} not exist in genome".format(chr))
        results = run_generator("samtools depth -r {}:{}-{} {}".format(chr, start, end, self.path))
        depth = []
        for line in results:
            depth.append(int(line.split()[2]))
        return round(float(sum(depth))/(end-start+1), 3)

    def bam_read_counts(self, bam, chr, star, end):
        pass