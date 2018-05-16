import re, os, time, sys
import numpy as np
from jinja2 import Template
from baseq.utils.runcommand import run_it, run_generator
from baseq.bed import BEDFILE
from multiprocessing.pool import ThreadPool
from baseq.bam.cigar import overlap

"""
reads_total/reads_mapped/reads_mapratio/
duplicate_reads/
match_length/match_bases/
region_length/region_mean_depth/region_total_bases/
"""

class BAMTYPE:
    """
    BAM File Handler, Based on samtools.
    While initiate, it read the path using samtools and will parse the headers.

    Usage:
        * Stats on enrichment quality.

    """
    def __init__(self, path, bedfile=""):
        if not os.path.exists(path):
            sys.exit("[error] Bam File {} Not Exists".format(path))
        lines = run_it("samtools view -H {}".format(path))
        #get the bam header
        header = [re.split("\t|:", x) for x in lines]
        chrs = [(x[2], int(x[4])) for x in header if len(x)>2 and x[1]=="SN"]
        #chromesome names
        self.chrs = {chr[0]:chr[1] for chr in chrs}
        self.path = path
        self.reads = 0

        #bam index file
        self.index = os.path.abspath(path+".bai")
        if not os.path.exists(self.index):
            self.index = ""

        #set bedfile
        if bedfile and os.path.exists(bedfile):
            self.bedfile = BEDFILE(bedfile)
        else:
            self.bedfile = None

    def get_columns(self, rows=10000, colIdx=6):
        """
        Read the bamfile using samtools, get the infors in the column<colIDx> and first <rows> of Rows.
        The columns of bam files are:

        #. header
        #. flags
        #. chromosome
        #. start
        #. mapping quality
        #. cigar

        The colIdx start from 1.
        ::
          BAMTYPE(path).get_columns(1000, 3)
          # ['chr1', 'chr1', ...]
        """
        cmd = Template("samtools view {{path}} | awk '{print ${{colIdx}}}' | head -n {{rows}}").render(path=self.path, colIdx=colIdx, rows=rows)
        lines = run_it(cmd)
        return lines

    def region_depth(self, chr, start, end, all=False):
        """
        Get the depth coverage of bases in the region.
        It will suitable for chromesome name like "chr1" and "1".

        :param all: Shall the bases with zero coverge be returned.

        Usage:
        ::
            BAMTYPE(path).region_depth("chr1", 1000, 2000, all=True)
            `return depth list [0,1,1,1,2,2,2,3,0]`
        """
        chr = str(chr)
        depth = []
        if chr in self.chrs:
            chr_real = chr
        elif chr.startswith("chr") and chr[3:] in self.chrs:
            chr_real = chr[3:]
        else:
            print("[error] {} not exist in genome".format(chr))
        if all:
            results = run_generator("samtools depth -a -r {}:{}-{} {}".format(chr_real, start, end, self.path))
        else:
            results = run_generator("samtools depth -r {}:{}-{} {}".format(chr_real, start, end, self.path))
        for line in results:
            depth.append(int(line.split()[2]))
        return depth

    def region_bed_depth(self, bedfile):
        """
        The depth of regions in an bed file.

        .. todo::

           This function is not suitable for sending spam e-mails.

        """
        results = run_generator("samtools depth -b {} {}".format(bedfile, self.path))
        depth = []
        for line in results:
            depth.append(int(line.split()[2]))
        return depth

    def stats_bam(self):
        """
        Read the bampath.stat, if not exists, perform the `samtools flagstat`
        The results will be:

        #. self.reads_total
        #. self.reads_mapped
        #. self.mapping_ratio
        """
        stat_file = self.path + ".stat"
        if os.path.exists(stat_file):
            print("[info] {} Exists".format(stat_file))
            with open(stat_file, "r") as file:
                lines = file.readlines()
        else:
            lines = run_it("samtools flagstat {}".format(self.path))
            try:
                with open(stat_file, "w") as file:
                    file.writelines("\n".join(lines))
            except:
                print("[info] Failed to write to {}".format(stat_file))

        infos = [x.split() for x in lines]
        total = int(infos[0][0])
        mapped = int(infos[4][0])
        self.reads_total = total
        self.reads_mapped = mapped
        self.mapping_ratio = round(float(mapped)/total, 3)
        print("[info] Reads Stats, Total {}, Mapped {}, Ratio {}".format(total, mapped, self.mapping_ratio))

    def stats_bases(self):
        """
        Stats on the mean match length for the top 100K reads in the bam file
        """
        from baseq.bam.cigar import match_length
        self.stats_bam()
        lines = self.get_columns(100000, 6)
        length = [match_length(line) for line in lines if 'M' in line]
        self.match_length = round(np.mean(length), 3)
        self.match_bases = self.match_length * self.reads_mapped
        print("[info] Mean matched read length: {}".format(self.match_length))

    def stats_duplicates(self):
        """
        Stats Duplication Rates from the top 1M reads;
        The duplication should be reflected in the flag
        """
        flags = self.get_columns(1000000, 2)
        dups = sum([1 for x in flags if int(x) & 1024])
        self.dup_ratio = round(float(dups)/1000000,3)
        print("[infos] Total {}, Dups {}, Ratio {}".format(1000000, dups, dups/1000000))

    def stats_regions(self):
        if self.bedfile:
            self.bedfile.stats()
        else:
            print("[info] No Bed File Exists")

    def stats_region_coverage(self, numbers=1000):
        """
        Check the enrichment quality.

        * Require a bedfile while initiating the class
        * Select <numbers> regions randomly
        * Use multithread pool to get the coverage depth of the regions
        * Stats on the ratio of 10X, 30X, 50X and 100X bases

        Usage:
        ::
            BAMTYPE("sample.bam", "panel.bed").stats_region_coverage(1000)
            The results will be save in object properies:
            self.mean_depth/self.pct_10X/..
        """
        intervals = self.bedfile.sampling(int(numbers))
        pool = ThreadPool(processes=10)
        results = []
        total_bases = 0

        for interval in intervals:
            total_bases += interval[2]-interval[1]+1
            results.append(pool.apply_async(self.region_depth, (interval[0], interval[1], interval[2])))

        pool.close()
        pool.join()

        region_depth = [r.get() for r in results]
        depth_flat = [item for sublist in region_depth for item in sublist]
        self.mean_depth = round(sum(depth_flat)/total_bases, 3)
        self.pct_10X = round(sum([1 for x in depth_flat if x>=10])/total_bases, 3)
        self.pct_30X = round(sum([1 for x in depth_flat if x>=30])/total_bases, 3)
        self.pct_50X = round(sum([1 for x in depth_flat if x>=50])/total_bases, 3)
        self.pct_100X = round(sum([1 for x in depth_flat if x>=100])/total_bases, 3)
        print("[info] StatBases {}, Mean {}, 10X {}, 30X {}, 50X {}, 100X {}".format(total_bases, self.mean_depth, self.pct_10X, self.pct_30X, self.pct_50X, self.pct_100X))

    def get_reads(self, chr, start, end):
        """
        Return The Reads that overlaps with region chrN:start-end.

        * Skip reads contains "N" cigar.
        """
        chr = str(chr)
        if chr in self.chrs:
            chr_real = chr
        elif chr.startswith("chr") and chr[3:] in self.chrs:
            chr_real = chr[3:]
        results = run_generator("samtools view {} {}:{}-{}".format(self.path, chr_real, start, end))
        reads = []
        for x in results:
            lines = x.split()
            flags = lines[5]
            #overlap_length = overlap(flags, lines[3], start, end)
            #if overlap_length<=30:
            if re.search("N", flags):
                print(flags)
                pass
            else:
                reads.append(lines[0:6])
        return reads

    def read_counts(self, chr, star, end):
        """
        Todo:
            * For module TODOs
        """
        pass