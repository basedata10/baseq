import subprocess, re, os

def check_bam(path):
    """
    chromesome name and lenght;
    bai;
    size;
    """
    res = {"chr":"", "bai":""}
    from baseq.utils.runcommand import run_it
    lines = run_it("samtools view -H {}".format(path))
    tt = [re.split("\t|:", x) for x in lines]
    chrs = [(x[2], int(x[4])) for x in tt if len(x)>2 and x[1]=="SN"]
    res["chr"] = chrs
    bai = path+".bai"
    if os.path.exists(bai):
        res["bai"] = bai
    return res

def length_distribution(path):
    pass

def bam_stats(path):
    pass

def bam_depth(path, interval):
    pass

def bam_read_counts(bam, chr, star, end):
    pass