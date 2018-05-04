from baseq.utils.runcommand import run_it

def region_reads(path, chr, start, end):
    command = "samtools depth {}".format(path)
    pass

def region_bases(path, chr, start, end):
    command = "samtools depth {}".format(path)
    pass

def multi_region_reads(path, regions):
    """
    region should be defined: chr/start/end
    """
    command = "samtools depth {}".format(path)
    pass