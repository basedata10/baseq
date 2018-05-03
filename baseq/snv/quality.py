import sys, os
from baseq.mgt.check import isExist



def enrich_saturation(bampath, interval, outpath):
    """ Check the unique number of reads with different levels of sequencing depth...
    """
    
    print("enrich_saturation analysis")

def parse_picard_wes_matrics(path):
    with open(path, 'r') as file:
        lines = file.readlines()
    names = lines[6].split()
    datas = lines[7].split()

    #stats
    res = {}
    res['Interval Size'] = int(datas[0])
    res['Mean Coverage'] = float(datas[1])
    res['PCT_30X'] = float(datas[18])

    #density plot...
    depth = [int(line.split()[0]) for line in lines[11:-1]]
    coverage = [int(line.split()[1]) for line in lines[11:-1]]
    coverage = [x*len(coverage)/sum(coverage) for x in coverage]
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.plot(depth, coverage)
    plt.ylim((0, 10))
    plt.savefig("./density.png")
    return res

