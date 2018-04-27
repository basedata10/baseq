import sys, os
from basematic.mgt.check import isExist

def QC_enrich(bampath, interval, outpath, fastMode=True):
    """
    Generate the
    :param bampath: input bam
    :param interval: genome regions
    :param outpath: output
    :param fastMode:
    :return:
        A density figure
    """
    import subprocess
    import pandas as pd
    import matplotlib as mpl
    mpl.use('Agg')

    path = bampath
    isExist(interval, "Interval bed file")
    isExist(bampath, "Bam file")

    if fastMode:
        print("Generating the bed file for chr15...")
        cmd_bed_chr1 = """less {} | awk '$1=="15" || $1=="chr15"' >{}.chr15 """.format(interval, interval)
        subprocess.check_call(cmd_bed_chr1, shell=True)
        interval = interval + ".chr15"

    depth_all_file = "{}.depth.all.txt".format(path)
    depth_bed_file = "{}.depth.bed.txt".format(path)

    try:
        #Depth for all position:
        if os.path.exists(depth_all_file):
            print("[info] {} exists, skip...".format(depth_all_file))
        else:
            cmd = "samtools depth -r 15 {} | awk '{{print $3}}' > {}.depth.all.txt".format(path, path)
            print("[info] Calculate the depth of bases...")
            print("[command] {}".format(cmd))
            subprocess.check_call(cmd, shell=True)

        #Depth for interval regions:
        if os.path.exists(depth_bed_file):
            print("[info] {} exists, skip...".format(depth_bed_file))
        else:
            cmd = "samtools depth -a -b {} {} | awk '{{print $3}}' > {}.depth.bed.txt".format(interval, path, path)
            print("[info] Calculate the depth of bases on the interval ...")
            print("[command] {}".format(cmd))
            subprocess.check_call(cmd, shell=True)
    except:
        sys.exit("[error] Failed to Execute the QC of Enrichment, please check samtools, bam or interval file")

    all = pd.read_table(depth_all_file)
    bed = pd.read_table(depth_bed_file)
    depth_all = sum(all.iloc[:,0])
    depth_bed = sum(bed.iloc[:,0])

    print(bed.iloc[:,0].mean(), bed.iloc[:,0].median(), depth_bed/depth_all)

    # ax = bed.iloc[10000:100000,0].plot.kde()
    # ax.set_xlim(0, 300)
    # fig = ax.get_figure()
    # fig.savefig("./figure.png")


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

