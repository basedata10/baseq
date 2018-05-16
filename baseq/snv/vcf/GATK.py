import pandas as pd
import re
import numpy as np

def vcf_stats(sample, vcfpath, min_depth=50):
    """
    Stats on the VCF from GATK
    ::
        vcf_stats("sample1", "path/to/vcf", min_depth=30)

    Return:
        A dict/json containing:
        Samplename/counts/mean_depth/GT_01/GT_11/MAF
        MAF is minor allel frequency.
    """

    print("[info] {} {}".format(sample, vcfpath))

    df = pd.read_table(vcfpath, comment="#", header=None,
                       names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sample'])

    infos = df['Sample'].tolist()
    info_splits = [re.split(":", x) for x in infos]
    df['DP'] = [int(x[2]) for x in info_splits]
    df['GT'] = [x[0] for x in info_splits]
    df['GQ'] = [int(x[2]) for x in info_splits]
    df['MAF'] = [round(min(z) / sum(z), 3) for z in [[int(y) for y in x[1].split(",")] for x in info_splits]]

    df = df.loc[df['DP']>=min_depth, :]
    df_mutation = pd.DataFrame(data=0, index=["A", "T", "C", "G"], columns=["A", "T", "C", "G"])

    for idx, row in df.iterrows():
        if row['REF'] in ["A", "T", "C", "G"] and row["ALT"] in ["A", "T", "C", "G"]:
            df_mutation.loc[row['REF'], row['ALT']] += 1

    MAF, bin_edges = np.histogram(df['MAF'], bins=50, range=(0, 0.5), normed=True)
    MAF = np.round(MAF, 3)
    print(MAF)

    stats = {
        "sample": sample,
        "counts": len(df['DP']),
        "mean_depth": round(sum(df['DP']) / len(df['DP']), 1),
        "GT_01": sum([1 for x in df['GT'] if x == "0/1"]),
        "GT_11": sum([1 for x in df['GT'] if x == "1/1"]),
        "MAF": MAF
    }

    return stats