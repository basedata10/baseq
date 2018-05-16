from baseq.mgt.config import get_config
import pandas as pd
import re, os, time

def read_gencode(genome, type):
    path = get_config("RNA_ref_"+genome, "gencode")
    path_hdf = path + ".hdf"
    if not os.path.exists(path_hdf):
        df = pd.read_table(path, comment="#", header=None, names=["chr", "source", "type", "start", "end", "dot", "strand", "dot2", "infos"])
        df.to_hdf(path_hdf, "gtf")
    else:
        df = pd.read_hdf(path_hdf, key="gtf")
    df = df.loc[df.type == type]
    infos = [re.split("\"", x) for x in df.infos.tolist()]
    df["gene"] = [x[9] for x in infos]
    df["transc"] = [x[3] for x in infos]
    df["exonid"] = [re.split("\s+|;", x[16])[3] for x in infos]
    df = df.drop(["infos", "dot", "dot2", "source", "type"], axis=1)
    return df

def get_isolated_genes(genome, min_Gap, minLength):
    r"""Get The Isolated Genes Which Do not Have Genes At upstream and downstream.
    :param genome: Genome Version Name like hg38.
    :param min_Gap: mininum distance to the upstream and downstream Gene.
    :param minLength: mininum length of the gene.
    :return: A int ...
    :rtype: list
    """
    path = get_config("RNA_ref_"+genome, "gencode")
    df = pd.read_table(path, comment="#", header=None, names=["chr", "source", "type", "start", "end", "dot", "strand", "dot2", "infos"])
    df = df.loc[df.type==type]
    infos = [re.split("\"", x) for x in df.infos.tolist()]
    df["gene"] = [x[9] for x in infos]
    df["transc"] = [x[3] for x in infos]
    df["exonid"] = [re.split("\s+|;", x[16])[3] for x in infos]
    print(df)