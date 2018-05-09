from baseq.mgt.config import get_config
import pandas as pd
import re

def read_gencode(genome, type):
    path = get_config("RNA_ref_"+genome, "gencode")
    df = pd.read_table(path, comment="#", header=None, names=["chr", "source", "type", "start", "end", "dot", "strand", "dot2", "infos"])
    df = df.loc[df.type==type]
    infos = [re.split("\"", x) for x in df.infos.tolist()]
    df["gene"] = [x[9] for x in infos]
    df["transc"] = [x[3] for x in infos]
    df["exonid"] = [re.split("\s+|;", x[16])[3] for x in infos]
    df = df.drop(["infos", "dot", "dot2", "source", "type"], axis=1)
    return df
