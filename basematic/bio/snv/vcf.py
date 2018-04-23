import os
import pandas as pd

def filter_annovar_result(inPath, outPath):
    df = pd.read_csv(inPath)
    csv_filter = df[df["1000G2015aug_ALL"]<0.5]
    return csv_filter

class GATKVCF:
    def __init__(self, path):
        pass

class Annovar_Filter:
    def __init__(self, path):
        self.path = path
        self.vcf = pd.read_csv(path)
        self.cols = list(self.vcf.columns)
        self.fields = [ "Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
            "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene", "1000g2015aug_all", "1000g2015aug_eas",
            "CLINSIG", "cytoBand"]

    def filter_1KGenome(self, maxRatio=0.05):
        fields = self.fields
        fields_idx = [self.cols.index(field) for field in fields if field in self.cols]
        filter_1 = pd.to_numeric(self.vcf["1000g2015aug_all"], errors='coerce') <= maxRatio
        filter_2 = pd.to_numeric(self.vcf["1000g2015aug_eas"], errors='coerce') <= maxRatio
        return self.vcf[filter_1 & filter_2].iloc[:, fields_idx]

class MergeMultiple:
    def __init__(self, files):
        pass

    def toTable(self):
        pass