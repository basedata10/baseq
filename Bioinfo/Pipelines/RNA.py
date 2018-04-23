import DataType as DT
from Pipeline.Template import ScriptTemplate

class DEGenes:
    name = DT.CHAR("Gene")
    pvalue = DT.FLOAT("pvalue")
    foldchange = DT.FLOAT("FoldChange")
    G1_mean = DT.FLOAT("Group1 Mean")
    G2_mean = DT.FLOAT("Group1 Mean")

from Bioinfo.FileTypes.Samples import SampleQC
from Bioinfo.Aligner.BWA import Hisat2
from Bioinfo.Aligner.cufflinks import Cufflinks

class RNAConfig(ScriptTemplate):
    template = """
        """

    def __init__(self, name = ""):
        super().__init__(name)
        self.root = DT.OBJECT("RNA-Seq")

        QC = SampleQC("Samples")
        hisat = Hisat2("Hisat2 快速基因组比对")
        cuff = Cufflinks("Cufflinks 估计RPKM")

        self.nested = {
            "QC" : QC,
            "hisat" : hisat,
            "cuff" : cuff
        }

        self.join(hisat.samfile, cuff.bamfile)
        self.join(QC.fq1, hisat.fq1)
        self.join(QC.fq2, hisat.fq2)