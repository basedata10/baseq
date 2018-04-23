import DataType as DT
from Pipeline.Template import ScriptTemplate
from Bioinfo.FileTypes.Samples import SampleQC
from Bioinfo.FileTypes.Sequences import SeqStats, QualityStats

class QCConfig(ScriptTemplate):
    template = """
        fq1 = ${readfile}("${fq1}", 4)
        seq1 = []
        qua1 = []
        for idx, lines in enumerate(fq1):
            if idx<=100000:
                seq1.append(lines[1].strip())
                qua1.append(lines[3].strip())
        R1_base = ${seqstats}(seq1)
        R1_qual = ${quality}(qua1)
        """

    def __init__(self, name = ""):
        super().__init__(name)
        self.root = DT.OBJECT(name)
        sample = SampleQC()
        self.fq1 = sample.fq1
        self.fq2 = sample.fq2

        self.nested = {
            "QC" : sample,
            "seqstats" : SeqStats(),
            "quality" : QualityStats()
        }

        self.export("R1_base", QualityStats().result)
        self.export("R1_qual", SeqStats().result)

    def viewResults(self):
        return None