import DataType as DT
from Bioinfo.FileTypes.Samples import SampleQC
from Bioinfo.Aligner.BWA import BWATemplate
from Bioinfo.SNV.Annovar import Annovar
from Bioinfo.FileTypes.VCF import VCFilter, VCFStats
from Pipeline.Template import ScriptTemplate

class NIPTconfig(ScriptTemplate):
    template = """
        """
    def __init__(self, name = ""):
        super().__init__(name)
        vcfStats = VCFStats()
        vcfFilter = VCFilter()
        self.root = DT.OBJECT(name)
        self.nested = {
            "AlignMother": BWATemplate("Mother"),
            "AlignBaby": BWATemplate("Baby"),
            "annovar" : Annovar("Annovar")
        }