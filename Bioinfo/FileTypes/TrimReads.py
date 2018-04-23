
import DataType as DT
from Pipeline.Template import ScriptTemplate
from Bioinfo.FileTypes.Samples import Sample
from Bioinfo.FileTypes.Sequences import SeqStats, QualityStats

class Trimomatic(ScriptTemplate):
    template = """
        java -jar ${path} SE ${infile} ${infile}.trim.fastq SLIDINGWINDOW:4:20 MINLEN:20
        """

    nested = {
    }

    def __init__(self, name = ""):
        super().__init__(name)
        root = DT.OBJECT(name)
        self.path = DT.CHAR("Trimomatic Jar", default="./trimmomatic-0.32.jar",parent=root)
        self.infile = DT.CHAR("In File", default="XX.fastq", parent=root)
        self.root = root