import DataType as DT
from Pipeline.Template import ScriptTemplate

intro = """
Cufflinks assembles transcripts, estimates their abundances, and tests for differential expression and
regulation in RNA-Seq samples. It accepts aligned RNA-Seq reads and assembles the alignments into a
parsimonious set of transcripts. Cufflinks then estimates the relative abundances of these transcripts based on
how many reads support each one, taking into account biases in library preparation protocols.
"""

link = """http://cole-trapnell-lab.github.io/cufflinks/"""

class Cufflinks(ScriptTemplate):
    template = """
        ${cuff} -o ${Dir} -p ${thread} -G ${GTF} ${bamfile}
        """

    def __init__(self, name=""):
        super().__init__(name)
        root = DT.OBJECT(name)
        self.cuff = DT.CHAR("cufflinks Path", default="cufflinks", parent=root)
        self.Dir = DT.CHAR("OutDir", default="./cufflinks", parent=root)
        self.GTF = DT.CHAR("Reference", desc="References, /path/to/gencode.gtf", parent=root)
        self.thread = DT.INT("Thread", default=8, parent=root)
        self.bamfile = DT.CHAR("Sorted Bam", default="aligned.sort.bam", parent=root)
        self.root = root