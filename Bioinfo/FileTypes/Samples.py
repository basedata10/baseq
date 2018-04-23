import DataType as DT
from Pipeline.Template import ScriptTemplate

class Sample:

    def __init__(self, name="Sample Infos"):
        self.root = DT.OBJECT("Sample")
        self.sample = DT.CHAR("Sample Name", parent=self.root, default="Sample")
        self.fastq1 = DT.CHAR("Fastq 1", parent=self.root)
        self.fastq2 = DT.CHAR("Fastq 2", desc="For Pair End Sequences", parent=self.root)
        self.root.name = name
        self.sample.name = name

class SampleQC(ScriptTemplate):

    template = """
        FASTQC ${fq1} ${fq2}
        CHECK IF THE FASTQ FILE EXISTS!
        """
    def __init__(self, name="Sample Infos"):
        super().__init__(name)
        root = DT.OBJECT(name)
        self.sample = DT.CHAR("Sample Name", parent=root, default="Sample")
        self.fq1 = DT.CHAR("Fastq 1", parent = root)
        self.fq2 = DT.CHAR("Fastq 2", desc="For Pair End Sequences", parent=root)
        self.root = root