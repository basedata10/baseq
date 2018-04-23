import DataType as DT
from Pipeline.Template import ScriptTemplate

class BamStats(ScriptTemplate):

    template = """
        samtools view ${path} ${chr}:${start}-${end}
        """
    def __init__(self, name = ""):
        super().__init__(name)
        config = DT.OBJECT("BamStats")
        self.path = DT.CHAR("Bam Path", "Should be a sorted bam containing an index", parent=config)
        self.chr = DT.CHAR("Chromosome", default="chr1", parent=config)
        self.start = DT.CHAR("Start", parent=config)
        self.end = DT.CHAR("End", parent=config)
        self.root = config
        self.out = DT.OBJECT("Bam Stats Output")

DT.OBJECT("Bam Depth")