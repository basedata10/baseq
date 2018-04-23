from Bioinfo.FileTypes.Samples import SampleQC
from Pipeline.Template import ScriptTemplate
import DataType as DT

class Salmon(ScriptTemplate):
    template = """
        ${salmon} quant -i ${index} -l A -1 ${fq1} -2 ${fq2} -p ${thread} -o ${outDir}
        """

    def __init__(self, name = "salmon"):
        super().__init__(name)
        root = DT.OBJECT(name)
        self.salmon = DT.CHAR("Salmon Path", default="salmon", parent=root)
        self.index = DT.CHAR("Salmon Index 索引", desc="Salmon index folder", parent=root)
        self.gene_map = DT.CHAR("Gene Map 基因映射", desc="Transcript to gene names", parent=root)
        self.thread = DT.INT("Thread", parent=root, default=8)
        self.outDir = DT.CHAR("OutDir", default="./", parent=root)
        self.root = root