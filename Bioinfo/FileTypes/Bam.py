from Pipeline.Template import ScriptTemplate
import DataType as DT


class SamCounter(ScriptTemplate):

    template = """
        COUNTING THE SAM FILE ......
        return the number of reads for each line in the bed file...
        ${samfile}
        ${bedfile}
        """
    def __init__(self, name="SamCounter"):
        super().__init__(name)
        root = DT.OBJECT(name)
        self.samfile = DT.CHAR("Sam/Bam Path", parent=root)
        self.bedfile = DT.CHAR("Bed Path", parent=root)
        self.root = root
        self.nested = {
            }

class ShowRangeBaseContent(ScriptTemplate):
    template = """
        samtools mpileup ${bam} -c ${chr}:${start}-${end}
        """
    pass