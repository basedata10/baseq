import DataType as DT
from Pipeline.Template import ScriptTemplate

class VCFStats:
    root = DT.OBJECT("VCF Stats")
    root, *stats = DT.newOBJ("VCF Stats", "Total", "Allel frequency")

class VCFilter:
    root, *other = DT.newOBJ("VCF Filter","Chr","Start","End","Genes")

class SimpleVCF(ScriptTemplate):
    template = """
        XXXX
        """
    root = DT.OBJECT("SimpleVCF")