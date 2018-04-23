import DataType as DT
from Pipeline.Template import ScriptTemplate

class VCFFilter:
    root = DT.OBJECT("Filter")

class Annovar(ScriptTemplate):
    template = """
     perl ${annovar} ${invcf} ${dbpath} -buildver hg19 -out ${outfile} -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_eas,snp138,ljb26_all,cosmic70  -operation g,r,r,f,f,f,f,f,f -nastring .
    """
    def __init__(self, name=""):
        super().__init__(name)
        self.root, self.annovar, self.dbpath, self.buildover, self.invcf, self.outfile = DT.newOBJ(name,
              "Table Annovar#table_annovar.pl",
              "Database Path#/path/to/annovar/humandb",
              "Build Version 版本#hg38", "InFile 输入#Input.vcf", "Outfile 输出路径#Annovar.txt"
              )