import DataType as DT
from Pipeline.Template import ScriptTemplate

class BarcodeSplitter(ScriptTemplate):
    template = """
        import re
        path_1 = "${path_1}"
        path_2 = "${path_2}"
        barcodes = '''${barcodes}'''
        barcodes = ${TextParser}(barcodes, 2)
        files = {}
        counts = {}
        for barcode in barcodes:
            files[barcode[1]] = open("Barcode_"+barcode[0]+".fq", 'w')
            counts[barcode[1]] = 0
        lines = ${filereader}(path_1, 4)
        for line in lines:
            barcode = re.split(":", line[0].strip())[-1].strip()
            if barcode in files:
                files[barcode].writelines(line)
                counts[barcode] += 1
        """

    nested = {
    }

    def __init__(self, name = ""):
        super().__init__(name)
        root = DT.OBJECT("BarcodeSplitter")
        self.path_1 = DT.CHAR("Read 1", parent=root, default="/Users/temp/Documents/chuangye/structure/Tests/Datas/PD10_HYY_M1_index1_ATCACG_L001_R1_001.fastq.gz")
        self.path_2 = DT.CHAR("Read 2", parent=root)
        self.outDir = DT.CHAR("Out Dir", parent=root, default="./")
        self.barcodes = DT.TEXT("Barcodes", desc= "Two columns", parent=root, default="CELL1 ATCACG")
        self.root = root