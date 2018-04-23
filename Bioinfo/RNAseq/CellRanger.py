import DataType as DT
from Pipeline.Template import ScriptTemplate
from Bioinfo.FileTypes.Samples import Sample, SampleQC

class CellRanger(ScriptTemplate):
    template = """
        mkdir ${sampleDir}
        cd ${sampleDir}
        ln -s ${fq1} ${sampleDir}/Sample1.1.fq.gz
        ln -s ${fq2} ${sampleDir}/Sample1.2.fq.gz
        ${cellranger} count --id=${sample} \
            --transcriptome=${reference} \
            --fastqs=./sampleFolder 
        """

    def __init__(self, name):
        super().__init__(name)
        root = DT.OBJECT(name)
        self.sample = DT.CHAR("Sample Name", "", default = "10X", parent = root)
        self.cellranger = DT.CHAR("Cellrange Path", "", parent = root,
                                 default = "/mnt/gpfs/Users/wufan/soft_inDrop/cellranger-2.0.1/cellranger")
        self.reference = DT.CHAR("Reference", "",
                                 default = "/mnt/gpfs/Users/wufan/p07_10X/reference/refdata-cellranger-GRCh38-1.2.0",
                                 parent = root)
        self.sampleDir = DT.CHAR("Path", "", default = "./", parent = root)
        self.root = root
        self.nested = {
           "sample" : SampleQC("Sample")
            }

class Versatil(ScriptTemplate):
    template = """
        #Extract the barcode
        cd ${outdir}
        python 00_get_cellbarcode.py ./_config.json ${fq1} CB.count.${sample}.pickle ${method} ${lowBC}
        #Stats barcode
        python 01_cellbarcode_stats.py _config.json CB.count.${sample}.pickle ./CB.stats.${sample}.json ${method} ${maxCell} ${minReads}
        #Split barcode
        python 02_split_barcode.py ${fq1} ${fq2} ./CB.stats.${sample}.json ./process_${sample} ${method}
        #Star alignment...
        ${star_path} --genomeLoad LoadAndRemove --genomeDir ${star_ref} --readFilesIn $1 --runThreadN 30 --outSAMunmapped Within --outSAMtype BAM Unsorted --outFileNamePrefix $2
        samtools sort -n -@ 30 $2Aligned.out.bam -o $2sort.bam
        #Tags
        python 04_Attach_Transcript.py STAR.sorted.bam ./UMIs.Barcode.txt Barcode
        #Levels...
        python 05_levels.py reads_indrop2 0.1 reads_indrop2/Counts_indrop2_1.txt reads_indrop2/Reads_indrop2_1.txt reads_indrop2/StatsUMI_indrop2_1.txt
        """
    def __init__(self, name):
        super().__init__(name)
        self.root = DT.OBJECT(name)
        self.method, self.ref, self.outdir = self.root.addF(["Method", "10X", "inDrop", "Drop-Seq"], "Reference#/Path/to/reference", "Out Path#./")
        self.sampleinfo, self.sample, self.fq1, self.fq2 = DT.newOBJ("Sample Information", "Name#sample", "Read1 Path#sample.1.fq.gz", "Read2 Path#sample.2.fq.gz")
        self.cell_barcode, self.maxCell, self.minReads, self.lowBC = DT.newOBJ("Cell Barcode Filtering", "%Max cell number#20000", "%Min Reads#1000", "%Lower Threshold#50")
        self.STAR, self.star_path, self.star_ref = DT.newOBJ("STAR Alignment", "STAR Path#STAR", "STAR Reference Index#/path/to/reference")
        self.root.add_node(self.sampleinfo)
        self.root.add_node(self.cell_barcode)
        self.root.add_node(self.STAR)