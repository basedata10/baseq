import DataType as DT
from Tool.Utils.Files import File, Folder
from Pipeline.Template import ScriptTemplate
from Bioinfo.FileTypes.Samples import Sample, SampleQC

class SamtoolsTemplate(ScriptTemplate):
    template = """
        samtools view -b -u -S {{insam}} > {{outBam}}.bam
        samtools sort -@ 8 {{outBam}}.bam -o {{outBam}}.sort.bam
        samtools index {{outBam}}.sort.bam
        rm {{insam}} {{outBam}}.bam
        """
    root = DT.OBJECT("Samtools Sort Config")
    samtools = DT.CHAR("Samtools Path", "samtools path", default = "samtools", parent = root)
    inSam = DT.CHAR("Sam File Path", "", parent = root)
    outBam = DT.CHAR("Bam File", "Prefix for bam file", parent = root)

class BWAConfig(ScriptTemplate):
    template = """
        ${bwa} mem -t ${thread} ${ref} ${fq1} ${fq2} > ${samfile}
        """

    def __init__(self, name):
        super().__init__(name)
        root = DT.OBJECT(name)
        self.root = root
        self.bwa = DT.CHAR("BWA Path 路径", default = "bwa", parent=root)
        self.ref = DT.CHAR("Reference", desc = "References For example, /path/to/genome.fa", parent=root)
        self.thread = DT.INT("Thread", default = 8, parent=root)
        self.samfile = DT.CHAR("samtools Path", default = "aligned.sam", parent=root)
        self.fq1, self.fq2 = root.addF("$Fastq1", "$Fastq2")

        folder = Folder("bwa")
        self.sam = File("bwa.sam", parent=folder)
        self.bam = File("bwa.sorted.bam", parent=folder)
        self.flagstat = File("flag.txt", parent=folder)

class BWATemplate(ScriptTemplate):
    template = """
        """
    def __init__(self, name):
        super().__init__(name)
        root = DT.OBJECT(name)
        QC = SampleQC("Sample")
        bwa = BWAConfig("BWA")
        self.nested = {
            "sample" : QC,
            "bwa": bwa
            }
        self.root = root
        self.join(QC.fq1, bwa.fq1)
        self.join(QC.fq2, bwa.fq2)

class Bowtie2(ScriptTemplate):
    template = """
        ${bowtie2} -x ${ref} -1 ${fq1} -2 ${fq2} --no-unal -p ${thread} -t -S ${samfile}
        """
    def __init__(self, name = "Bowtie2 Alignment"):
        super().__init__(name)
        sample = SampleQC("Sample")
        self.nested = {
            "sample": sample
            }
        root = DT.OBJECT(name)
        self.bowtie2 = DT.CHAR("Bowtie2 Path", default = "bowtie2", parent=root)
        self.ref = DT.CHAR("Reference", desc = "References, /path/to/genome.fa", parent=root)
        self.thread = DT.INT("Thread", default = 8, parent=root)
        self.maxLines = DT.INT("Maxlines", default = 5000000, parent=root)
        self.fq1 = DT.CHAR("fq1")
        self.fq2 = DT.CHAR("fq2")
        self.join(sample.fq1, self.fq1)
        self.join(sample.fq2, self.fq2)

        folder = Folder("bowtie2")
        self.sam = File("bowtie2.sam", parent=folder)
        self.bam = File("bowtie2.sorted.bam", parent=folder)
        self.flagstat = File("flag.txt", parent=folder)
        self.root = root


class Hisat2(ScriptTemplate):
    template = """
        ${hisat2} -x ${index} -1 ${fq1} -2 ${fq2} -S ${samfile}
        """
    def __init__(self, name = ""):
        super().__init__(name)
        root = DT.OBJECT(name)
        self.hisat2 = DT.CHAR("Hisat2 Path", default = "hisat2", parent=root)
        self.index = DT.CHAR("Index", desc = "Index for the Hisat2", parent=root)
        self.samfile = DT.CHAR("Sam Files", desc = "OutSamFiles", parent=root)
        self.fq1 = DT.CHAR("Fastq 1", desc = "", parent=root)
        self.fq2 = DT.CHAR("Fastq 2", desc = "", parent=root)
        self.root = root