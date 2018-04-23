from string import Template
from Tool.Utils.Files import File, Folder
from Bioinfo.Aligner.Script.Aligner import *

from Pipeline import ShellScript

bwa = """
    ${bwa} mem -t ${thread} ${ref} ${fq1} ${fq2} > ${samfile}
    """

class BWA(ShellScript):
    def __init__(self, bwa, thread, fq1, fq2, ref, samfile):
        pass