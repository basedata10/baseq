import DataType as DT
from DataType import *
from Pipeline.Template import ScriptTemplate

def Stats():
    seqs = "{}"
    seqcounts = max(len(seqs), 1)
    seqlen = len(seqs[0])
    content = {}
    for base in ['A', 'T', 'C', 'G']:
        content[base] = [0] * seqlen
    for seq in seqs:
        for idx, base in enumerate(seq):
            if base in ['A', 'T', 'C', 'G'] and idx < seqlen:
                content[base][idx] += 1
    for base in ['A', 'T', 'C', 'G']:
        content[base] = [float(x) / seqcounts for x in content[base]]
    return content

class SeqStats(ScriptTemplate):
    script = """
        def ${name}(seqs):
            seqcounts = max(len(seqs), 1)
            seqlen = len(seqs[0])
            content = {}
            for base in ['A', 'T', 'C', 'G']:
                content[base] = [0]*seqlen
            for seq in seqs:
                for idx, base in enumerate(seq):
                    if base in ['A', 'T', 'C', 'G'] and idx<seqlen:
                        content[base][idx] += 1
            for base in ['A', 'T', 'C', 'G']:
                content[base] = [float(x)/seqcounts for x in content[base]]
            return content
        """
    def __init__(self, name="seqstats"):
        super().__init__(name)
        self.subTemplate(self.script, name=self.ID)
        self.root, = newOBJ("Sequence Base Stats")
        self.result, A, T, C, G = newDF("Base Content", "%3A", "%3T", "%3C", "%3G")

class QualityStats(ScriptTemplate):
    script = """
        def ${name}(seqs):
            seqlen = len(seqs[0])
            quality = [0] * seqlen
            for seq in seqs:
                for idx, base in enumerate(seq):
                    if idx < seqlen:
                        quality[idx] += ord(base) - 33
            return [int(q/len(seqs)) for q in quality]
        """

    def __init__(self, name="qualStats"):
        super().__init__(name)
        self.subTemplate(self.script, name=self.ID)
        self.root, = newOBJ("Seq Quality")
        self.result = DT.INT("Quality", isArray=True)

class MostFrequent(ScriptTemplate):
    pass