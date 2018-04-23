import DataType as DT
from Pipeline.Template import ScriptTemplate
from Bioinfo.FileTypes.Samples import Sample, SampleQC

class AlignerMaster(ScriptTemplate):
    template = """
        {% if ${method} == "BWA" %}
            BWA 使用中哈哈哈哈
        {% endif %}

        {% if ${method} == "BOWTIE2" %}
            BOWTIE2 使用中哈哈哈哈
        {% endif %}
        """
    def __init__(self, name):
        super().__init__(name)
        root = DT.OBJECT(name)
        self.method = DT.Selection("Software", options=["BWA", "BOWTIE2"])
        self.root = root