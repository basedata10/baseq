#Define the datatypes in the result...
from DataType import *
QC_base, *bases = newDF("Base Content", "A", "T", "C", "G")
QC_score, *quals = newDF("Base Score", "Score")
Align_stats, *stats = newOBJ("Alignment", "Total", "Mapped", "Ratio")
CNV_Figure, *figures = newOBJ("CNV Plots", "200K Bin Size", "1M Bin Size")
CNV_Table, *CNV = newDF("CNV List", "Chr", "start", "end", "Copy Number", "Cytoband")

#Define the template
from HTML import Boxes
template = Boxes("double", rowcount=1, W=900)
template.nodes = [
    QC_base.viewLineChart("", bases),
    QC_score.viewBarChart("", quals),
    Align_stats.viewLevel(),
    CNV_Table.viewTable()
]
print(template.toPage())

