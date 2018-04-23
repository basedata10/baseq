import DataType as DT

#The Enrichment Quality
config_enrichQC = DT.OBJECT("Enrichment QC")
DT.CHAR("Bam File", "Should be a sorted Bam file.", parent=config_enrichQC)
DT.CHAR("Bed File", "The exome or panel bed file.", parent=config_enrichQC)

#Run manage Log File?
run_datas = DT.OBJECT("Run Enrichment")

#Result Maybe File?
monitor_enrich = DT.OBJECT("富集质量展示")
DT.FLOAT("Total总数", parent = monitor_enrich)
DT.FLOAT("Aligned比对率", parent = monitor_enrich)
DT.FLOAT("Padding比率", parent = monitor_enrich)
DT.FLOAT("Region比率", parent = monitor_enrich)
DT.FLOAT("Region平均深度", parent = monitor_enrich)
DT.FLOAT("Region深度", isArray=True, parent = monitor_enrich)
DT.FLOAT("Padding深度", isArray=True, parent = monitor_enrich)

#Monitor Stratagies
statagies = DT.OBJECT("展示方式")

#Run Logic For Config
