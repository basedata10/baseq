# import re
# from sqlalchemy import *
# from sqlalchemy.orm import *
# from sqlalchemy.sql.expression import ColumnClause
# from sqlalchemy.sql import table, column, select, update, insert
#
# from sqlalchemy import create_engine
# engine = create_engine('mysql+mysqldb://geneuser:geneuser002@192.168.2.1:3306/genes')
# metadata = MetaData(engine)
#
# geneinfo = Table('info_hg19', metadata,
#     Column('geneid', Integer, primary_key = True),
#     Column('symbol', String(40), unique=True, index=True),
#     Column('synonyms', String(200)),
#     Column('dbxrefs', String(300)),
#     Column('chrom', String(20), index=True),
#     Column('cytoband', String(80), index=True),
#     Column('desc', String(200)),
#     Column('type', String(100)),
#     Column('start', Integer, index=True),
#     Column('end', Integer, index=True),
#     Column('abs_start', BigInteger, index=True),
#     Column('abs_end', BigInteger, index=True),
#     )
#
# metadata.create_all(engine)
# conn = engine.connect()
#
# gene_pos = {}
# with open("../gencode.v19.pos", 'r') as infile:
#    for line in infile:
#        tags = line.strip().split("\t")
#        gene_pos[tags[3]] = tags[0:3]
#
# chr_abspos = {}
# with open("../genome_chr_absolute_hg19.txt") as infile:
#     for line in infile:
#         tags = re.split(r'\s+', line.strip())
#         print tags
#         chr_abspos[tags[1][3:]] = int(tags[2])
#
# print chr_abspos
#
# with open("../Homo_sapiens.gene_info", 'r') as infile:
#     for line in infile:
#         if line[0]=="#":
#             continue
#         tags = line.strip().split("\t")
#         u = dict(geneid=int(tags[1]),
#                  symbol=tags[2],
#                  synonyms=tags[4],
#                  dbxrefs=tags[5],
#                  chrom=tags[6],
#                  cytoband=tags[7],
#                  desc=tags[8],
#                  type=tags[9])
#         if tags[2] in gene_pos:
#             pos = gene_pos[tags[2]]
#             if pos[0][3:] == tags[6]:
#                 print "SUCCESS", pos
#                 u["start"] = int(gene_pos[tags[2]][1])
#                 u["end"]   = int(gene_pos[tags[2]][2])
#                 u["abs_start"] = int(gene_pos[tags[2]][1])+chr_abspos[tags[6]]
#                 u["abs_end"] = int(gene_pos[tags[2]][2])+chr_abspos[tags[6]]
#                 i = geneinfo.insert()
#                 try:
#                     r1 = conn.execute(i, **u)
#                 except:
#                     print "DUPLICATIONSSS!"
#             else:
#                 print "FAILED", pos, tags[6]
#         else:
#             print "DO NOT EXISTS", tags[2]