from sqlalchemy import *
from sqlalchemy.orm import *
from sqlalchemy.sql.expression import ColumnClause
from sqlalchemy.sql import table, column, select, update, insert

from sqlalchemy import create_engine
engine = create_engine('mysql+mysqldb://geneuser:geneuser002@192.168.2.1:3306/genes')
metadata = MetaData(engine)

geneinfo = Table('rifs', metadata,
    Column('id', Integer, autoincrement=True, primary_key = True),
    Column('taxid', Integer),
    Column('geneid', Integer, index=True),
    Column('pubmedid', Integer),
    Column('text', String(1000)))

metadata.create_all(engine)
conn = engine.connect()

index = 0
bulk = []

with open("../generifs_basic", 'r') as infile:
    for line in infile:
        if line[0]=="#":
            continue
        tags = line.strip().split("\t")

        u = dict(taxid=int(tags[0]),
                 geneid=int(tags[1]),
                 pubmedid=int(tags[2].split(",")[0]),
                 text=tags[4])

        bulk.append(u)
        if len(bulk)==100:
            conn.execute(geneinfo.insert(), bulk)
            print("SUCCESS", str(index))
            bulk = []
        index += 1

if len(bulk)>0:
    conn.execute(geneinfo.insert(), bulk)
    print("SUCCESS", str(index))