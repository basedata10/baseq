# Softs...
import os
from sqlalchemy import *
from basematic.file.DataBase import DataBase

class Jobs:
    def __init__(self):
        dir = os.path.expanduser("~")
        self.db = DataBase(os.path.join(dir, "basematic.sqlite3"))
        self.table = self.db.table("jobs", [
            Column('id', Integer, primary_key=True),
            Column('name', String),
            Column('path', Integer),
            Column('start', Integer),
            Column('status', Integer),
            Column('infos', String)
        ])

    def add(self, name="", path="", infos="", status=""):
        self.table.insert().execute(name = name, path=path, status=status, infos=infos)

    def list(self):
        for row in self.table.select().execute():
            print(row)

    def query(self):
        pass

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
Base = declarative_base()

if __name__ == "__main__":

    jobs = Jobs()
    jobs.add("haha", "haha", "haha", "hahe")
    jobs.list()