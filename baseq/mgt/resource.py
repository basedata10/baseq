import os
import pandas as pd
from Bioinfo.DataBase import DataBase
from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

Base = declarative_base()

class LocalResource:

    def __init__(self):
        engine = DataBase(os.path.join(os.path.expanduser("~"), "baseq.sqlite3")).db
        Base.metadata.create_all(engine)
        self.DBSession = sessionmaker(bind=engine)

class Resource(Base):
    __tablename__ = 'resource'
    name = Column(String(), primary_key=True)
    path = Column(String())
    infos = Column(String())

class ResourceModel(LocalResource):
    def __init__(self):
        super().__init__()

    def add(self, name="", path="", infos=""):
        session = self.DBSession()
        new_user = Resource(name=name, path=path, infos=infos)
        session.merge(new_user)
        session.commit()
        session.close()

    def add_lists(self, name="", path="", infos=""):
        session = self.DBSession()
        new_user = Resource(name=name, path=path, infos=infos)
        session.merge(new_user)
        session.commit()
        session.close()

    def list(self):
        session = self.DBSession()
        for row in session.query(Resource).all():
            print(row.id, row.name, row.path, row.infos)
        session.close()

    def update(self):
        pass

    def delete(self):
        pass

    def list_name(self, names=[]):
        session = self.DBSession()
        result = []
        for name in names:
            count = session.query(Resource).filter(Resource.name == name).count()
            if count == 0:
                self.add(name = name)
        for row in session.query(Resource).filter(Resource.name.in_(names)).all():
            result.append([row.name, row.path])
        session.close()
        return pd.DataFrame(result, columns=["name", "path"])

if __name__ == "__main__":
    soft = ResourceModel()
    soft.add("bwa", "/path/to/bwa", "")
    print(soft.list_name(['bwa']))