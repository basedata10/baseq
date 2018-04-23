import sqlite3
import os
import sqlalchemy
from sqlalchemy import *
from basematic.file.Folder import EnsurePath
from sqlalchemy.orm import sessionmaker

class DBcols:
    def __init__(self):
        self.cols = []

    def ID(self, name):
        self.cols.append(Column(name, sqlalchemy.Integer, primary_key=True))

    def String(self, name, maxLen=-1):
        if maxLen:
            self.cols.append(Column(name, sqlalchemy.String(maxLen)))
        else:
            self.cols.append(Column(name, sqlalchemy.Integer))

    def Int(self, name):
        self.cols.append(Column(name, sqlalchemy.Integer))

class TableModel:
    def __init__(self, table):
        self.table = table

    def add(self, **kargs):
        self.table.insert().execute(**kargs)

    def list(self):
        return list(self.table.select().execute())

class DataBase:
    def __init__(self, path="test.sqlite3", db="sqlite"):
        self.path = path
        self.db = create_engine('sqlite:///' + self.path)
        self.db.echo = False
        self.metadata = MetaData(self.db)
        # create a configured "Session" class create a Session
        Session = sessionmaker(bind=self.db)
        self.session = Session()

    def table(self, name, columns):
        table = Table(name, self.metadata, *columns)
        try:
            table.create()
            print("Table Created.")
        except:
            print("Table Already Exists")
        return table

    def test(self):
        users = Table('users', self.metadata,
                      Column('user_id', Integer, primary_key=True),
                      Column('name', String(40)),
                      Column('age', Integer),
                      Column('password', String),
                  )
        try:
            users.create()
        except:
            pass

        i = users.insert()
        i.execute(name='Mary', age=30, password='secret')
        i.execute({'name': 'John', 'age': 42},
                  {'name': 'Susan', 'age': 57},
                  {'name': 'Carl', 'age': 33})

        i.execute(name="some name")

        rs = users.select().execute()
        for row in rs:
            print(row.name, 'is', row.age, 'years old')

if __name__ == "__main__":
    DB = DataBase(path = "test.sqlite3")
    cols = [
            Column('user_id', Integer, primary_key=True),
            Column('name', String(40)),
            Column('age', Integer),
            Column('password', String)]

    user = DB.table("user", cols)
    #Use table methods
    user.insert().execute(name='Mary 123', age=30, password='secret')

    #Use TableModel
    user_model = TableModel(user)
    user_model.add(name='James', age=33)
    print(user_model.list())