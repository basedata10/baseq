import pandas as pd

def groupfile_check():
    pass

def groupcompare_check():
    pass

class group_compare:
    """
    Groups that should be compared...
    """
    def __init__(self, path):
        self.cols = ["group1", "group2"]
        self.sep = "csv"
        self.path = path

    def read(self):
        return pd.read_table(self.path, sep='\s+')

class sample_group:
    """
    sample and its group title ...
    """
    def __init__(self, path):
        self.cols = ["sample", "group"]
        self.sep = "csv"
        self.path = path

    def check(self):
        pass

    def read(self):
        return pd.read_table(self.path, sep='\s+')