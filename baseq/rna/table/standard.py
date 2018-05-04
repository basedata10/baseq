"""
Input:
    Single Table Path
    Multiple Table Path
    Groups Information
    Group Compare Table

A standard EXPTABLE file should be:
    csv, comma seprated;
    first column is gene;
    other columns is sample;

Some operations of a standard EXPTable
    Write Simplify Table
    Write HDF5

Table Vs Table
    Compare between the EXPTables

"""

class EXPTable:
    def __init__(self, path=[], groups="", group_comp=""):
        pass

    def read_exp_table(self):
        pass

    def mean(self):
        pass

    def corrs(self):
        pass