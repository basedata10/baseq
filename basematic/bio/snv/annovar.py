from basematic.mgt.check import isExist
import os

class Annovar:
    def __init__(self, script, db, version="hg38"):
        self.script = script
        self.db = db
        self.version = version

    def help_download(self):
        cmd = """
            The download command could be:
            annotate_variation.pl -buildver hg38 -downdb 1000g2015aug ./db_1000/
            annotate_variation.pl -buildver hg38 -downdb 1000g2015aug ./db_1000/
            """
        print(cmd)

    def check_script(self):
        isExist(os.path.join(self.script, "annotate_variation.pl"), "annotate_variation.pl")
        isExist(os.path.join(self.script, "convert2annovar.pl"), "convert2annovar.pl")
        isExist(os.path.join(self.script, "table_annovar.pl"), "table_annovar.pl")
        print("[Success] The Annovar Scripts Exists ...")

    def check_database(self):
        dbs = ['refGene','genomicSuperDups','esp6500siv2_all',
               #'1000G2015aug_ALL', '1000G2015aug_EAS',
               'exac03', 'avsnp147', 'dbnsfp30a', 'clinvar_20170130', 'cosmic70', 'dbscsnv11', 'cytoBand']
        for db in dbs:
            isExist(os.path.join(self.db, "{}_{}.txt".format(self.version, db)), "Annovar DataBase {}".format(db))
        print("[Success] All the Annovar Database is valid")

    def generate_annovar_script(self):
        self.check_script()
        self.check_database()