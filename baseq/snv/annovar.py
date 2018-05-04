from baseq.mgt.check import isExist
import os
from baseq.mgt import get_config, run_cmd

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

annovar_cmd_script = """
{0}/convert2annovar.pl -format vcf4 {1}_raw_snps.vcf >{1}_snps.avinput
{0}/table_annovar.pl {1}_snps.avinput {2} -buildver hg38 -out {1} -remove -protocol \
refGene,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,avsnp147,dbnsfp30a,clinvar_20170130,cosmic70,dbscsnv11,cytoBand \
-operation g,f,f,f,f,f,f,f,f,f,r -nastring . -csvout"""
def run_annovar(sample,genome,run=True):
    annovar = get_config("Annovar","annovar")
    annovar_ref = get_config("Annovar","annovar_db_hg38")
    annovar_cmd = annovar_cmd_script.format(annovar, sample, annovar_ref)
    if run:
        run_cmd("convert vcf file to aninput format","".join(annovar_cmd))
    return annovar_cmd

