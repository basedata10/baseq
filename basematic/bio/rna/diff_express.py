from subprocess import call
from basematic.setting import r_script_dir
import os, sys, json
from basematic.mgt import get_config, run_cmd

def deseq2(tpmfile, countfile, groupfile, comparefile, outpath):
    """ Run DNACopy.R file ...
    input:
        tmp file
        count file
        group file: tell the group name for each
            samplename/groups/
        compare file: which groups should be compared...
            compare_name/group1/group2
        output path
    output:
        under the output path, for each

    """
    Rscript = get_config("RNA", "deseq")
    script = os.path.join(r_script_dir, "DESeq2.R")
    cmd = "{} {} {} {} {} {} {}".format(Rscript, script, tpmfile, countfile, groupfile, comparefile, outpath)
    if not os.path.exists(outpath):
        os.mkdir(outpath)
        print("[info] Create OutDir {}".format(outpath))
    try:
        run_cmd("Haha...", cmd)
    except:
        sys.exit("[error] Failed to run the Normalize Rscript ...")


def pack_DeSeq2_result(groupfile, comparefile, outpath):
    import pandas as pd
    with open(comparefile, "r") as file:
        lines = file.readlines()
    groups = ["_".join(line.split()) for line in lines][1:]

    import xlsxwriter
    workbook = xlsxwriter.Workbook(os.path.join(outpath, 'DESeq.xlsx'))
    workbook.formats[0].set_font_size(12)
    workbook.formats[0].set_font_name('arial')
    format_main = workbook.add_format({'bold': False, 'font_size': 12, 'font_name': 'arial'})
    format_header = workbook.add_format({'bold': True, 'font_size': 15, 'font_name': 'arial'})

    for group in groups:

        path = os.path.join(outpath, "genes_up.{}.txt".format(group))
        print(group, path)
        table = pd.read_table(path, sep=" ")
        data = json.loads(table.to_json(orient="split"))
        colnames = data["columns"]
        qcpage = workbook.add_worksheet(group)
        qcpage.merge_range('A1:E1', 'Up Regulated', format_header)
        qcpage.write('A2', 'Gene', format_main)
        qcpage.write('B2', 'Log2FC', format_main)
        qcpage.write('C2', 'Pvalue', format_main)
        qcpage.write('D2', colnames[2], format_main)
        qcpage.write('E2', colnames[3], format_main)
        for idx_row, gene in enumerate(data['index']):
            qcpage.write(idx_row + 2, 0, gene, format_main)
            for idx_col in range(4):
                qcpage.write(idx_row + 2, idx_col+1, data['data'][idx_row][idx_col], format_main)

        path = os.path.join(outpath, "genes_down.{}.txt".format(group))
        print(path)
        table = pd.read_table(path, sep=" ")
        data = json.loads(table.to_json(orient="split"))
        colnames = data["columns"]
        qcpage.merge_range('G1:K1', 'Down Regulated', format_header)
        qcpage.write('G2', 'Gene', format_main)
        qcpage.write('H2', 'Log2FC', format_main)
        qcpage.write('I2', 'Pvalue', format_main)
        qcpage.write('J2', colnames[2], format_main)
        qcpage.write('K2', colnames[3], format_main)
        for idx_row, gene in enumerate(data['index']):
            qcpage.write(idx_row + 2, 6, gene, format_main)
            for idx_col in range(4):
                qcpage.write(idx_row + 2, idx_col+7, data['data'][idx_row][idx_col], format_main)

    #build an Excel...
    workbook.close()