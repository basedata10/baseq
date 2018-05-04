
script = """
more {} | awk '$3=="gene"' | awk '{print $10, $16}' | perl -ne '@a=split(/\.|"/,$_);print "$a[1]\t$a[4]\n"' >{}
"""

def gencode_get_ENSG_genename(path, outfile):
    pass