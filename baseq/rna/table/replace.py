
def ensg_to_genename(namefile, table, outteble):
    import pandas as pd
    df_name = pd.read_table(namefile, header=0, names=["ENSG", "GENE"])
    ensg_gene = {}
    for index, row in df_name.iterrows():
        ensg_gene[row['ENSG']] = row['GENE']
    print(ensg_gene)

    df_gene = pd.read_csv(table, index_col=0)
    ganenames = []
    for index, row in df_gene.iterrows():
        if index in ensg_gene:
            ganenames.append(ensg_gene[index])
        else:
            ganenames.append(index)
    df_gene.insert(0, 'gene', ganenames)
    df_gene.to_csv(outteble, index=False)

def mean_gene_expression(table, outtable):
    pass