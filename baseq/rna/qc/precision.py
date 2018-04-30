# (c) 2018 <friedpine@gmail.com>
# This file is part of Baseq

def group_similarities(groups, qcinfos, exp_table):
    '''
    report the similarities between the samples in the group...
    input:
        TPM table
        Group Design
    return:
        FC2genes: the total number of FC>2 genes in a sample compared to the average;
        MeanCorrs: the mean pairwise correlation
    '''
    print("[info] Group similarity")