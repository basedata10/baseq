import sys, os, time, re, json
from subprocess import Popen, PIPE, call
from baseq.mgt import get_config, run_cmd
from baseq.fastq.sample_file import check_sample_files
import pandas as pd

def run_build_hisat_index():
    pass

script = """
cd {0}
{1} -x {2} -1 {3} -2 {4} -p 8 -S hisat2_align.sam
{5} view -b -u -S hisat2_align.sam > hisat2_align.bam
{5} sort -@ 8 hisat2_align.bam -o hisat2_sorted.bam
{5} index hisat2_sorted.bam
{5} flagstat hisat2_sorted.bam > flatstat.txt
rm hisat2_align.sam hisat2_align.bam
{6} -q -o ./ -p 8 -G {7} hisat2_sorted.bam
"""
script1 = """
cd {0}
{1} -x {2} -U {3} -p 8 -S hisat2_align.sam
{4} view -b -u -S hisat2_align.sam > hisat2_align.bam
{4} sort -@ 8 hisat2_align.bam -o hisat2_sorted.bam
{4} index hisat2_sorted.bam
{4} flagstat hisat2_sorted.bam > flatstat.txt
rm hisat2_align.sam hisat2_align.bam
{5} -q -o ./ -p 8 -G {6} hisat2_sorted.bam
"""

def run_hisat(fq1, fq2, genome, outdir, run=True):
    hisat = get_config("RNA", "hisat")
    samtools = get_config("RNA", "samtools")
    cufflinks = get_config("RNA", "cufflinks")
    hisat_ref = get_config("RNA_ref_"+genome, "hisat_index")
    cufflinks_anno = get_config("RNA_ref_"+genome, "cufflinks_anno")

    # Run hisat, samtools and cufflinks
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("[info] Create outdir in: {}".format(outdir))
    if fq1 and fq2 and os.path.exists(fq1) and os.path.exists(fq2):
        hisat_pipe_cmd = script.format(outdir, hisat, hisat_ref, fq1, fq2, samtools, cufflinks, cufflinks_anno)
    elif fq1 and os.path.exists(fq1):
        hisat_pipe_cmd = script1.format(outdir, hisat, hisat_ref, fq1, samtools, cufflinks,cufflinks_anno)
    else:
        pass

    if run:
        run_cmd("hisat and cufflinks analysis","".join(hisat_pipe_cmd))
    return hisat_pipe_cmd

def run_multiple_hisat(path, genome, outdir):
    samples = check_sample_files(path)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("[info] Create outdir in: {}".format(outdir))
    script_lists = []
    script_main_path = os.path.join(outdir, "work.sh")
    qsub_main_path = os.path.join(outdir, "qsub_work.sh")
    for sample in samples:
        name = sample[0]
        path = os.path.join(outdir, name)
        script_path = os.path.join(path, "work.sh")
        if not os.path.exists(path):
            os.mkdir(path)
            print("[info] Create outdir for sample {} in: {}".format(name, path))
        else:
            print("[info] Outdir for sample {} exists in: {}".format(name, path))
        if sample[2]:
            script_cmd = run_hisat(sample[1], sample[2], genome, path, False)
        else:
            script_cmd = run_hisat(sample[1], "", genome, path, False)
        with open(script_path, "w") as file:
            file.writelines("#!/bin/bash"+"\n"+script_cmd+"\n")
            print("[info] work script written in {}".format(script_path))
        script_lists.append(script_path)

    #write the main script
    with open(script_main_path, "w") as file:
        file.writelines("\n".join(["bash " + x for x in script_lists]))
    with open(qsub_main_path, "w") as file:
        file.writelines("\n".join(["qsub -cwd -l vf=8g,p=8 " + x for x in script_lists]))
    print("[info] Main script written in {}".format(script_main_path))

# Build FPKM and QC ...
def build_FPKM_table(processdir, samplefile, outpath):
    samples = check_sample_files(samplefile)
    sample_names = [sample[0] for sample in samples]
    tpm_file_path = outpath + "fpkm.txt"
    qc_file_path = outpath + "qc.txt"
    tpm = {}
    qc = ["\t".join(['sample', 'reads', 'mapped', 'ratio', 'genecounts'])]
    for sample in sample_names:
        #build TPM table
        sample_TPM = []
        salmon_gene_path = os.path.join(processdir, sample, 'genes.fpkm_tracking')
        if not os.path.exists(salmon_gene_path):
            print("[info] File not exists for {}".format(sample))
            continue
        with open(salmon_gene_path, 'r') as infile:
            infile.readline()
            for line in infile:
                infos = re.split("\t", line)
                gene = infos[4]
                print(gene, float(infos[9]))
                sample_TPM.append(float(infos[9]))
                if not gene in tpm:
                    tpm[gene] = [float(infos[9])]
                else:
                    tpm[gene].append(float(infos[9]))

        #Genes Detected
        genes_TPM_1 = sum([1 for x in sample_TPM if x>=1])

        #build QC data
        metainfo = os.path.join(processdir, sample, 'flatstat.txt')
        with open(metainfo, "r") as file:
            qc_info = file.readlines()
            datas = [x.split(" ") for x in qc_info]
            qc_sample = [sample, int(datas[0][0]), int(datas[4][0]), float(int(datas[4][0]))/int(datas[0][0]), genes_TPM_1]
            qc.append("\t".join([str(x) for x in qc_sample]))

    #Write TPM
    tpm_file = open(tpm_file_path, "w")
    tpm_file.write("\t".join(["gene"] + sample_names) + "\n")
    for gene in tpm:
        sum_tpm_genes = sum(tpm[gene])
        if sum_tpm_genes > 0:
            tpm_file.write("\t".join([gene] + [str(x) for x in tpm[gene]]) + "\n")
    tpm_file.close()
    print("[info] FPKM file: {}".format(tpm_file_path))

    #Write QC Table
    qc_file = open(qc_file_path, "w")
    qc_file.writelines("\n".join(qc))
    qc_file.close()