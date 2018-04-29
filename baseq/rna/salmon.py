import sys, os, time, re, json
from subprocess import Popen, PIPE, call
from baseq.mgt import get_config, run_cmd
from baseq.fastq.sample_file import check_sample_files

import pandas as pd

def run_build_salmon_index():
    pass

def run_salmon(fq1, fq2, genome, outdir):
    print(fq1, fq2, genome, outdir)
    salmon = get_config("RNA", "salmon")
    salmon_ref = get_config("RNA_ref_"+genome, "salmon_index")
    gene_map = get_config("RNA_ref_"+genome, "gene_map")

    # Run salmon
    if fq1 and fq2 and os.path.exists(fq1) and os.path.exists(fq2):
        salmon_cmd = [salmon, 'quant', '-i', salmon_ref, '-l A', '-1', fq1, '-2', fq2, '-p 8', '-g', gene_map, '-o', outdir]
    elif fq1 and os.path.exists(fq1):
        salmon_cmd = [salmon, 'quant', '-i', salmon_ref, '-l A', '-r', fq1, '-p 8', '-g', gene_map, '-o', outdir]
    else:
        pass

    run_cmd("Salmon Quantification", " ".join(salmon_cmd))

def run_multiple_salmons(path, genome, outdir):
    samples = check_sample_files(path)
    print(samples)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("[info] Create outdir in: {}".format(outdir))
    script_lists = []
    script_main_path = os.path.join(outdir, "work.sh")
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
            script = "baseq-RNA run_salmon -1 {} -2 {} -g {} -d {}".format(sample[1], sample[2], genome, path)
        else:
            script = "baseq-RNA run_salmon -1 {} -g {} -d {}".format(sample[1], genome, path)
        with open(script_path, "w") as file:
            file.writelines(script+"\n")
            print("[info] work script written in {}".format(script_path))
        script_lists.append(script_path)
    #write the main script
    with open(script_main_path, "w") as file:
        file.writelines("\n".join(["bash " + x for x in script_lists]))
    print("[info] Main script written in {}".format(script_main_path))


# Build TPM and QC ...
def build_tpm_table(processdir, samplefile, outpath):
    samples = check_sample_files(samplefile)
    sample_names = [sample[0] for sample in samples]
    tpm_file_path = outpath + "tpm.txt"
    count_file_path = outpath + "count.txt"
    qc_file_path = outpath + "qc.txt"
    tpm = {}
    count = {}
    qc = ["\t".join(['sample', 'reads', 'mapped', 'ratio', 'genecounts'])]
    for sample in sample_names:
        #build TPM table
        sample_TPM = []
        salmon_gene_path = os.path.join(processdir, sample, 'quant.genes.sf')
        with open(salmon_gene_path, 'r') as infile:
            infile.readline()
            for line in infile:
                infos = re.split("\t", line)
                gene = infos[0]
                sample_TPM.append(float(infos[3]))
                if not gene in tpm:
                    tpm[gene] = [float(infos[3])]
                    count[gene] = [float(infos[4])]
                else:
                    tpm[gene].append(float(infos[3]))
                    count[gene].append(float(infos[4]))

        #Genes Detected
        genes_TPM_1 = sum([1 for x in sample_TPM if x>=1])

        #build QC data
        metainfo = os.path.join(processdir, sample, 'aux_info', 'meta_info.json')
        with open(metainfo, "r") as file:
            qc_info = json.load(file)
            qc_sample = [sample, qc_info["num_processed"], qc_info["num_mapped"], qc_info["percent_mapped"], genes_TPM_1]
            qc.append("\t".join([str(x) for x in qc_sample]))

    #Write TPM
    tpm_file = open(tpm_file_path, "w")
    tpm_file.write("\t".join(["gene"] + sample_names) + "\n")
    for gene in tpm:
        sum_tpm_genes = sum(tpm[gene])
        if sum_tpm_genes > 0:
            tpm_file.write("\t".join([gene] + [str(x) for x in tpm[gene]]) + "\n")
    tpm_file.close()

    #Write Counts
    count_file = open(count_file_path, "w")
    count_file.write("\t".join(["gene"] + sample_names) + "\n")
    for gene in count:
        sum_tpm_genes = sum(count[gene])
        if sum_tpm_genes > 0:
            count_file.write("\t".join([gene] + [str(x) for x in count[gene]]) + "\n")
    count_file.close()

    #Write QC Table
    qc_file = open(qc_file_path, "w")
    qc_file.writelines("\n".join(qc))
    qc_file.close()