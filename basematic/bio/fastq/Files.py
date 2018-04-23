import os, sys

def check_infiles(multiple="", name="", fq1="", fq2=""):
    samples = []
    if multiple:
        print("[info] use multiple mode, the name, fq1, fq2")
        if os.path.exists(multiple):
            print("[info] Process Samples Infos in File {}".format(multiple))
            with open(multiple, 'r') as infile:
                lines = infile.readlines()
            for line in lines:
                info = line.split()
                sample = info[0]
                fq1 = info[1]
                fq2 = info[2]
                if os.path.exists(fq1):
                    if os.path.exists(fq2):
                        print("[info] '{}' is Paired End Sample, fq1: '{}' fq2:'{}'".format(sample, fq1, fq2))
                        samples.append([sample, os.path.abspath(fq1), os.path.abspath(fq2)])
                    else:
                        print("[info] '{}' is Single End Sample, fq1: '{}'".format(sample, fq1))
                        samples.append([sample, os.path.abspath(fq1), ''])
                else:
                    print("[Exit] No valid file for {}".format(sample))
        else:
            sys.exit("[Error] Sample Info File do not exist {}".format(multiple))

    else:
        sample = name
        if os.path.exists(fq1):
            if os.path.exists(fq2):
                print("[info] Paired End Sample, name:{} fq1:{} fq2:{}".format(sample, fq1, fq2))
                samples.append([sample, os.path.abspath(fq1), os.path.abspath(fq2)])
            else:
                print("[info] Single End Sample, name:{} fq1:{}".format(sample, fq1))
                samples.append([sample, os.path.abspath(fq1), ''])
        else:
            print("[Warning] No valid file for {}".format(sample))

    return samples

def List_Fastq_InFiles(path):
    pass