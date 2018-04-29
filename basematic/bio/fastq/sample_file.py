import os, sys, re

def check_sample_files(samplefile="", name="", fq1="", fq2=""):
    """Check sample and fastq paths from samplefile and fqs"""
    samples = []
    if samplefile:
        print("[info] use multiple mode, the name, fq1, fq2")
        if os.path.exists(samplefile):
            print("[info] Process Samples Infos in File {}".format(samplefile))
            with open(samplefile, 'r') as infile:
                lines = infile.readlines()
            for line in lines:
                info = line.split()
                #Less than 2 columns
                if len(info)<2:
                    print("[warning] Line do not contains sample ...")
                sample = info[0]
                fq1 = info[1]
                if len(info)>2:
                    fq2 = info[2]
                else:
                    fq2 = ""
                if fq1 and os.path.exists(fq1):
                    if fq2 and os.path.exists(fq2):
                        print("[info] '{}' is Paired End Sample, fq1: '{}' fq2:'{}'".format(sample, fq1, fq2))
                        samples.append([sample, os.path.abspath(fq1), os.path.abspath(fq2)])
                    else:
                        print("[info] '{}' is Single End Sample, fq1: '{}'".format(sample, fq1))
                        samples.append([sample, os.path.abspath(fq1), ''])
                else:
                    print("[Exit] No valid file for {}".format(sample))
        else:
            sys.exit("[Error] Sample Info File do not exist {}".format(samplefile))
    #Check Single Sample...
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

def list_fastq_files(sampleDir, writeFile):
    files = [file for file in os.listdir(sampleDir) if os.path.isfile(os.path.join(sampleDir, file))]
    out = []
    samples = {}
    print("[info] The file name should be <sample>_ANY_INFOS_1.fq.gz or <sample>_ANY_INFOS_2.fq.gz")
    for file in files:
        try:
            infos = file.split("_")
            if infos[-1] in ["1.fq.gz", "2.fq.gz", "1.fastq.gz", "2.fastq.gz", "R1.fq.gz", "R2.fq.gz",
                             "1.clean.fq.gz", "2.clean.fq.gz"]:
                out.append([infos[0], infos[-1], file])
        except:
            pass

    print("[info] Detect {} files in path {}".format(len(out), sampleDir))
    print("[info] Build SE/PE samples infos from file lists")
    for sample in out:
        name = sample[0]
        path = sample[2]
        if name not in samples:
            samples[name] = [path]
        else:
            samples[name].append(path)

    lines = []
    for name in sorted(samples.keys()):
        print("\t".join([name] + sorted(samples[name])))
        lines.append("\t".join([name] + sorted([os.path.join(sampleDir, p) for p in samples[name]])))

    print("[info] ##############################")
    print("[info] ##### DECTED {} samples ######".format(len(lines)))
    print("[info] ##############################")
    try:
        with open(writeFile, 'w') as file:
            file.writelines("\n".join(lines))
        print("[info] Write the sample infos to {}".format(writeFile))
    except:
        sys.exit("[error] Failed to Save the sample files")

    return out