from baseq.utils.file_reader import read_file_by_lines

def filter_fastq_pair_by_sequence(fastq1, fastq2, seqfile, samplename):
    #read seqfile or template file
    with open(seqfile, "r") as file:
        seqs = [line.strip() for line in file.readlines()]
    templates = [seq for seq in seqs if len(seq)>5]
    print("[info] The sequence is : {}".format(templates))

    buffers1 = []
    buffers2 = []
    file1 = open(samplename + ".1.fq", 'wb')
    file2 = open(samplename + ".2.fq", 'wb')
    inlines1 = read_file_by_lines(fastq1, 1000 * 1000 * 1000, 4)
    inlines2 = read_file_by_lines(fastq2, 1000 * 1000 * 1000, 4)
    count = 0
    valid = 0
    for line1 in inlines1:
        line2 = next(inlines2)
        count += 1
        if count % 1000000 == 0:
            print("[info] {}M reads for {}".format(count / 1000000, samplename))
        if count % 100000 == 1:
            file1.write("".join(buffers1).encode('utf-8'))
            file2.write("".join(buffers2).encode('utf-8'))
            buffers1 = []
            buffers2 = []

        line_skip = 0
        for seq in templates:
            if seq in line1[1] or seq in line2[1]:
                line_skip = 1
                continue

        if line_skip:
            continue

        valid += 1
        buffers1.append("".join(line1))
        buffers2.append("".join(line2))

    file1.writelines(buffers1)
    file2.writelines(buffers2)
    file1.close()
    file2.close()
    return [samplename, count, valid]