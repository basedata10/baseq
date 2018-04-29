def write_buffer(buffer, filehandles):
    for key in buffer.keys():
        filehandles[key].writelines(buffer[key])
        buffer[key] = []

def write_and_close(buffer, filehandles):
    for key in buffer.keys():
        filehandles[key].writelines(buffer[key])
        filehandles[key].close()

def filter_polyAT(fastq1, fastq2, outname):
    """
    filter reads with polyA30 or polyT30 sequences...
    """
    #Iterating the reads...
    from baseq.utils.file_reader import read_file_by_lines
    buffers1 = []
    buffers2 = []
    # file1 = gzip.open(outname + ".1.fq.gz", 'wb')
    # file2 = gzip.open(outname + ".2.fq.gz", 'w')
    file1 = open(outname + ".1.fq.gz", 'wb')
    file2 = open(outname + ".2.fq.gz", 'wb')

    inlines1 = read_file_by_lines(fastq1, 1000 * 1000 * 1000, 4)
    inlines2 = read_file_by_lines(fastq2, 1000 * 1000 * 1000, 4)
    count = 0
    seqA = "".join(['A']*30)
    seqT = "".join(['T']*30)
    print(seqA, seqT)
    for line1 in inlines1:
        line2 = next(inlines2)
        count += 1
        if count % 100000 == 1:
            print("[info] 10K reads")
            file1.write("".join(buffers1).encode('utf-8'))
            file2.write("".join(buffers2).encode('utf-8'))
            buffers1 = []
            buffers2 = []
        if seqA in line1[1] or seqT in line1[1] or seqA in line2[1] or seqT in line2[1]:
        #if re.search('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA|TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', line1[1]) or re.search('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA|TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', line2[1]):
            pass
        else:
            buffers1.append("".join(line1))
            buffers2.append("".join(line2))

    file1.writelines(buffers1)
    file2.writelines(buffers2)
    file1.close()
    file2.close()