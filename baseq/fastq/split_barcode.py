import os

def write_buffer(buffer, filehandles):
    for key in buffer.keys():
        filehandles[key].writelines(buffer[key])
        buffer[key] = []

def write_and_close(buffer, filehandles):
    for key in buffer.keys():
        filehandles[key].writelines(buffer[key])
        filehandles[key].close()

def split_barcode(barcode_file, fastq, outdir, suffix):
    """
    barcode_file: tsv: samplename barcode_string...
    """
    with open(barcode_file, 'r') as file:
        barcodes = file.readlines()
    bc_files = {}
    bc_buffers = {}
    for bc in barcodes:
        info = [x.strip() for x in bc.split()]
        if len(info) == 2:
            sample = info[0]
            barcode = info[1]
            print("[info] The barcode for {} is {}".format(info[0], info[1]))
            path = os.path.join(outdir, sample) + "." + suffix
            print("[info] File for '{}' is '{}'".format(sample, path))
            bc_files[barcode] = open(path, "w")
            bc_buffers[barcode] = []

    #Iterating the reads...
    from baseq.utils.file_reader import read_file_by_lines
    inlines = read_file_by_lines(fastq, 1000 * 1000 * 1000, 4)
    count = 0
    for line in inlines:
        count += 1
        barcode = line[0].strip().split(":")[-1]
        if barcode in bc_files:
            bc_buffers[barcode].append("".join(line))
        if count % 1000000 == 1:
            write_buffer(bc_buffers, bc_files)

    write_and_close(bc_buffers, bc_files)