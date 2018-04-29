import subprocess

def read_file_by_lines(filepath, maxLines, linecount):
    counter = 0
    if filepath.endswith("gz") or filepath.endswith("gzip") or filepath.endswith("gz2"):
        reading = subprocess.Popen(["gunzip", "-c", filepath], stdout=subprocess.PIPE, bufsize=1000000)
        infile = reading.stdout
        while True:
            counter = counter + 1
            if counter > maxLines:
                break
            data = [infile.readline().decode('utf8') for i in range(linecount)]
            if data[0] == "":
                break
            yield data
    else:
        infile = open(filepath, 'r')
        while True:
            counter = counter + 1
            if counter > maxLines:
                break
            data = [infile.readline() for i in range(linecount)]
            if data[0] == "":
                break
            yield data

class Batch_Writer:
    def __init__(self, filepath, batchSize=100000):
        self.outfile = open(filepath, 'wt')
        self.index = 0
        self.lineBuffer = []
        self.batchSize = batchSize

    def add_line(self, line):
        self.index += 1
        if self.index % self.batchSize == 0:
            self.outfile.writelines(self.lineBuffer)
            self.lineBuffer = []
        else:
            self.lineBuffer.append(line)

    def on_finish(self):
        self.outfile.writelines(self.lineBuffer)