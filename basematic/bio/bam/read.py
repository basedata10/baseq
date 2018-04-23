import subprocess

def bowtie2_unique_pos(filepath):
    reading = subprocess.Popen(["samtools", "view", filepath], stdout=subprocess.PIPE, bufsize=1000000)
    infile = reading.stdout
    while True:
        counter = counter + 1
        data = infile.readline().decode('utf8')
        if data[0] == "":
            break
        yield data