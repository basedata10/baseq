import subprocess

def run_it(command):
    reading = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, bufsize=1000000)
    infile = reading.stdout
    return [line.decode('utf8').strip() for line in infile.readlines()]

def run_generator(command):
    reading = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, bufsize=1000000)
    infile = reading.stdout
    for line in infile.readlines():
        yield line.decode('utf8').strip()