import subprocess

def run_it(command):
    reading = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, bufsize=1000000)
    infile = reading.stdout
    return [line.decode('utf8').strip() for line in infile.readlines()]