from subprocess import call
import sys

def run_cmd(name, cmd):
    print("[run] Command: {}".format(cmd))
    try:
        call(cmd, shell=True)
        print("[info] {} complete.".format(name, cmd))
    except:
        sys.exit("[error] {} exit with error.".format(name))