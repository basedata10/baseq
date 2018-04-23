import subprocess
from jinja2 import Template

def parseBash(path):
    with open(path, 'r') as infile:
        lines = infile.readlines()
    for idx, line in enumerate(lines):
        if line.startswith("##"):
            command = line[2:].rstrip()
            echo_start = '''echo "$(date +"%m-%d %T") {} START" >>log.txt\n'''.format(command)
            echo_end = '''echo "$(date +"%m-%d %T") {} SUCSESS" >>log.txt\n'''.format(command)
            check_success = Template(
                '''rc=$?; if [[ $rc != 0 ]]; then echo "$(date +"%m-%d %T") {{cmd}} Failed, Exit the Script" >>log.txt; exit $rc; fi''').render(
                cmd=command)
            lines[idx] = echo_start
            lines[idx + 1] = lines[idx + 1] + check_success + "\n" + echo_end
    return "".join(lines)

def run_cmd(cmd):
    try:
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    except:
        return []
    result = [x.decode("utf-8").strip() for x in process.stdout.readlines()]
    return result