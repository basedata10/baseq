import os
from basematic.Operate.Bash import parseBash
from basematic.setting import bash_script_dir

class GATK:
    def __init__(self, configs):
        self.configs = configs

    def check(self):
        pass

    def generate_script(self):
        path = os.path.join(bash_script_dir, "GATK.sh")
        template = parseBash(path)
        return "\n".join(["{}={}".format(key, self.configs[key]) for key in self.configs]) + "\n" + template