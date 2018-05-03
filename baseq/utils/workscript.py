import os

def write_bash(dir, script_paths, name="work.sh"):
    path = os.path.join(dir, name)
    with open(path, "w") as file:
        file.writelines("\n".join(["bash " + x for x in script_paths]))
    print("[info] Main script written in {}".format(path))

def write_bash_qsub(dir, script_paths, name="work_qsub.sh"):
    path = os.path.join(dir, name)
    with open(path, "w") as file:
        file.writelines("\n".join(["qsub -cwd -l vf=8g,p=8 " + x for x in script_paths]))
    print("[info] Main script written in {}".format(path))