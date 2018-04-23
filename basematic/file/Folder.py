from jinja2 import Template
from basematic.exec.bash import run_cmd
import pandas as pd
import os, sys

def getHomeDir():
    home = os.path.expanduser("~")
    return home

def EnsurePath(path, name =""):
    if not os.path.exists(path):
        print("[info] {} Folder do not exists, try to create : {}".format(name, path))
        try:
            os.mkdir(path)
            print("[info] {} Path created : {}".format(name, path))
        except:
            sys.exit("[error] {} created Failed".format(name))
    else:
        print("[info] {} Path Already Exists".format(name))

def WriteFile(path, content="", info="Content", append=False):
    try:
        with open(path, 'w') as file:
            file.writelines(content)
            print("[info] {} Successfully Write to path {}".format(info, path))
    except:
        sys.exit("[error] {} Failed to writing to file {}".format(info, path))

def List_Folder(path):
    try:
        cmd = Template("ls -lh {{path}} | awk '{print $1, $5, $9}'").render(path = path)
        files = run_cmd(cmd)
        results = []
        for file in files:
            infos = file.split()
            try:
                if infos[0][0] == "d":
                    results.append(["folder", infos[2], ""])
                else:
                    results.append(["file", infos[2], infos[1]])
            except:
                pass
    except:
        print('[error]: Failed to List the Folders ...')
        results = []
    return pd.DataFrame(results, columns=["Type", "Name", "Size"]).sort_values(by=['Type'], ascending=False)

if __name__ == "__main__":
    print(List_Folder("./").to_html(index=False))