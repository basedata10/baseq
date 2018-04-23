import os

HOME = os.path.expanduser("~")
DIR = os.path.join(HOME, "basematic")
Config = os.path.join(DIR, "config.ini")

INSTALLED_DIR = os.path.dirname(os.path.abspath(__file__))
Config_template = os.path.join(INSTALLED_DIR, "config.ini")

Rdir = os.path.join(INSTALLED_DIR, "ScriptR", "src")
bash_script_dir = os.path.join(INSTALLED_DIR, "ScriptShell")
perl_script_dir = os.path.join(INSTALLED_DIR, "src", "perl")
r_script_dir = os.path.join(INSTALLED_DIR, "src", "r")