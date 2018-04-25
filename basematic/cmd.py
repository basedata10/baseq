import click, os, sys

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)

def cli():
    import subprocess
    cmd = 'eval "$(_BASEMATIC_COMPLETE=source basematic)"'
    subprocess.call(cmd, shell=True)
    click.echo("Welcome to Basematic.")

@cli.command(short_help = "Show the status of jobs")
def init():
    import basematic.setting as setting
    from shutil import copyfile
    if os.path.exists(setting.DIR):
        print("[info] Config Folder Exits")
    else:
        print("[info] Config Folder Do Not Exits, make one")
        try:
            os.mkdir(setting.DIR)
            print("[info] Success")
        except:
            sys.exit("[error] Failed to Create the new Dir")
    if os.path.exists(setting.Config):
        print("[info] Config File Exits {}".format(setting.Config))
    else:
        print("[info] Config File Exits Dose not Exist, Create one from template")
        copyfile(setting.Config_template, setting.Config)
        print("[info] The config file is ready: {}".format(setting.Config))

@cli.command(short_help = "Show the status of jobs")
@click.option('--port', '-p', default='./', help='port')
@click.argument('name')
def looker(name, port):
    click.echo('Start Using CNV Analysis')

@cli.command(short_help="Show the status of jobs")
@click.option('--port', '-p', default='./', help='port')
def serve(port):
    import sys, subprocess
    import socket
    ipAdd = socket.gethostbyname(socket.gethostname())
    print("Visit: {}:{}".format(ipAdd, port))
    if sys.version_info[0] < 3:
        subprocess.call("python -m SimpleHTTPServer {}".format(port), shell=True)
    else:
        subprocess.call("python -m http.server {}".format(port), shell=True)

@cli.command(short_help="mount_disk/")
@click.argument("command")
def notes(command):
    print("lsblk -f")
    print("mount /dev/sdb2 /mnt/path")