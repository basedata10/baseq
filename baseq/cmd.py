import click, os, sys
from baseq.mgt.command import run_cmd
from baseq.docs import docs
from baseq.mgt.history import add_record

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.pass_context
def cli(ctx):
    add_record(dir)
    print(ctx.obj, dir(ctx))
    print(ctx.command_path, ctx.token_normalize_func, ctx.args, ctx.lookup_default)

@cli.command(name='init', short_help = "init the configure files")
def init():
    import baseq.setting as setting
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

@cli.command(short_help = "Print The Docs")
def docs():
    print(docs)

@cli.command(short_help = "Command History")
@click.option('--current_dir', '-c', is_flag=True)
def hist(current_dir):
    from baseq.mgt.history import history_cmd
    if current_dir:
        print("Current Dir Is...", os.getcwd())
    else:
        history_cmd()

@cli.command(short_help = "Monitor the Status Of Server")
@click.argument('name')
def monitor(name):
    pass

@cli.command(short_help = "Show the status of jobs")
@click.option('--port', '-p', default='./', help='port')
@click.argument('name')
def looker(name, port):
    pass

@cli.command(short_help = "List The Samples and Files")
@click.argument('path')
@click.argument('filepattern')
def list_sample_files(path, filepattern):
    from jinja2 import Template
    cmd = Template(r"""ls -l {{path}}/*/*{{filepattern}}* | awk '{print $9}' | perl -ne '@a=split(/\//,$_);print "$a[-2]\t$_"' """).render(path=path, filepattern=filepattern)
    print(cmd)
    run_cmd("", cmd)

@cli.command(short_help="Show the status of jobs")
@click.option('--port', '-p', default='7878', help='port (7878)')
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

@cli.command(short_help="Deploy the package to mnt path...")
@click.argument("command", nargs=-1)
def deploy(command):
    print("[info] Deploy the packages to the mnt server...")

@cli.command(short_help="PPT...")
def PPT():
    from baseq.utils.ppt import generate_PPT
    generate_PPT()

@cli.command(short_help="Pack folder into tar.gz")
@click.argument("path")
def pack(path):
    name = os.path.basename(path)
    run_cmd("PACK", "tar -zcvf {}.tar.gz {}".format(name, path))
    print("[info] {}.tar.gz is generated...".format(name))

@cli.command(short_help="tools for CSV files, run 'baseq csv' for help")
@click.option('--subject', '-s', default='A Email', help='Subject')
@click.option('--message', '-m', default='The message', help='Message')
@click.option('--attches', '-a', multiple=True, default='', help='Message')
def email(subject, message, attches):
    from baseq.utils.email_send import Client
    print(subject, message, attches)
    Client().send_mail(subject, message, attches)