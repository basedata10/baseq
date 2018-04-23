import basematic
import click
import os


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.option('--name', default='', help='Who are you?')

def cli(name):
    pass

@cli.command(short_help = "Show the status of jobs")
@click.argument('path')
def list_fastq(path):
    click.echo('Start Using CNV Analysis', path)

@cli.command(short_help = "Show the Softs")
def softs():
    click.echo('Show the softs')

@cli.command(short_help = "Show the status of jobs")
@click.argument('path')
def list_reference(path):
    click.echo('Start Using CNV Analysis')

@cli.command(short_help="Qsub the Jobs In The File")
@click.option('--thread', '-t', default='1', help='Threads used for each command (1)')
@click.option('--mem', '-m', default='2G', help='Memory (2G)')
@click.argument('script_path')

def QSUB(script_path, thread, mem):
    print(script_path, thread, mem)

@cli.command(short_help="Serve the current Dir")
@click.option('--port', default='8777', help='The Port')
def server(port):
    try:
        from http.server import test
        test(port = int(port))
    except:
        import SimpleHTTPServer
        import SocketServer
        PORT = int(port)
        Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
        Handler.extensions_map.update({
            '.webapp': 'application/x-web-app-manifest+json',
        })
        httpd = SocketServer.TCPServer(("", PORT), Handler)
        httpd.serve_forever()