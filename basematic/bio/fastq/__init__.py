import click, os, sys

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)

def cli():
    click.echo("Basematic Fastq...")

@cli.command(short_help="Quality Control")
@click.option('--samplefile', '-m', default='', help='sample file lists')
@click.option('--fq1', '-1', default='', help='fastq1')
@click.option('--fq2', '-2', default='', help='fastq2')
@click.option('--out', '-o', default='./qc.txt', help='file path (./qc.txt)')
def QC(samplefile, fq1, fq2, out):
    import pandas as pd
    print("[info] Qualtity control of the fastq files ...")
    print("[info] The result file will be write to {}".format(out))
    from basematic.bio.fastq.quality import fastq_basecontent_quality
    from .samplefile import check_sample_files
    result = []

    if samplefile:
        print("[info] Read the samples infos the file...")
        samples = check_sample_files(samplefile)
        print(samples)
    else:
        if fq1 and os.path.exists(fq1):
            result.append(fastq_basecontent_quality(fq1))
        else:
            sys.exit("[error] fastq1 is not a valid path, exit...")
        if fq2:
            if os.path.exists(fq2):
                result.append(fastq_basecontent_quality(fq2))
            else:
                sys.exit("[error] fastq2 is not a valid path, exit...")
    #pd.concat(result).to_csv(out, sep="\t")

#sampling from a fastq file ...
@cli.command(short_help="Quality Control")
@click.option('--path', '-d', default='./', help='File Path')
@click.option('--ratio', '-d', default='1', help='ratio of reads')
@click.option('--reads', '-d', default='1', help='read counts')
@click.option('--tofasta', '-d', default='1', help='read counts')
@click.option('--path', '-d', default='./', help='File Path')
def sampling(path):
    from basematic.bio.fastq.quality import fastq_basecontent_quality
    print(fastq_basecontent_quality(path))

@cli.command(short_help="List all the Fastq Files in Path (not include the subdir)")
@click.argument("path")
@click.option('--subdir', '-s', default='', help='The files are in the subdirectory')
@click.option('--write', '-f', default='./samples.txt', help='Write sample infos to this path (./samples.fqs.txt)')
def list_samples(path, write, subdir):
    print("[info] List fastq files in {}".format(path))
    from .samplefile import list_fastq_files
    list_fastq_files(path, write)