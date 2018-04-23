import click, sys
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)

def cli():
    click.echo("Basematic-SNP Pipeline Started")

@cli.command(short_help="Generate Pipeline")
@click.option('--multiple', '-m', default='', help="Tab seprated file: name, fq1, fq2, if set, the name, fq1, fq2 options is disabled")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--fq1', '-1', default='', help='fastq1 path')
@click.option('--fq2', '-2', default='', help='fastq2 path (optional)')
@click.option('--dir', '-d', default='./', help='Output Path (./)')
@click.option('--interval', '-i', default='', help='Interval file, chrN:start_end')
@click.option('--job', default='local', help='job submit mode', type=click.Choice(['local', 'qsub', 'bsub']))

def run_all(multiple, name, fq1, fq2, interval, dir, job):
    from basematic.mgt.config import configManager
    from basematic.bio.fastq import Files as FQfiles
    from basematic.bio.snv.gatk import GATK
    from basematic.utils.clean import cleanStr
    from basematic.file.Folder import EnsurePath, WriteFile
    from jinja2 import Template
    import os, sys

    outDir = os.path.abspath(dir)

    #Check configurations
    from basematic.mgt.check import isExist
    print("[info] Check the config files")

    cfgs = configManager().get_section("SNV")
    for cfg in cfgs.data:
        isExist(cfgs.data[cfg], cfg)
    script_template = GATK(cfgs.data).generate_script()

    #Check Annovar
    from basematic.bio.snv.annovar import Annovar
    cfg_annovar = configManager().get_section("Annovar")
    annovar = Annovar(cfg_annovar.get("Annovar"), cfg_annovar.get("annovar_db_hg38"), "hg38")

    annovar.check_script()

    #Check Intervals
    interval = os.path.abspath(interval)
    EnsurePath(interval, "Interval")
    EnsurePath(outDir, "OutPut")

    #Check Samples
    samples = FQfiles.check_infiles(multiple, name, fq1, fq2)
    if len(samples) == 0:
        sys.exit("[error] No valid samples")
    else:
        print("[info] {} samples detected".format(len(samples)))

    print('[info] Generate Scripts ...')
    script_paths = []

    for sample in samples:
        name = cleanStr(sample[0])
        fq1 = sample[1]
        fq2 = sample[2]
        sample_dir = os.path.join(outDir, name)
        EnsurePath(sample_dir, "Sample '{}' Result".format(name))
        script = Template(script_template).render(outdir = sample_dir, name=name, fq1=fq1, fq2=fq2, interval=interval)
        script_path = os.path.join(sample_dir, "work_GATK_{}.sh".format(name))
        WriteFile(script_path, script, "Script of {}".format(name))
        script_paths.append(script_path)

    if job == "local":
        work = "\n".join(["bash {}".format(line) for line in script_paths])
    elif job == "qsub":
        work = "\n".join(["qsub -cwd -l vf=10g {}".format(line) for line in script_paths])
    elif job == "bsub":
        work = "\n".join(["bsub {}".format(line) for line in script_paths])
    else:
        pass

    work_path = os.path.join(outDir, "work.sh")
    WriteFile(work_path, work, "Main Script".format(sample))
    print("[Success] Generate Work Script in {}".format(work_path))
    print("[Finish] You can run or submit by running: bash {}".format(work_path))

@cli.command(short_help="Run Annovar")
@click.option('--vcf', '-i', default='', help='Input VCF File')
@click.option('--out', '-o', default='./Annovar.vcf', help='输出VCF')
@click.option('--genome', '-g', default='', help='Genome Version')
def run_annovar(vcf, out, genome):
    from basematic.mgt import config
    from basematic.bio.snv.annovar import Annovar
    cfgs = config.configManager().get_section("Annovar")
    print(cfgs)

@cli.command(short_help="Filter the VCF")
@click.argument("path")
@click.argument("outpath")
def filter(path, outpath):
    from basematic.bio.snv.vcf import Annovar_Filter
    click.echo('Filter the VCF...')
    res = Annovar_Filter(path).filter_1KGenome()
    res.to_csv(outpath,sep="\t")

@cli.command(short_help="Check the enrichment quality, input: bam, interval and outpath")
@click.argument("bampath")
@click.argument("interval")
@click.argument("outpath")
def QC_enrich(bampath, interval, outpath):
    from .quality import QC_enrich
    QC_enrich(bampath, interval, outpath)