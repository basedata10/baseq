import click
from baseq.snv import cli

@cli.command(short_help="Generate Pipeline")
@click.option('--multiple', '-m', default='', help="Tab seprated file: name, fq1, fq2, if set, the name, fq1, fq2 options is disabled")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--fq1', '-1', default='', help='fastq1 path')
@click.option('--fq2', '-2', default='', help='fastq2 path (optional)')
@click.option('--dir', '-d', default='./', help='Output Path (./)')
@click.option('--interval', '-i', default='', help='Interval file, chrN:start_end')
@click.option('--job', default='local', help='job submit mode', type=click.Choice(['local', 'qsub', 'bsub']))

def run_all(multiple, name, fq1, fq2, interval, dir, job):
    from baseq.mgt.config import configManager
    from baseq.fastq import sample_file as FQfiles
    from baseq.snv.gatk import GATK
    from baseq.utils.clean import cleanStr
    from baseq.utils.folder import ensure_path, write_file
    from jinja2 import Template
    import os, sys

    outDir = os.path.abspath(dir)

    #Check configurations
    from baseq.mgt.check import isExist
    print("[info] Check the config files")

    cfgs = configManager().get_section("SNV")
    for cfg in cfgs.data:
        isExist(cfgs.data[cfg], cfg)
    script_template = GATK(cfgs.data).generate_script()

    #Check Annovar
    from baseq.snv.annovar import Annovar
    cfg_annovar = configManager().get_section("Annovar")
    annovar = Annovar(cfg_annovar.get("Annovar"), cfg_annovar.get("annovar_db_hg38"), "hg38")

    annovar.check_script()

    #Check Intervals
    interval = os.path.abspath(interval)
    ensure_path(interval, "Interval")
    ensure_path(outDir, "OutPut")

    #Check Samples
    samples = FQfiles.check_sample_files(multiple, name, fq1, fq2)
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
        ensure_path(sample_dir, "Sample '{}' Result".format(name))
        script = Template(script_template).render(outdir = sample_dir, name=name, fq1=fq1, fq2=fq2, interval=interval)
        script_path = os.path.join(sample_dir, "work_GATK_{}.sh".format(name))
        write_file(script_path, script, "Script of {}".format(name))
        script_paths.append(script_path)

    # if job == "local":
    #     work = "\n".join(["bash {}".format(line) for line in script_paths])
    # elif job == "qsub":
    #     work = "\n".join(["qsub -cwd -l vf=10g {}".format(line) for line in script_paths])
    # elif job == "bsub":
    #     work = "\n".join(["bsub {}".format(line) for line in script_paths])
    # else:
    #     pass

    work_path = os.path.join(outDir, "work.sh")
    write_file(work_path, work, "Main Script".format(sample))
    print("[Success] Generate Work Script in {}".format(work_path))
    print("[Finish] You can run or submit by running: bash {}".format(work_path))

@cli.command(short_help="Run Annovar")
@click.option('--vcf', '-i', default='', help='Input VCF File')
@click.option('--out', '-o', default='./Annovar.vcf', help='输出VCF')
@click.option('--genome', '-g', default='', help='Genome Version')
def run_annovar(vcf, out, genome):
    from baseq.mgt import config
    cfgs = config.configManager().get_section("Annovar")
    print(cfgs)

@cli.command(short_help="Filter the VCF")
@click.argument("path")
@click.argument("outpath")
def filter(path, outpath):
    from baseq.snv.vcf import Annovar_Filter
    click.echo('Filter the VCF...')
    res = Annovar_Filter(path).filter_1KGenome()
    res.to_csv(outpath,sep="\t")


@cli.command(short_help="Check the enrichment quality, input: bam, interval and outpath")
@click.argument("bampath")
@click.argument("interval")
@click.argument("outpath")
def enrich_saturation(bampath, interval, outpath):
    from .quality import enrich_saturation
    enrich_saturation(bampath, interval, outpath)

#Prepare Web Datas
@cli.command(short_help="Generating the datas for web view: baseq.io/viewcnv")
@click.option('--path', '-p', default='', help="The path to the Process Folder")
@click.option('--vcf_annovar', '-v', default='', help="The path to the Process Folder")
def web_data(path, vcf_annovar):
    from baseq.utils.buildresult import pack_web_datas
    data = pack_web_datas()

    #quality stats...
    from .quality import parse_picard_wes_matrics
    res = parse_picard_wes_matrics(path)
    data.add_image("depth_density", "./density.png")
    for k in res:
        data.data[k] = res[k]

    #SNV Table
    from baseq.snv.vcf import Annovar_Filter
    click.echo('[info] Filter the VCF...')
    res = Annovar_Filter(vcf_annovar).filter_1KGenome()

    data.data["VCF"] = res.to_json(orient='split')
    print(data.data.keys())
    import json
    with open("./SNV.json", "w") as file:
        json.dump(data.data, file)


@cli.command(short_help="bwa alignment")
@click.option('--fq1', '-1', default='', help='fastq1 path')
@click.option('--fq2', '-2', default='', help='fastq2 path (optional)')
@click.option('--name','-n',default='',help='prefix of bamfile')
@click.option('--genome','-g',default='',help='Species hg19 or mm10/mm38')
def run_bwa(fq1,fq2,name,genome):
    from .gatk import run_alignment
    run_alignment(fq1,fq2,name,genome)
    print("[info] alignment complete")



@cli.command(short_help="mark duplicates")
@click.option('--name', '-n', default='', help='prefix of bamfile')
def run_markdup(name):
    from .gatk import run_markdup
    run_markdup(name)
    print("[info] mark duplicates complete")

@cli.command(short_help="Base Recalibrator")
@click.option('--name','-n',default='',help='prefix of bamfile')
@click.option('--genome','-g',default='',help='Species hg38/hg19 or mm10/mm38')
def run_bqsr(name,genome):
    from .gatk import bqsr
    bqsr(name,genome)
    print("[info]base recalibrator complete")

@cli.command(short_help="call variants")
@click.option('--name','-n',default='',help='prefix of bamfile')
@click.option('--genome','-g',default='',help='Species hg19 or mm10/mm38')
def run_callvar(name,genome):
    from .gatk import run_callvar
    run_callvar(name,genome)
    print("[info]call variants complete")

@cli.command(short_help="select variants")
@click.option('--name','-n',default='',help='prefix of bamfile')
@click.option('--genome','-g',default='',help='Species hg19 or mm10/mm38')
def run_selectvar(name,genome):
    from .gatk import selectvar
    selectvar(name,genome)
    print("[info] select variants complete")


@cli.command(short_help="annovar annotation")
@click.option('--name','-n',default='',help='prefix of bamfile')
@click.option('--genome','-g',default='',help='Species hg19 or mm10/mm38')
def run_annovar(name,genome):
    from .annovar import run_annovar
    run_annovar(name,genome)
    print("[info] annovar annotation complete")