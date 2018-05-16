import click, os, sys
from baseq.snv import cli
from baseq.snv.qc.cmd import *

@cli.command(short_help="Check the enrichment quality, input: bam, interval and outpath")
@click.argument("bampath")
@click.argument("interval")
@click.argument("outpath")
def enrich_saturation(bampath, interval, outpath):
    pass

@cli.command(short_help="Bwa alignment")
@click.option('--fq1', '-1', default='', help='fastq1 path')
@click.option('--fq2', '-2', default='', help='fastq2 path (optional)')
@click.option('--name', '-n', default='', help='prefix of bamfile')
@click.option('--thread', '-t', default='8', help='Thread for BWA')
@click.option('--genome','-g', default='', help='Species hg19 or mm10/mm38')
@click.option('--outfile','-o', default='', help='bamfile name')
def run_bwa(fq1, fq2, name, genome, outfile, thread):
    from baseq.snv.gatk import alignment
    alignment(fq1, fq2, name, genome, outfile, thread)
    print("[info] alignment complete")

@cli.command(short_help="mark duplicates")
@click.option('--bamfile', '-b', default='', help='bamfile for gatk analysis')
@click.option('--markedbam', '-m', default='', help='marked bamfile')
def run_markdup(bamfile, markedbam):
    from baseq.snv.gatk import run_markdup
    run_markdup(bamfile, markedbam)
    print("[info] Mark duplicates complete")

@cli.command(short_help="Base Recalibrator")
@click.option('--markedbam', '-m',default='', help='bam file with mark duplicates')
@click.option('--genome', '-g', default='', help='Species hg38/hg19 or mm10/mm38')
@click.option('--keepdup', '-f', default='', help='Set Yes to keep duplicate reads')
@click.option('--bqsrbam', '-q',default='', help='output bam file after base recalibrator')
def run_bqsr(markedbam, bqsrbam, genome, keepdup):
    from baseq.snv.gatk import bqsr
    if keepdup:
        bqsr(markedbam, bqsrbam, genome, disable_dup_filter=True)
    else:
        bqsr(markedbam, bqsrbam, genome)
    print("[info] Base recalibrator complete")

@cli.command(short_help="call variants")
@click.option('--bqsrbam', '-q',default='', help='bam file with base recalibrator')
@click.option('--genome', '-g',default='', help='Species hg19 or mm10/mm38')
@click.option('--rawvcf', '-r',default='', help='output vcf file include snp and indel')
def run_callvar(bqsrbam, rawvcf, genome):
    from baseq.snv.gatk import run_callvar
    run_callvar(bqsrbam, rawvcf, genome)
    print("[info] Call variants complete")

@cli.command(short_help="select variants ")
@click.option('--rawvcf', '-r', default='', help='vcf file include snp and indel')
@click.option('--selectvcf', '-s', default='', help='output file include snp only')
@click.option('--filtervcf', '-f', default='', help='output vcf file after filtering low quality snp')
@click.option('--genome', '-g', default='', help='Species hg19 or mm10/mm38')
def run_selectvar(rawvcf, selectvcf, filtervcf, genome):
    from baseq.snv.gatk import selectvar
    selectvar(rawvcf, selectvcf, filtervcf, genome)
    print("[info] Select variants complete")

@cli.command(short_help="annovar annotation")
@click.option('--filtervcf', '-f', default='', help='vcf file after filtering low quality snp')
@click.option('--genome','-g', default='', help='Species hg19 or mm10/mm38')
@click.option('--dir', '-d', default='', help='folder for output files')
@click.option('--name', '-n', default='', help='prefix of annovar annotation file')
def run_annovar(filtervcf,dir,name,genome):
    from baseq.snv.annovar import run_annovar
    annovarfile = os.path.join(dir, "{}.avinput",format(name))
    run_annovar(filtervcf,annovarfile,name,genome)
    print("[info] Annovar annotation complete")

@cli.command(short_help="run gatk pipeline")
@click.option('--config', '-c', default = '', help = 'config file path')
@click.option('--genome', '-g', default = '', help = 'Species hg19/hg38')
@click.option('--name', '-n', default = '', help = 'sample name')
@click.option('--fq1', '-1', default = '', help = 'fastq1 path')
@click.option('--fq2', '-2', default = '', help = 'fastq2 path(optional)')
@click.option('--bamfile', '-b', default = '', help = 'bamfile for gatk analysis(When fq1,fq2 not set)')
@click.option('--dir', '-d', default = './', help = 'folder for output files')
def run_gatkpipe(config, dir, fq1, fq2, bamfile, name, genome):

    if config:
        if not os.path.exists(config):
            sys.exit("[error] Config file not exists")
        os.environ["BASEQCFG"] = config

    outdir = os.path.abspath(os.path.join(dir, name))
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    from baseq.snv.gatk import alignment, run_markdup, run_callvar, selectvar, bqsr
    from baseq.snv.annovar import run_annovar

    bam_marked = os.path.join(outdir, "{}.marked.bam".format(name))
    bam_bqsr = os.path.join(outdir, "{}.marked.bqsr.bam".format(name))
    vcf_raw = os.path.join(outdir, "{}.raw.snp.indel.vcf".format(name))
    vcf_select = os.path.join(outdir, "{}.raw.snp.vcf".format(name))
    vcf_filter = os.path.join(outdir, "{}.filter.snp.vcf".format(name))
    annovarfile = os.path.join(outdir, "{}.snp.avinput".format(name))

    if fq1 :
        bamfile = os.path.join(outdir, "{}.bam".format(name))
        alignment(fq1, fq2, name, genome, bamfile)
        run_markdup(bamfile, bam_marked)
    if bamfile:
        run_markdup(bamfile, bam_marked)

    bqsr(bam_marked, bam_bqsr, genome)
    run_callvar(bam_bqsr, vcf_raw, genome)
    selectvar(vcf_raw, vcf_select, vcf_filter, genome)
    run_annovar(vcf_filter, annovarfile, name, genome)

def run_multi_gatkpipe():
    pass

@cli.command(short_help="mutect2 analysis")
@click.option('--genome', '-g', default = '', help = 'Species hg19/hg38')
@click.option('--normalbam', '-N', default = '', help = 'bam file of normal sample')
@click.option('--tumorbam', '-T', default = '', help = 'bam file of tumor sample')
@click.option('--normalname', '-n', default = 'normal', help = 'name of normal sample')
@click.option('--tumorname', '-t', default = 'tumor', help = 'name of tumor sample')
@click.option('--path', '-o', default = '', help = 'path of output files')
@click.option('--pon', '-p', default = '', help = 'panel of normal samples')
@click.option('--list', '-l', default = '', help = 'list of normal bamfile and tumor bamfile,normalname,normalbam,tumorname,tumorbam')
@click.option('--germline', '-G', default = '', help = 'population gnomAD')

def run_mutect2(genome, normalbam, normalname, tumorbam, tumorname, list, path, pon, germline):
    from baseq.snv.gatk import mutect2
    if list:
        with open(list, 'r') as file:
            lines = file.readlines()
        lists = [line.strip().split() for line in lines]

    else:
        lists = [[normalname,normalbam,tumorname,tumorbam]]
    import multiprocessing as mp
    pool = mp.Pool(processes=10)
    results = []
    for list in lists:
        vcffile = os.path.join(path,"{}_{}.vcf.gz".format(list[0],list[2]))
        results.append(pool.apply_async(mutect2, (genome, list[0], list[1], list[2], list[3], vcffile, pon, germline)))
    pool.close()
    pool.join()
    [x.get() for x in results]
    print("[info] mutect2 analysis complete")




@cli.command(short_help="mutect2 analysis")
@click.option('--resource', '-r', default = '', help = 'common biallelic vcf files')
@click.option('--tumorbam', '-T', default = '', help = 'bam file of tumor sample')
@click.option('--somaticvcf', '-s', default = '', help = 'mutect call')
@click.option('--path', '-o', default = '', help = 'path of output files')
@click.option('--tumorname', '-t', default = '', help = 'tumor sample name')
def filter_mutect_call(tumorname, tumorbam, somaticvcf, resource, path):
    from baseq.snv.gatk import get_filter_table,filter_mutect_vcf
    gps_table = os.path.join(path,"{}_gps.table".format(tumorname))
    calcontam_table = os.path.join(path, "{}_calcontam.table".format(tumorname))
    filter_call = os.path.join(path, "{}_filter.vcf.gz".format(tumorname))

    get_filter_table(tumorbam, resource, gps_table, calcontam_table)
    filter_mutect_vcf(somaticvcf, calcontam_table, filter_call)


@cli.command(short_help="test for filter_mutect_vcf")
@click.option('--calcontam_table', '-c', default = '', help = 'table')
@click.option('--somaticvcf', '-s', default = '', help = 'mutect call')
@click.option('--path', '-o', default = '', help = 'path of output files')
def filter_mutect_vcf(somaticvcf, calcontam_table, path):
    from baseq.snv.gatk import filter_mutect_vcf
    filter_call = os.path.join(path, "somatic_filter.vcf.gz")
    filter_mutect_vcf(somaticvcf, calcontam_table, filter_call)


@cli.command(short_help="create PoN")
@click.option('--path', '-p', default = '', help = 'path of the output files')
@click.option('--genome', '-g', default = '', help = 'Species hg19/hg38')
@click.option('--list', '-l', default = '', help = 'list of normal bamfile and tumor bamfile,normalname,normalbam,tumorname,tumorbam')
@click.option('--interval', '-L', default = '', help = 'interval list')
def create_pon(genome, list, path, interval):
    from baseq.snv.gatk import create_pon
    create_pon(genome, list, path, interval)

@cli.command(short_help="creat args file for PoN")
@click.option('--path', '-p', default = '', help = 'path of all normal-call vcf files')
def listofvcf(path):
    from baseq.snv.gatk import listofvcf
    listofvcf(path)

@cli.command(short_help="create panel of normals")
@click.option('--outpath', '-o', default = '', help = 'path of output file')
@click.option('--ponpath', '-p', default = '', help = 'path of all normal-call vcf files')
def run_createsomatic_pon(ponpath,outpath):
    from baseq.snv.gatk import run_createsomatic_pon
    run_createsomatic_pon(ponpath,outpath)


