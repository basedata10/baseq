import click, os, sys

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)

def cli():
    click.echo("baseq Fastq...")

#sampling from a fastq file ...
@cli.command(short_help="Print Docs")
def doc():
    from .docs import print_doc
    print_doc()

@cli.command(short_help="Quality Control")
@click.option('--samplefile', '-m', default='', help='sample file lists')
@click.option('--fq1', '-1', default='', help='fastq1')
@click.option('--fq2', '-2', default='', help='fastq2')
@click.option('--out', '-o', default='./qc.txt', help='file path (./qc.txt)')
def QC(samplefile, fq1, fq2, out):
    print("[info] Qualtity Static of the Fastq Files ...")
    print("[info] The result file will be write to {}".format(out))
    from baseq.fastq.quality import fastq_basecontent_quality
    from .sample_file import check_sample_files
    result = []
    samples = check_sample_files(samplefile, "sample", fq1, fq2)
    print(samples)

    import xlsxwriter
    workbook = xlsxwriter.Workbook('QC.xlsx')
    workbook.formats[0].set_font_size(12)
    workbook.formats[0].set_font_name('arial')
    format_main = workbook.add_format({'bold': False, 'font_size': 12, 'font_name': 'arial'})
    format_header = workbook.add_format({'bold': True, 'font_size': 15, 'font_name': 'arial'})
    #prepare Page...
    qcpage = workbook.add_worksheet("Report")
    qcpage.set_column('D:D', 40)
    qcpage.set_column('E:E', 40)
    qcpage.write('A1', 'Sample', format_header)
    qcpage.write('B1', 'MeanQuality', format_header)
    qcpage.write('C1', 'BiasIndex', format_header)
    qcpage.write('D1', 'BasePlot', format_header)
    qcpage.write('E1', 'QualityPlot', format_header)

    #build the Excel...
    for idx, sample in enumerate(samples):
        print(idx, sample)
        result = fastq_basecontent_quality(sample[0], sample[1])
        qcpage.set_row(idx+1, 120)
        qcpage.write(idx+1, 0, sample[0], format_main)
        qcpage.write(idx+1, 1, result[2], format_main)
        qcpage.write(idx+1, 2, result[3], format_main)
        qcpage.insert_image(idx+1, 3, result[0], {"x_scale":0.7, "y_scale":0.7, 'x_offset': 5, 'y_offset': 5})
        qcpage.insert_image(idx+1, 4, result[1], {"x_scale":0.7, "y_scale":0.7, 'x_offset': 5, 'y_offset': 5})

    workbook.close()

#sampling from a fastq file ...
@cli.command(short_help="Quality Control")
@click.option('--path', '-d', default='./', help='File Path')
@click.option('--ratio', '-d', default='1', help='ratio of reads')
@click.option('--reads', '-d', default='1', help='read counts')
@click.option('--tofasta', '-d', default='1', help='read counts')
@click.option('--path', '-d', default='./', help='File Path')
def sampling(path):
    from baseq.fastq.quality import fastq_basecontent_quality
    print(fastq_basecontent_quality(path))

@cli.command(short_help="List all the Fastq Files in Path (not include the subdir)")
@click.argument("path")
@click.option('--subdir', '-s', default='', help='The files are in the subdirectory')
@click.option('--save', '-f', default='./samples.txt', help='Write sample infos to this path (./samples.fqs.txt)')
def list_samples(path, save, subdir):
    print("[info] List fastq files in {}".format(path))
    from .sample_file import list_fastq_files
    list_fastq_files(path, save)

@cli.command(short_help="Split the fastq into diffrent barcodes")
@click.argument("barcode")
@click.argument("fastq")
@click.argument("outprefix")
@click.option('--suffix', '-s', default='.fq', help='')
def split_barcode(barcode, fastq, outprefix, suffix):
    print("[info] Split the barcodes...")
    from .split_barcode import split_barcode
    split_barcode(barcode, fastq, outprefix, suffix)

@cli.command(short_help="Split the fastq into diffrent barcodes")
@click.option('--samplefile', '-m', default='', help='sample file lists')
@click.option('--thread', '-t', default=8, help='sample file lists')
@click.option('--name', '-s', default='sample', help='')
@click.option('--seqfile', default='', help='')
@click.option('--fq1', '-1', default='', help='fastq1')
@click.option('--fq2', '-2', default='', help='fastq2')
def filter_polyAT(samplefile, seqfile, fq1, fq2, name, thread):
    print("[info] Filter the Reads with polyA/polyT...")
    from .filter_reads import filter_fastq_pair_by_sequence
    from baseq.fastq.sample_file import check_sample_files
    samples = check_sample_files(samplefile, fq1, fq2)
    from concurrent.futures import ThreadPoolExecutor
    pool = ThreadPoolExecutor(int(thread))
    print("[info] Using the Multiple Threads: {}".format(thread))
    for sample in samples:
        pool.submit(filter_fastq_pair_by_sequence, sample[1], sample[2], seqfile, sample[0])