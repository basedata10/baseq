CNV
=====
baseq-CNV is a toolkit to infer and visualize copy number from high-throughput DNA
sequencing data. It is designed for use with Whole Genome Sequencing datas for
both bulk and single cell experiments.

The copy number is based on the reads counts per genomic region.
The region are predefined to exclude and discount the low complexity parts.

- Reads Alignment using Bowtie2
- Bin Counting for unique mapped reads
- Normalize by GC content
- CBS for reducing noise

Result
-------

.. csv-table::
   :header: "Name", "Description"
   :widths: 20, 30

   "sample.bowtie.bam", "Aligned bam file"
   "sample.bin.counts.txt", "The counts of reads for each bin in the dynamic_bin_file"
   "sample.CNV_plot_[size].png", "CNV plot figure for each bin-size"
   "sample.GC.png", "GC content datas"


Pipeline
----------
The total pipeline
::
    baseq-CNV run_pipeline ./Tn5_S1.fq.gz -g hg19

Alignment
~~~~~~~~~~~
We use bwa to alignment.
::
    baseq-CNV align -1 Tn5_S1.fq.gz -r 4000000 -g hg19 -t 10

BinCounting
~~~~~~~~~~~~~
According to dynamicbin ... The command is..
::
    baseq-CNV bincount -g hg19 -i ./sample.bam -o normbincounts.txt

Normalize
~~~~~~~~~
Normalize the raw read counts.
::
    baseq-CNV normalize -g hg19 -i ./bincounts.txt -o bincounts_norm.txt

CBS
~~~~~
Segmentation
::
    baseq-CNV cbs -i ./bincounts_norm.txt  -o ./out.file

.. image:: http://p8v379qr8.bkt.clouddn.com/CNV_normalize.png

Plot
~~~~~
Plot genomic...
::
    baseq-CNV plotgenome -i ./bincounts_norm.txt -c ./out.file

Config
-------

.. code-block:: sh

    [CNV]
    bowtie2 = /mnt/gpfs/Database/softs/anaconda2/bin/bowtie2
    samtools = /mnt/gpfs/Database/softs/anaconda2/bin/samtools

    [CNV_ref_hg19]
    bowtie2_index = /mnt/gpfs/Database/ref/hg19/hg19
    dynamic_bin = /mnt/gpfs/Users/zhangxiannian/basematic/cnv/hg19.dynabin.txt

Quality Control
---------------
Alignment inforamtion and MAD

- Alignment: Total reads, mapping ratio
- MAD : Median Absolute Deviations, indicates the technical noise level of the sample.

Dynamic Bins
---------------

*Dynamic Bin*: can be downloaded from github_
::
    datas containing columns.

APIs
-------
.. automodule:: baseq.cnv
    :members:

bincount
~~~~~~~~~
.. automodule:: baseq.cnv.bincount
    :members:

normalize
~~~~~~~~~~~
.. automodule:: baseq.cnv.normalize
    :members:





.. _github: https://pypi.org/project/MarkupSafe/