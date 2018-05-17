CNV
=====


Config
--------------

.. code-block:: sh

    [CNV]
    bowtie2 = /mnt/gpfs/Database/softs/anaconda2/bin/bowtie2
    samtools = /mnt/gpfs/Database/softs/anaconda2/bin/samtools

    [CNV_ref_hg19]
    bowtie2_index = /mnt/gpfs/Database/ref/hg19/hg19
    dynamic_bin = /mnt/gpfs/Users/zhangxiannian/basematic/cnv/hg19.dynabin.txt

.. code-block:: sh
    baseq-CNV ...

Dynamic Bins
---------------

*Dynamic Bin*: can be downloaded from github_
::
    datas containing columns.

.. _github: https://pypi.org/project/MarkupSafe/