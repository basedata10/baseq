#!/usr/bin/env bash

Ref=/mnt/gpfs/Database/ref/hg19/hg19
seqfile=$1
Outdir=$2
Cut=$3
sample=$4
scriptdir=$5
mkdir $Outdir

/mnt/gpfs/Database/bin/bowtie2 -p 18 -5 $Cut -x $Ref -U $seqfile > $Outdir/$sample.bowtie2.sam 2>> $Outdir/bowtie2.LOG

perl $scriptdir/uniq2bin.50k.pl $Outdir/$sample.bowtie2.sam $scriptdir

/mnt/gpfs/Database/bin/Rscript $scriptdir/cbs.copy.R $scriptdir $Outdir $sample.bowtie2.sam $Outdir

echo -e "$sample\t$Outdir/CNV/$sample.bowtie2.sam.hg19.50k.k50.varbin.data.copynumber.txt" >$Outdir/CNV/temp.txt

/mnt/gpfs/Database/bin/Rscript $scriptdir/plot.figures.R $scriptdir $Outdir/CNV/temp.txt $Outdir/CNV/$sample.genome_CNV.png

cd $Outdir/CNV

/mnt/gpfs/Database/bin/zip ./$sample.chr_pngs.zip ./*png

#curl "http://192.168.2.11:3000/change_state?sample=$sample&action=status&info=Finished"
#curl "http://192.168.2.11:3000/change_state?sample=$sample&action=image&info=/$sample/CNV/$sample.genome_CNV.png"
#curl "http://192.168.2.11:3000/change_state?sample=$sample&action=result_path&info=/$sample/CNV/$sample.chr_pngs.zip"