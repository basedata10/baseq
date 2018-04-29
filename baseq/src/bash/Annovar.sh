#!/usr/bin/env bash

##annotation
$annovar/convert2annovar.pl -format vcf4 ${sample}_raw_snps.vcf > ${sample}_snps.avinput
$annovar/table_annovar.pl ${sample}_snps.avinput $annovar_db_hg38 -buildver hg38 -out ${sample} -remove \
            -protocol refGene,genomicSuperDups,esp6500siv2_all,ALL.sites.2015_08,EAS.sites.2015_08,exac03,avsnp147,dbnsfp30a,clinvar_20170130,cosmic70,\
            dbscsnv11,cytoBand -operation g,r,f,f,f,f,f,f,f,f,f,r -nastring . -csvout
