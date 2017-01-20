#!/bin/bash

#paste <(paste -d '-' <(cat S1_I1.snp.vcf|grep -v '#'| cut -f1-2,4-5) |tr -s "\t" "-") <(cat S1_I1.snp.vcf|grep -v '#'| cut -f10 |awk -F ':' '{print $7"\t"$6}') > S1_I1.snp.txt
#paste <(paste -d '-' <(cat S1_I2.snp.vcf|grep -v '#'| cut -f1-2,4-5) |tr -s "\t" "-") <(cat S1_I2.snp.vcf|grep -v '#'| cut -f10 |awk -F ':' '{print $7"\t"$6}') > S1_I2.snp.txt
#paste <(paste -d '-' <(cat S1_I3.snp.vcf|grep -v '#'| cut -f1-2,4-5) |tr -s "\t" "-") <(cat S1_I3.snp.vcf|grep -v '#'| cut -f10 |awk -F ':' '{print $7"\t"$6}') > S1_I3.snp.txt
#paste <(paste -d '-' <(cat S1_I1.indel.vcf|grep -v '#'| cut -f1-2,4-5) |tr -s "\t" "-") <(cat S1_I1.indel.vcf|grep -v '#'| cut -f10 |awk -F ':' '{print $7"\t"$6}') > S1_I1.indel.txt
#paste <(paste -d '-' <(cat S1_I2.indel.vcf|grep -v '#'| cut -f1-2,4-5) |tr -s "\t" "-") <(cat S1_I2.indel.vcf|grep -v '#'| cut -f10 |awk -F ':' '{print $7"\t"$6}') > S1_I2.indel.txt
#paste <(paste -d '-' <(cat S1_I3.indel.vcf|grep -v '#'| cut -f1-2,4-5) |tr -s "\t" "-") <(cat S1_I3.indel.vcf|grep -v '#'| cut -f10 |awk -F ':' '{print $7"\t"$6}') > S1_I3.indel.txt

#vcf=$*
cd ./VarScan
for n in $(ls *.vcf)
do
	echo $n
	paste <(paste -d '-' <(cat $n |grep -v '#'| cut -f1-2,4-5) |tr -s "\t" "-") <(cat $n | grep -v '#'| cut -f10 | awk -F ':' '{print $7"\t"$6}') > $n.txt
done






