#!/usr/bin/bash

trimmomatic=/media/netdisk246/Bioinformatics/software/Trimmomatic-0.36/trimmomatic-0.36.jar
fq1=ND180807091_L2_L8_1.fq.gz
fq2=ND180807091_L2_L8_2.fq.gz
RGID=L2L7L8 #一般用Lane ID代替
library=ND180807091 #测序文库编号
sample=GS03615 #样本ID
fq_file_name=ND180807091


mkdir ../cleanfq

# 使用Trimmomatic对原始数据进行质控
time java -jar ${trimmomatic} PE \
	$fq1 $fq2 \
	../cleanfq/${fq_file_name}.paired.1.fq.gz ../cleanfq/${fq_file_name}.unpaired.1.fq.gz \
	../cleanfq/${fq_file_name}.paired.2.fq.gz ../cleanfq/${fq_file_name}.unpaired.2.fq.gz \
	ILLUMINACLIP:/media/netdisk246/Bioinformatics/software/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:True \
	SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50 && echo "** fq QC done **"