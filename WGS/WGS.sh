#!/usr/bin/bash

# 执行命令：
# sh WGS.sh {read_1.fq.gz} {read_2.fq.gz} {RGID} {lib} {样本ID} {输出文件夹}

# 软件
trimmomatic=/media/netdisk246/Bioinformatics/software/Trimmomatic-0.36/trimmomatic-0.36.jar
bwa=/media/netdisk246/Bioinformatics/software/bwa_0.7.12/bwa
#samtools=/media/netdisk246/Bioinformatics/software/samtools/samtools
gatk=/media/netdisk246/Bioinformatics/software/gatk-4.0.4.0/gatk

# 数据库
reference=/media/netdisk246/Bioinformatics/database/GATK/library/b37/human_g1k_v37_decoy.fasta
bundle=/media/netdisk246/Bioinformatics/database/GATK/library/b37

# 对数据库建立索引
# $gatk IndexFeatureFile --feature-file $bundle/xxx.vcf

# shell执行参数
fq1=../Rawdata/NKHS180119045-1A_1.fq.gz
fq2=../Rawdata/NKHS180119045-1A_2.fq.gz
RGID=NKHS180119045-1A #一般用Lane ID代替
library=NKHS180119045 #测序文库编号
sample=linqiao #样本ID
outdir=../output #输出目录的路径

# 设置目录
outdir=${outdir}/${sample}

# 获得fq前缀名
# 假设fq文件名格式为*.1.fq.gz
fq_file_name=`basename $fq1`
fq_file_name=${fq1_file_name%%.1.fq.gz}

# 新建文件夹
if [ ! -d $outdir/cleanfq ]
then mkdir -p $outdir/cleanfq
fi

if [ ! -d $outdir/bwa ]
then mkdir -p $outdir/bwa
fi

if [ ! -d $outdir/gatk ]
then mkdir -p $outdir/gatk
fi

# 使用Trimmomatic对原始数据进行质控
time java -jar ${trimmomatic} PE \
	$fq1 $fq2 \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz $outdir/cleanfq/${fq_file_name}.unpaired.1.fq.gz \
	$outdir/cleanfq/${fq_file_name}.paired.2.fq.gz $outdir/cleanfq/${fq_file_name}.unpaired.2.fq.gz \
	ILLUMINACLIP:/media/netdisk246/Bioinformatics/software/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:True \
	SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50 && echo "** fq QC done **"

# 使用bwa mem进行比对
#bwa index $reference
time $bwa mem -t 8 -M -Y -R "@RG\tID:$RGID\tPL:ILLUMINA\tPU:$PU\tLB:$library\tSM:$sample" $reference \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz \
	$outdir/cleanfq/${fq_file_name}.paired.2.fq.gz \
	| samtools view -Sb - > $outdir/bwa/${sample}.bam && \
	echo "** BWA MEM done **"

# samtools排序
#time samtools sort -@ 4 -m 4G -O bam -o $outdir/bwa/${sample}.sorted.bam && \
time samtools sort -@ 4 -m 4G $outdir/bwa/${sample}.bam $outdir/bwa/${sample}.sorted && \
	echo "** sorted raw bam file done **"

# time samtools index $outdir/bwa/${sample}.sorted.bam && \
# 	echo "** sorted bam index done **"

# 标记重复序列###gatk这一步开始出现问题
gatk MarkDuplicates \
	-I $outdir/bwa/${sample}.sorted.bam \
	-O $outdir/bwa/${sample}.sorted.marked.bam \
	-M $outdir/bwa/${sample}.marked_metrics.txt \
	--REMOVE_DUPLICATES false && \
	echo "** mark duplicates done **"

time samtools index $outdir/bwa/$outdir/bwa/${sample}.sorted.markdup.bam && \
	echo "** mark duplicates bam index done **"

# BQSR
time gatk BaseRecalibrator \
	-R $reference \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	--known-sites $bundle/1000G_phase1.indels.b37.vcf \
	--konwn-sites $bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
	--konwn-sites $bundle/dbsnp_150.b37.vcf \
	-O $outdir/bwa/${sample}.sorted.md.recal.table && \
	echo "** ${sample} BQSR done"

time gatk ApplyBQSR \
	--bqsr-recal-fle $outdir/bwa/${sample}.sorted.md.recal.table \
	-R $reference \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	-O $outdir/bwa/${sample}.BQSR.bam && \
	echo "** Apply BQSR done **"

time samtools index $outdir/bwa/${sample}.BQSR.bam && \
	echo "** BQSR index done **"

# 输出vcf
time gatk HaplotypeCaller \
	-R $reference \
	-I $outdir/bwa/${sample}.BQSR.bam \
	-O $outdir/gatk/${sample}.HC.vcf.gz && \
	echo "** ${sample}.HC.vcf.gz done **"

# VQSR校正,首先是SNP
time gatk VariantRecalibrator \
	-R $reference \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	-resourse:hapmap,known=false,training=true,truth=true,prior=15.0 $bundle/hapmap_3.3.b37.vcf \
	-resourse:omini,known=false,training=true,truth=false,prior=12.0 $bundle/omni.vcf \
	-resourse:1000G,known=false,training=true,truth=false,prior=10.0 $bundle/1000G.vcf \
	-resourse:dbsnp_150,known=true,training=false,truth=false,prior=6.0 $bundle/dbsnp_3.3.b37.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	-rscriptFile $outdir/gatk/${sample}.HC.snps.plots.R \
	--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
	-O $oudir/gatk/${sample}.snp.recal

time gatk ApplyVQSR \
	-R $reference \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	--ts_filter_level 99.0 \
	--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
	-recalFile $oudir/gatk/${sample}.snp.recal \
	-mode SNP \
	-O $outdir/gatk/snp.vcf.gz
#然后是indel mode
time gatk VariantRecalibrator \
	-R $reference \
	-input $outdir/gatk/snp.vcf.gz \
	-resourse:mills,known=true,training=true,truth=true,prior=12.0 $bundle/Mills.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode INDEL \
	--max-gaussians 6 \
	-rscriptFile $outdir/gatk/${sample}.HC.indels.plots.R \
	--tranches-file $outdir/${sample}.HC.indels.tranches \
	-O $oudir/gatk/${sample}.snp.indels.recal

time gatk ApplyVQSR \
	-R $reference \
	-input $outdir/gatk/snp.vcf.gz \
	--ts_filter_level 99.0 \
	--tranches-file $outdir/${sample}.HC.snps.tranches \
	-recalFile $oudir/gatk/${sample}.snp.recal \
	-mode INDEL \
	-O $outdir/gatk/VQSR.vcf.gz