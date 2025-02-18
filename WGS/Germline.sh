#!/usr/bin/bash
######################################################
RGSM=GS03615

######################################################
gatk=/media/gsadmin/vd2/tmp/software/GenomeAnalysisTK/gatk-4.0.4.0/gatk

# 数据库
reference=/media/gsadmin/vd2/tmp/database/GATK/library/b37/human_g1k_v37_decoy.fasta
bundle=/media/gsadmin/vd2/tmp/database/GATK/library/b37
mkdir ../Vcf
####
# time $gatk HaplotypeCaller \
# 	-R $reference \
# 	-I ../Bam/${RGSM}.final.bam \
# 	-O ../Vcf/${RGSM}.HC.vcf.gz && \
# 	echo "** ${sample}.HC.vcf.gz done **"

# VQSR校正
time $gatk VariantRecalibrator \
	-R $reference \
	-V ../Vcf/${RGSM}.HC.vcf.gz \
	-resource hapmap,known=false,training=true,truth=true,prior=15.0:$bundle/hapmap_3.3_b37_pop_stratified_af.vcf.gz \
	-resource omni,known=false,training=true,truth=false,prior=12.0:$bundle/1000G_omni2.5.b37.vcf \
	-resource 1000G,known=false,training=true,truth=false,prior=10.0:$bundle/1000G_phase1.snps.high_confidence.b37.vcf \
	-resource dbsnp,known=true,training=false,truth=false,prior=6.0:$bundle/dbsnp_138.b37.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	--rscript-file ../Temp/${RGSM}.HC.snps.plots.R \
	--tranches-file ../Temp/${RGSM}.HC.snps.tranches \
	-O ../Temp/${RGSM}.snp.recal

time $gatk ApplyVQSR \
	-R $reference \
	-V ../Vcf/${RGSM}.HC.vcf.gz \
	-ts-filter-level 99.0 \
	--tranches-file ../Temp/${RGSM}.HC.snps.tranches \
	--recal-file ../Temp/${RGSM}.snp.recal \
	-mode SNP \
	-O ../Temp/${RGSM}.snp.vcf.gz
#--max-gaussians 4用于设定Gaussians（clusters of variants that have similar properties）的数目，即减少聚类的组数，从而使得每个组的变异位点数目达到要求
time $gatk VariantRecalibrator \
	-R $reference \
	-V ../Temp/${RGSM}.snp.vcf.gz \
	-resource mills,known=true,training=true,truth=true,prior=12.0:$bundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode INDEL \
	--max-gaussians 4 \
	-rscript-file ../Temp/${RGSM}.HC.indels.plots.R \
	--tranches-file ../Temp/${RGSM}.HC.indels.tranches \
	-O ../Temp/${RGSM}.indels.recal

time $gatk ApplyVQSR \
	-R $reference \
	-V ../Temp/${RGSM}.snp.vcf.gz \
	-ts-filter-level 99.0 \
	--tranches-file ../Temp/${RGSM}.HC.indels.tranches \
	--recal-file ../Temp/${RGSM}.indels.recal \
	-mode INDEL \
	-O ../Vcf/${RGSM}.VQSR.vcf.gz
