##################
RGSM=GS03615
RGID=L2L7L8
RGPU=ND180807091
RGLB=ND180807091
sample1=../cleanfq/ND180807091.paired.1.fq.gz
sample2=../cleanfq/ND180807091.paired.2.fq.gz
#################
bwa=/media/netdisk246/Bioinformatics/software/bwa-0.7.17/bwa
gatk=/media/netdisk246/Bioinformatics/software/GenomeAnalysisTK/gatk-4.0.4.0/gatk
# dataset
reference=/media/netdisk246/Bioinformatics/database/GATK/library/b37/human_g1k_v37_decoy.fasta
indel1=/media/netdisk246/Bioinformatics/database/GATK/library/b37/1000G_omni2.5.b37.vcf
indel2=/media/netdisk246/Bioinformatics/database/GATK/library/b37/1000G_phase1.indels.b37.vcf

mkdir ../Temp
mkdir ../Bam

#bwa index $reference

# mapping
bwa mem -t 8 $reference $sample1 $sample2 > ../Temp/${RGSM}.sam

# sorting
samtools view -Sb ../Temp/${RGSM}.sam > ../Temp/${RGSM}.bam
samtools sort -@ 8 ../Temp/${RGSM}.bam ../Temp/${RGSM}.sorted

# mark duplicates N /

gatk MarkDuplicates \
	-I ../Temp/${RGSM}.sorted.bam \
	-O ../Temp/${RGSM}.marked_dup.bam \
	-M ../Bam/${RGSM}_marked_dup.txt \
	--REMOVE_DUPLICATES false

# add head info
gatk AddOrReplaceReadGroups \
	-I ../Temp/${RGSM}.marked_dup.bam \
	-O ../Temp/${RGSM}.addhead.bam \
	--RGID $RGID \
	--RGLB $RGLB \
	--RGPL illumina \
	--RGPU $RGPU \
	--RGSM $RGSM

# make index
gatk BuildBamIndex \
	-I ../Temp/${RGSM}.addhead.bam


# indel region realign
gatk BaseRecalibrator \
	-I ../Temp/${RGSM}.addhead.bam \
	-R $reference \
	--known-sites $indel1 \
	--known-sites $indel2 \
	-O ../Temp/${RGSM}_recal_data.table

# commit 
gatk ApplyBQSR \
	-R $reference \
	-I ../Temp/${RGSM}.addhead.bam \
	--bqsr-recal-file ../Temp/${RGSM}_recal_data.table \
	-O ../Bam/${RGSM}.final.bam

echo "**task done**"
