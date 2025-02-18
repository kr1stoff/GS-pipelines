#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use autodie;

my ($bamfile,$outdir,$regionfile,$Plot,$help);

GetOptions(
		"i:s"=>\$bamfile,
		"o:s"=>\$outdir,
		"r:s"=>\$regionfile,
		"h"=>\$help,
		"plot"=>\$Plot,
		);

my $usage=<<USAGE;
usage:perl $0
		-i <bamfile>
		-r [region file]
		-o <outdir>
		-plot 
		-h help
USAGE

die $usage if (!$bamfile || $help || !$outdir);

my $Initial_bases_on_wg=0;
my $Initial_bases_on_target=0;
my $Total_reads=0;
my $Total_yield=0;
my $Total_Q20_bases=0;
my $Total_Q30_bases=0;
my $Average_qual=0;
my $Effective_sequences_on_target=0;
my $Total_effective_yield=0;
my $Average_read_length=0;
my $Mapped_to_genome=0;
my $Average_sequencing_depth=0;
my $Average_sequencing_depth_on_target=0;
my $sambamba="/media/gsadmin/vd1/heguangliang/pipeline/WGS/bin/sambamba";
print "read bam\n";
`mkdir -p $outdir` unless -e $outdir;
open BAM,"$sambamba view -h $bamfile 2>/dev/null| ";
while(<BAM>)
{
	if(m/^\@SQ/){
		chomp;
		$Initial_bases_on_wg += $1 if m/LN:(\d+)/;
	}else{
		next if m/^\@/;
		chomp;
		my @t=split (/\t/,$_);
		$Total_reads++;
		my @qual_s=split "",$t[10];
		$Total_yield += length($t[10]);
		foreach my $i(@qual_s){
			my $q=ord($i)-33;
			$Average_qual += $q;
			$Total_Q20_bases++ if $q >= 20;
			$Total_Q30_bases++ if $q >= 30;
		}
		unless($t[1] & 0x4){
			$Mapped_to_genome++;
			$Total_effective_yield += length($t[9]);
		}
	}
}
close BAM;
$Average_qual /= $Total_yield;
$Average_read_length=$Total_effective_yield/$Mapped_to_genome;
$Average_sequencing_depth=$Total_effective_yield/$Initial_bases_on_wg;
my $name=basename($bamfile);
my $sample=(split /\./,$name)[0];
open STAT,">$outdir/information.xls" or die $!;
print STAT "Sample\t$sample\n";
print STAT "Total reads\t$Total_reads\n";
printf STAT "Total yield(Mb)\t%.2f\n",$Total_yield/1000000;
printf STAT "Average read length(bp)\t%.2f\n",$Average_read_length;
printf STAT "Average base quality\t%.2f\n",$Average_qual;
printf STAT "Fraction of Q30 bases\t%.2f%%\n",100*$Total_Q30_bases/$Total_yield;
printf STAT "Fraction of Q20 bases\t%.2f%%\n",100*$Total_Q20_bases/$Total_yield;
printf STAT "Total mapped reads\t$Mapped_to_genome\n";
printf STAT "Mappping rate\t%.1f%%\n",100*$Mapped_to_genome/$Total_reads;
printf STAT "Initial bases on whole genome(Mb)\t%.2f\n",$Initial_bases_on_wg/1000000;
printf STAT "Total mapped yield(Mb)\t%.2f\n",$Total_effective_yield/1000000;
printf STAT "Fraction of effective bases\t%.1f%%\n",100*$Total_effective_yield/$Total_yield;
printf STAT "Average sequencing depth on whole genome\t%.2f\n",$Average_sequencing_depth;
my %hash=();
print "read region\n";
if($regionfile){
	my $Base_covered_on_target=0;
	`cut -f1,2,3 $regionfile|bedtools sort -i stdin|bedtools merge -i stdin >$outdir/regionfile.bed`;
	open REG,"$outdir/regionfile.bed" or die $!;
	while(<REG>){
		chomp;
		next if(/^$/);
		my @info=split (/\t/,$_);
		$hash{$info[0]}+=$info[2]-$info[1];
		$Initial_bases_on_target+=$info[2]-$info[1];
	}
	close(REG);
	`awk '{print \$1"\t"\$2+1"\t"\$3+1}' $outdir/regionfile.bed >$outdir/regionfile.txt`;
	`$sambamba depth base -t 8 -L $outdir/regionfile.txt $bamfile 2>/dev/null|grep -v REF >$outdir/target.depth`;
	$Effective_sequences_on_target=`awk '{total+=\$3};\$3 != 0 {cov+=1};END{print total"\t"cov}' $outdir/target.depth`;
	chomp($Effective_sequences_on_target);
	($Effective_sequences_on_target,$Base_covered_on_target)=split(/\t/,$Effective_sequences_on_target);
	my $Fraction_of_effective_bases_on_target=$Effective_sequences_on_target/$Total_effective_yield;
	my $Mapped_to_target=`$sambamba view -L $regionfile $bamfile 2>/dev/null | wc -l`;
	chomp($Mapped_to_target);	
	my $Fraction_of_mapping_on_target = $Mapped_to_target/$Mapped_to_genome;
	$Average_sequencing_depth_on_target=$Effective_sequences_on_target/$Initial_bases_on_target;
	my $Coverage_of_target_region=$Base_covered_on_target/$Initial_bases_on_target;

	my $cutoff1=int(0.2*$Average_sequencing_depth_on_target);
	my $cutoff2=int($Average_sequencing_depth_on_target);
	my $tmp1=`awk '\$3 >=20 {total1++};\$3 >=10 {total2++};\$3 >=4 {total3++};\$3 >=$cutoff1 {total4++};\$3 >=$cutoff2 {total5++};END{print total1"\t"total2"\t"total3"\t"total4"\t"total5}' $outdir/target.depth`;
	chomp($tmp1);
	my @info1;
	@info1=split /\t/,$tmp1;
	$info1[0]=0 unless($info1[0] =~ /\d+/);
	$info1[1]=0 unless($info1[1] =~ /\d+/);
	$info1[2]=0 unless($info1[2] =~ /\d+/);
	$info1[3]=0 unless($info1[3] =~ /\d+/);
	$info1[4]=0 unless($info1[4] =~ /\d+/);
	my $Fraction_of_target_covered_with_at_least_20x=$info1[0]/$Initial_bases_on_target;
	my $Fraction_of_target_covered_with_at_least_10x=$info1[1]/$Initial_bases_on_target;
	my $Fraction_of_target_covered_with_at_least_4x=$info1[2]/$Initial_bases_on_target;
	my $Fraction_of_target_covered_with_at_least_02xAverage_sequencing_depth_on_target=$info1[3]/$Initial_bases_on_target;
	my $Fraction_of_target_covered_with_at_least_Average_sequencing_depth_on_target=$info1[4]/$Initial_bases_on_target;

	printf STAT "Total reads mapped on target\t$Mapped_to_target\n";
	printf STAT "Fraction of mapped reads on target\t%.1f%%\n",100*$Fraction_of_mapping_on_target;
	printf STAT "Effective sequences on target(Mb)\t%.2f\n",$Effective_sequences_on_target/1000000;
	printf STAT "Initial bases on target\t$Initial_bases_on_target\n";	
	printf STAT "Fraction of effective bases on target\t%.1f%%\n",100*$Fraction_of_effective_bases_on_target;
	printf STAT "Average sequencing depth on target\t%.2f\n",$Average_sequencing_depth_on_target;
	printf STAT "Base covered on target\t$Base_covered_on_target\n";
	printf STAT "Coverage of target region\t%.1f%%\n",100*$Coverage_of_target_region;
	printf STAT "Fraction of target covered with at least 20x\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_20x;
	printf STAT "Fraction of target covered with at least 10x\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_10x;
	printf STAT "Fraction of target covered with at least 4x\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_4x;
	printf STAT "Uniformity\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_02xAverage_sequencing_depth_on_target;
	printf STAT "Fraction of target covered with at least Average sequencing depth on target\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_Average_sequencing_depth_on_target;
}else{
	`$sambamba depth base -t 8 $bamfile 2>/dev/null|grep -v REF >$outdir/wg.depth`;
	my $cutoff1=int(0.2*$Average_sequencing_depth);
	my $cutoff2=int($Average_sequencing_depth);
	my $tmp1=`awk '\$3 >=20 {total1++};\$3 >=10 {total2++};\$3 >=4 {total3++};\$3 >=$cutoff1 {total4++};\$3 >=$cutoff2 {total5++};END{print total1"\t"total2"\t"total3"\t"total4"\t"total5}' $outdir/wg.depth`;
	chomp($tmp1);
	my @info1;
	@info1=split /\t/,$tmp1;
	$info1[0]=0 unless($info1[0] =~ /\d+/);
	$info1[1]=0 unless($info1[1] =~ /\d+/);
	$info1[2]=0 unless($info1[2] =~ /\d+/);
	$info1[3]=0 unless($info1[3] =~ /\d+/);
	$info1[4]=0 unless($info1[4] =~ /\d+/);
	my $Fraction_of_covered_with_at_least_20x=$info1[0]/$Initial_bases_on_wg;
	my $Fraction_of_covered_with_at_least_10x=$info1[1]/$Initial_bases_on_wg;
	my $Fraction_of_covered_with_at_least_4x=$info1[2]/$Initial_bases_on_wg;
	my $Fraction_of_covered_with_at_least_02xAverage_sequencing_depth=$info1[3]/$Initial_bases_on_wg;
	my $Fraction_of_covered_with_at_least_Average_sequencing_depth=$info1[4]/$Initial_bases_on_wg;
	printf STAT "Fraction of whole genome covered with at least 20x\t%.1f%%\n",100*$Fraction_of_covered_with_at_least_20x;
	printf STAT "Fraction of whole genome covered with at least 10x\t%.1f%%\n",100*$Fraction_of_covered_with_at_least_10x;
	printf STAT "Fraction of whole genome covered with at least 4x\t%.1f%%\n",100*$Fraction_of_covered_with_at_least_4x;
	printf STAT "Uniformity\t%.1f%%\n",100*$Fraction_of_covered_with_at_least_02xAverage_sequencing_depth;
	printf STAT "Fraction of whole genome covered with at least Average sequencing depth\t%.1f%%\n",100*$Fraction_of_covered_with_at_least_Average_sequencing_depth;
}
close STAT;

print "generate plot information\n";
if($regionfile){
	open DB,"$outdir/target.depth";
}else{
	open DB,"$outdir/wg.depth";
}
open DF,">$outdir/depth_frequency.xls";
my %depth=();
while(<DB>)
{
	chomp;
	my @tmp=split;
	$depth{$tmp[2]}++;
}
close(DB);
my $maxCov=0;
	
foreach my $depth (sort {$a<=>$b} keys %depth)
{
	next if($depth==0);
	my $per;
	if($regionfile){$per=$depth{$depth}/$Initial_bases_on_target}else{$per=$depth{$depth}/$Initial_bases_on_wg}
	$maxCov = $per if($per > $maxCov);
	print DF "$depth\t$per\t$depth{$depth}\n";
}
close(DF);
	
open CU,">$outdir/cumu.xls";
print CU "Depth\tTRPercent\n";
my @depth= sort {$a<=>$b} keys %depth;

foreach my $depth1 (sort {$a<=>$b} keys %depth)
{
	my $tmp=0;
	next if($depth1==0);
	foreach my $depth2 (@depth)
	{
		if($depth2 >= $depth1)
		{
			$tmp+=$depth{$depth2};
		}
	}
	if($regionfile){$tmp = $tmp/$Initial_bases_on_target}else{$tmp = $tmp/$Initial_bases_on_wg}
	print CU "$depth1\t$tmp\n";
}
close(CU);

print "draw pictures\n";
if($Plot)
{
	histPlot($outdir,"$outdir/depth_frequency.xls");
	cumuPlot($outdir,"$outdir/cumu.xls");
}

my %dep=();
my %cov=();
if($regionfile){
	open IN,"$outdir/target.depth";
}else{
	open IN,"$outdir/wg.depth";
}
while (<IN>)
{
	chomp;
	my ($chr,$dd)=(split /\t/,$_)[0,2];
	$cov{$chr}++;
	$dep{$chr}+=$dd;
}
close IN;
#if($regionfile){`rm  $outdir/target.depth`}else{`rm $outdir/wg.depth`}

sub cumuPlot {
        my ($outdir, $dataFile) = @_;
        my $figFile = "$outdir/cumuPlot.pdf";
		my $average_depth;
		if($regionfile){
			$average_depth=$Average_sequencing_depth_on_target
		}else{
			$average_depth=$Average_sequencing_depth
		}
        my $Rline=<<Rline;
		library(Cairo)
		rt <- read.table("$dataFile",header=T)
		average_depth <- $Average_sequencing_depth_on_target
		max_lim=2*average_depth
		mi=log10(average_depth)%/%1
		x=rt\$Depth[rt\$Depth<=max_lim]
		y=100*rt\$TRPercent[rt\$Depth<=max_lim]
        CairoPNG(file="$outdir/cumuPlot.png",width = 960, height = 720)
        plot(x,y,col="red",type='l', lwd=3, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
		xpos <- seq(10**mi,max_lim,by=10**mi)
		ypos <- seq(0,100,by=20)
		axis(side=1, xpos, tcl=0.2, labels=FALSE)
		axis(side=2, ypos, tcl=0.2, labels=FALSE)
		mtext("Cumulative sequencing depth",side=1, line=2.5, at=median(xpos), cex=1.5 )
		mtext("Fraction of target bases (%)",side=2, line=2.5, at=median(ypos), cex=1.5 )
		mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1)
		mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1)
		dev.off()

        pdf(file="$figFile",w=8,h=6)
        opar <- par()
		par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="red",type='l', lwd=3, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
		xpos <- seq(10**mi,max_lim,by=10**mi)
		ypos <- seq(0,100,by=20)
		axis(side=1, xpos, tcl=0.2, labels=FALSE)
		axis(side=2, ypos, tcl=0.2, labels=FALSE)
		mtext("Cumulative sequencing depth",side=1, line=2.5, at=median(xpos), cex=1.5 )
		mtext("Fraction of target bases (%)",side=2, line=2.5, at=median(ypos), cex=1.5 )
		mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1)
		mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1)
        par(opar)
        dev.off()
Rline
        open (ROUT,">$figFile.R");
        print ROUT $Rline;
        close(ROUT);

        system("Rscript $figFile.R");
}


sub histPlot {
	my ($outdir, $dataFile) = @_;
	my $figFile = "$outdir/histPlot.pdf";
	my $Rline=<<Rline; 
	library(Cairo)
	rt <- read.table("$dataFile")
	average_depth <- $Average_sequencing_depth_on_target
	max_lim=2
	min_lim=0.1
	rt\$V1=rt\$V1/average_depth
	y=c()
	y[1]=sum(rt\$V2[rt\$V1<0.1])
	y[21]=sum(rt\$V2[rt\$V1>=2])
	for(i in 2:20){y[i]=sum(rt\$V2[rt\$V1>=0.1*(i-1) & rt\$V1<0.1*i])}
	y=y*100
	y_max_lim=0
	for(i in seq(5,50,by=5)){if(i>max(y)){y_max_lim=i;break}}
	y[y==0]="NA"
	x=seq(0,2,by=0.1)
	CairoPNG(file="$outdir/histPlot.png",width = 960, height = 720)
	plot(x,y,col="blue",type='h', lwd=8, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,y_max_lim),xlim=c(0,2))
	xpos <- seq(0,2,by=0.1)
	ypos <- seq(0,y_max_lim,by=1)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth / Average sequencing depth on target",side=1, line=3, at=median(xpos), cex=2 )
	mtext("Fraction of target bases (%)",side=2, line=2, at=median(ypos), cex=2 )
	mtext(c("<0.1",xpos[2:20],">2"), side=1, las=1, at=xpos, line=0.5, cex=1)
	mtext(ypos, side=2, las=1, at=ypos, line=0.5, cex=1)
	dev.off()

	pdf(file="$figFile",w=8,h=6)
	opar <- par()
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=8, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,y_max_lim),xlim=c(0,2))
	xpos <- seq(0,2,by=0.1)
	ypos <- seq(0,y_max_lim,by=1)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth / Average sequencing depth on target",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of target bases (%)",side=2, line=2, at=median(ypos), cex=2 )
	mtext(c("<0.1",xpos[2:20],">2"), side=1, las=1, at=xpos, line=0.5, cex=1)
	mtext(ypos, side=2, las=1, at=ypos, line=0.5, cex=1)
	par(opar)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("Rscript $figFile.R");
}
