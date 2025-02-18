#!/usr/bin/env perl

my ($gtf) = @ARGV;

my %hash = ();
open GTF, $gtf;
while (<GTF>) {
	chomp;
	my @arr = split /\t/;
	next if $arr[2] ne 'exon';
	my ($trans) = $arr[8] =~ /transcript_name \"(.+?)\"/;
	my $len = $arr[4] - $arr[3];
	$hash{$trans}{$arr[3]} = [$arr[0], $arr[4], $arr[6]]; 
}
close GTF;

foreach my $x (keys %hash) {
	#print qq{$x};
	my ($start, $end) = (0, 0);
	my ($chr, $strand);
	my $cnt = 0;
 	my @exon_start = ();
	my @exon_lens = ();
	foreach my $y (sort {$a <=> $b} keys %{$hash{$x}}) {
		$cnt++;
		$start = $y if $cnt == 1;
		$end   = $hash{$x}{$y}->[1] if  $hash{$x}{$y}->[1] > $end;
		$chr   = $hash{$x}{$y}->[0];
		$strand = $hash{$x}{$y}->[2];
		my $len =  $hash{$x}{$y}->[1] - $y;
		push @exon_lens, $len;
		push @exon_start, $y - $start;
	}
	my $lens  = join ",", @exon_lens;
	my  $init = join ",", @exon_start;
	print qq{$chr\t$start\t$end\t$x\t0\t$strand\t$start\t$end\t0\t$cnt\t$lens\t$init\n};
}
