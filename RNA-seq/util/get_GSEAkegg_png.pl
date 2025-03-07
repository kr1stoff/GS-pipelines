#!usr/bin/perl
use File::Basename;
use LWP::Simple;
use Parallel::ForkManager;

my ($gene_exp, $enrich_xls) = @ARGV;


my %meta    = ();


extract_regulation();



sub extract_regulation
{
	open EXO, $gene_exp or die "Can't open $gene_exp\n";
	while (<EXO>) {
		chomp;
		my @arr = split /\t/;
		my @val = qw/Up Down/;
		next if not $arr[$#arr] ~~ @val;
		#print qq{$arr[$#arr]\n};
		$meta{$arr[0]} = $arr[$#arr];
	}
	close EXO;
}


open EXO, $enrich_xls;
while (<EXO>) {
	chomp;
	next if /^ID/;
	my @arr = split /\t/;
	next if $arr[5]  >0.05;
	foreach my $x (split /\//,$arr[10]) {
		print qq{$arr[0]\t$x\t$meta{$x}\n};
	} 
}
close EXO;
