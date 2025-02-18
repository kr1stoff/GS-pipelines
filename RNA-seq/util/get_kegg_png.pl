#!usr/bin/perl
use File::Basename;
use LWP::Simple;
use Parallel::ForkManager;

my ($gene_exp, $enrich_xls, $db) = @ARGV;


my %meta    = ();


my %hash = ();
open TXT, $db;
while (<TXT>) {
	chomp;
	my @arr = split /\t/;
	$hash{$arr[2]} = $arr[0];
}
close TXT;



open EXO, $gene_exp or die "Can't open $gene_exp\n";
while (<EXO>) {
	chomp;
	my @arr = split /\t/;
	my @val = qw/Up Down/;
	next if not $arr[$#arr] ~~ @val;
	my ($id) =  $arr[0] =~ /(\w+)/;
	#print qq{$id\n};
	$meta{$id} = $arr[$#arr];
}
close EXO;



open EXO, $enrich_xls;
while (<EXO>) {
	chomp;
	next if /^ID/;
	my @arr = split /\t/;
	next if $arr[4]  >0.05;
	my @urls = ();
	foreach my $x (split /\//,$arr[7]) {
		my $reg  = lc($meta{$hash{$x}});
		my $col  = $reg eq 'up' ? 'red' : 'blue';
		my $val  = qq{$x%09$col};
		push @urls, $val;
	}
	my $line = join "+", @urls;
	print qq{"http://www.genome.jp/kegg-bin/show_pathway?$arr[0]/$line"\n};
}
close EXO;
