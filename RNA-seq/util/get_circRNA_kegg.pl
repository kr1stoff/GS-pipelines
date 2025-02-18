#!/usr/bin/env perl

my ($xls) = @ARGV;





open EXO, $xls;
while (<EXO>) {
	chomp;
	next if /^ID/;
	my @arr = split /\t/;
	next if $arr[4]  >0.05;
    my @urls = ();
    foreach my $x (split /\//,$arr[7]) {
        my $col  =  'red';
        my $val  = qq{$x%09$col};
        push @urls, $val;
    }
    my $line = join "+", @urls;
    print qq{"http://www.genome.jp/kegg-bin/show_pathway?$arr[0]/$line"\n};
}
close EXO;

