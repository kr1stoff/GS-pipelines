#!/usr/bin/env perl

my ($fastq) = @ARGV;


my $gc_cnt = 0;
my $q20_cnt = 0;
my $q30_cnt = 0;
my $total_reads = 0;
my $total_bases = 0;


open FASTQ, $fastq or die "can't open $fastq!\n";
while (my $readid = <FASTQ>) {
    my $seq     = <FASTQ>;
    my $comment = <FASTQ>;
    my $quality = <FASTQ>;

    chomp($readid);
    chomp($comment);
    chomp($seq);
    chomp($quality);


    $total_reads += 1;
    $total_bases += length($seq);



    my $cnt = $seq =~ tr/GCgc//;

    $gc_cnt += $cnt;

    foreach my $x (split //, $quality) {
        $q20_cnt += 1 if ord($x) - 33 >= 20;
        $q30_cnt += 1 if ord($x) - 33 >= 30;
    }

}
close FASTQ;

print qq{$total_reads\t$total_bases\t$q20_cnt\t$q30_cnt\t$gc_cnt\n};