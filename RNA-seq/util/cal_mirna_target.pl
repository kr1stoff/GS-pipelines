my ($db, $id) = @ARGV;

my %meta = ();

my $cnt = 0;
open EXO, $id or die "Can't open $id!\n";
while (<EXO>) {
        chomp;
        $meta{$_} = "-";
}
close EXO;

my $line_cnt = 0;
open EXO, $db or die "Can't open $db!\n";
while (<EXO>) {
        # body...
        chomp;
        $line_cnt++;
        print qq{$_\n} if $line_cnt == 1;
        my @arr = split /\t/;
	#print qq{$arr[0]\n};
        next if not exists $meta{$arr[0]};
        print qq{$_\n};
}
close EXO;
