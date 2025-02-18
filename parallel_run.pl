#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd 'abs_path';
use threads;

my ($lines, $max) = (1, 12);
GetOptions(
	"l:i" => \$lines,
	"m:i" => \$max,
);

my $usage = qq{
Usage: perl $0 [options] in.sh

Options:
-l	# of lines of each job [1]
-m	# of jobs to run concurrently [8]
};

die $usage unless(@ARGV);

my $dir = abs_path($ARGV[0]) . '.workdir';
`mkdir -p $dir`;

#read input file and split
my $line_count = 0;
my $job_count = '0001';
my $content = '';
my (@script, @running_script, @running_script_name, @logFiles);
while(<>){
	next if(/^\s*$/ || /^\s*#/);
	$line_count++;
	$content .= $_;
	if($line_count % $lines == 0){
		&output_to($content, "$dir/$job_count.sh");
		$job_count++;
		$content = '';
	}
}
&output_to($content, "$dir/$job_count.sh") if($content);

#run split scripts
for(my $i=0; $i<@script; $i++){
	&join_all if(@running_script >= $max);
	my $tmp =  async { return system("bash $script[$i] >$script[$i].log 2>&1"); };
	push @logFiles,"$script[$i].log";
	push @running_script, $tmp;
	push @running_script_name, $script[$i];
}
&join_all if(@running_script);

## sub routines ##
sub output_to{
	my ($content, $dest) = @_;
	push @script, $dest;
	open OUT, ">$dest" or die $!;
	print OUT $content;
	close OUT;
}

sub join_all{
	system("echo -n 'start running @running_script_name at ';date +'%Y-%m-%d %H:%M:%S'");
	for(@running_script){ $_->join() }
	(@running_script, @running_script_name) = ();
	system("echo -n 'finish at ';date +'%Y-%m-%d %H:%M:%S'");
	for my $f(@logFiles){my $str=`cat $f`;print STDERR "###### Log information for $f #####\n",$str,"###### Log information for $f #####\n"}
}
