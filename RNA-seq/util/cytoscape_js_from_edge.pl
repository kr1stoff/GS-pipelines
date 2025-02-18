#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($edge, $out_dir,  $help);

GetOptions(
    'edge|e=s'    => \$edge,
    'out_dir|o=s' => \$out_dir,
    'help|h!'     => \$help
);

if ($help or not $edge or not $out_dir) {
    usage();
    exit();
}

sub usage
{
my $help =<<EOF;
Usage: perl cytoscape_js.pl  -edge edges.txt -out_dir out_dir
    -edge     -e  edges txt  file [required]
    -out_dir  -o  output dir [required]
    -help     -h  print help message
EOF

print $help;

}

my %nodes_attr = ();
my %edges_attr = ();


download_js();
parse_edges();
write_html();
html_2_png();

sub download_js
{

    system qq{mkdir -p $out_dir/js} if not -d qq{$out_dir/js};
    system qq{wget -q -O $out_dir/js/jquery.min.js     "http://cdn.bootcss.com/jquery/1.11.2/jquery.min.js"};
    system qq{wget -q -O $out_dir/js/cytoscape.min.js  "http://cdn.bootcss.com/cytoscape/2.3.16/cytoscape.min.js"};
}


sub parse_edges
{
    my $cnt = 1;
    my $line_cnt = 0;
    open TXT, $edge or die "Can't open $edge!\n";
    while (<TXT>) {
        chomp;
        $line_cnt++;
        next if $line_cnt == 1;
        my @arr = split /\t/;
        if (not exists $nodes_attr{$arr[0]}) {
            $nodes_attr{$arr[0]} = $cnt;
            $cnt++;
        }
        if (not exists $nodes_attr{$arr[1]}) {
            $nodes_attr{$arr[1]} = $cnt;
            $cnt++;
        }
        my $key = qq{$arr[0]\t$arr[1]};
        $edges_attr{$key} = "";
    }
    close TXT; 
}

sub write_html
{

open SAVE, qq{>$out_dir/cytoscape.html};
my @datas = ();

foreach my $x (keys %nodes_attr) {
    my $id   = $nodes_attr{$x};

    my $line = qq{               {data : {id : '$id', name : '$x'}}};
    push @datas, $line;
}

my $nodes_info = join ",\n", @datas;




my @edge_datas = ();
foreach my $x (keys %edges_attr) {
    my ($source, $target) = split /\t/, $x;
    my $source_id   = $nodes_attr{$source};
    my $target_id   = $nodes_attr{$target};


    my $line = qq{              {data : {source : '$source_id', target : '$target_id'}}};
    push @edge_datas, $line;
}
my $edges_info = join ",\n", @edge_datas;




my $help =<<EOF;
<!DOCTYPE html>
<html>
<head>
    <title>Learning Cytoscape.js</title>
    <style type="text/css">
        /* cytoscape graph */
        #cy {
            height: 1200px;
            width: 1800px;
            background-color: #f9f9f9;
        }
    </style>
    <script src="js/jquery.min.js"></script>
    <script src="js/cytoscape.min.js"></script>
    <script>
        \$(function(){
            cytoscape({
              container: document.getElementById('cy'),
              style: [
                { selector: 'node', 
                  css: { 'content': 'data(name)', 'shape': 'ellipse', 'background-color': 'red'}
                }
                                                        
              ],
              elements: {
                nodes:[\n$nodes_info\n],
                edges:[\n$edges_info\n]
              },
              layout: { name: 'cose'} 
            });
        });
    </script>
</head>
<body>
    <div id="cy"></div>
</body>
</html>  
EOF

print SAVE $help;
close SAVE;
}

sub html_2_png
{
my $txt =<<EOF;
var page = require('webpage').create();
page.open('$out_dir/cytoscape.html', function(status) {
    setTimeout(function(){
        page.render('$out_dir/network.png');
        phantom.exit();
    }, 1000000);
});
EOF

open SAVE, qq{>$out_dir/js/html2png.js};
print SAVE qq{$txt\n};
close SAVE;

system qq{/home/genesky/software/phantomjs/2.1.1/bin/phantomjs $out_dir/js/html2png.js\n};
}
