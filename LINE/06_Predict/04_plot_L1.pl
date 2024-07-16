#! /usr/bin/env perl 

use strict;
use warnings;

my $outpath = $ARGV[0];
my $sub = $ARGV[1];
my $ver = $ARGV[2];
my $cutoff = $ARGV[3];
my $masterpath = $ARGV[4];
my $hg = $ARGV[5];

my $cand_file = $outpath.'/'.$sub.'/retro_v'.$ver.'_1/LINE/cutoff'.$cutoff.'.candidate.calls';
my ($line, @temp, @data, %plotted);
(open (CAND, "<$cand_file")) || die "cannot open the candidate file $cand_file\n";

<CAND>;
while ($line = <CAND>)
   {
    chomp($line);
    next if ($line !~ /ok/);
    @data = split("\t", $line);
    next if ( ($data[7] < 0.4) || ($data[18] < 0.4) );
    next if (exists($plotted{$data[0]}));
    $plotted{$data[0]} = 1;
    @temp = split('_', $data[0]); 
    print "$temp[0]\t$temp[1]\t$temp[2]\t$data[2]\n";
    if ($data[2] eq '-strand')
       {
        system("$masterpath/visual/RetroVis.sh -i $sub -t LINE -f L1HS -c $temp[0] -d $temp[1] -e $temp[2] -r $ver -s 0 -p $outpath -m $masterpath -g $hg");
       }
    else
       {
        system("$masterpath/visual/RetroVis.sh -i $sub -t LINE -f L1HS -c $temp[0] -d $temp[1] -e $temp[2] -r $ver -s 1 -p $outpath -m $masterpath -g $hg");
       }
   }

