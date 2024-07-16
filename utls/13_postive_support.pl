#! /usr/bin/perl -w

use strict;

my $sub = $ARGV[0];
my $path = $ARGV[1];
my $retro = $ARGV[2];
my $TE = 'ALU';
my $sr_file = $path.'/'.$retro.'/'.$TE.'/'.$sub.'.sr2.pred.filter1.summary';
my $pe_file = $path.'/'.$retro.'/'.$TE.'/'.$sub.'.pe2.pred.filter1.summary';
my $sr_all  = $path.'/'.$retro.'/'.$sub.'.sr.dedup.discover';
my $pe_all  = $path.'/'.$retro.'/'.$TE.'/'.$sub.'.'.$TE.'.nodup.sites';
my $sr_out  = $path.'/'.$retro.'/'.$TE.'/'.$sub.'.sr.pred.posreads';
my $pe_out  = $path.'/'.$retro.'/'.$TE.'/'.$sub.'.pe.pred.posreads';

my ($r1, $line, @data);
my (%sr, %pe);

(open (SRF, "<$sr_file")) || die "cannot open the sr file $sr_file\n";
(open (PEF, "<$pe_file")) || die "cannot open the pe file $pe_file\n";
(open (SRA, "<$sr_all ")) || die "cannot open the sr all \n";
(open (PEA, "<$pe_all ")) || die "cannot open the pe all \n";
(open (SRO, ">$sr_out ")) || die "cannot open the sr all \n";
(open (PEO, ">$pe_out ")) || die "cannot open the sr all \n";

while ($line = <SRF>)
   {
    chomp($line);
    @data = split("\t", $line);
    if ($data[8] eq '1')
       {
        $sr{$data[3]} = 1;
       }
   }
while ($line = <PEF>)
   {
    chomp($line);
    @data = split("\t", $line);
    if ($data[9] eq '1')
       {
        $pe{$data[3]} = 1;
       }
   }

while ($line = <SRA>)
   {
    chomp($line);
    @data = split("\t", $line);
    if ( ($data[3] =~ /Alu/) && (exists($sr{$data[4]})) )
       {
        print SRO "$line\n";
       }
   }

while ($line = <PEA>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    if ( ($data[3] =~ /Alu/) && (exists($pe{$r1.':'.$data[4]})) )
       {
        print PEO "$line\n";
       }
   }

