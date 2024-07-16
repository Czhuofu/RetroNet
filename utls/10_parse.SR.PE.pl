#! /usr/bin/perl -w

use strict;

my $in_file = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[3].'/'.$ARGV[0].'.'.$ARGV[3].'.SR.PE.calls.test';
my $ou_file = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[3].'/'.$ARGV[0].'.'.$ARGV[3].'.SR.PE.calls';

my ($line, @data, @temp, @single);
my ($total, $pe, $sr, $support);
(open (IN, "<$in_file")) || die "cannot open the in file\n";
(open (OU, ">$ou_file")) || die "cannot open the ou file\n";

while ($line = <IN>)
   {
    chomp($line);
    @data = split("\t", $line);
    @temp = split(/\;/, $data[3]);
    $sr = $pe = $total = 0;
    foreach $support(@temp)
       {
        next if (!$support);
        @single = split(/\,/, $support);
        if ($single[0] eq 'pe')
           {
            $pe += $single[1];
           }
        else
           {
            $sr += $single[1];
           }
       }
   $total = $sr + $pe;
   print OU "$data[0]\t$data[1]\t$data[2]\t$ARGV[3]\t$total\t$sr\t$pe\t$data[3]\n";
  }
close(IN);

