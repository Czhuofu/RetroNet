#! /usr/bin/perl -w

use strict;

my $sample = $ARGV[0];
my $path   = $ARGV[1];
my $retro  = $ARGV[2];
my $TE     = $ARGV[3];

my $call_file = $path.'/'.$retro.'/'.$TE.'/'.$sample.'.'.$TE.'.PE.nodup.calls';
my $read_file = $path.'/'.$retro.'/'.$TE.'/'.$sample.'.'.$TE.'.novel.sites';
my $loc_file  = $path.'/'.$retro.'/'.$TE.'/'.$sample.'.'.$TE.'.PE.refine.calls';

my ($line, @data);
my ($chr, $start, $end, $te, $count, $c1, $c2, $c3, $c4, $read, $sign, $call, $cord);
my (%signs, @direction);
my $OVER = 50;
my $TSD  = 25;
(open (READ, "<$read_file")) || die "cannot open the read file $read_file\n";
(open (CALL, "<$call_file")) || die "cannot open the call file $call_file\n";
(open (LOC, ">$loc_file")) || die "cannot open the out file $loc_file\n";

#11      38699654        38699754        Alu     HWI-ST1158:42:C0E1TACXX:8:2202:9860:103910      +       -       89.74
$count = 0;
while ($line = <READ>)
   {
    chomp($line);
    @data = split("\t", $line);
    $count ++;
    if ($data[5] eq '+')
       {
        $direction[$count] = $data[2] - $TSD;
       }
    else
       {
        $direction[$count] = -1*$data[1] - $TSD;
       }
   }

while ($line = <CALL>)
   {
    chomp($line);
    ($chr, $start, $end, $te, $count, @data) = split("\t", $line);
    next if ($count < 1);
    $c1 = $start;
    $c2 = $end;
    %signs = ();
    foreach $read(@data)
        {
         $cord = abs($direction[$read]);
         $sign = ($direction[$read] > 0) ? '+' : '-';
         $signs{$sign} = 1;
         if ($sign eq '+')
            {
             $start = ($start < $cord) ? $cord : $start;
            }
         else
            {
             $end = ($end < $cord) ? $end : $cord;
            }
        }
    $c3 = ($start <= $end) ? $start : $end;
    $c4 = ($start > $end) ? $start : $end;
    if ($start <= ($OVER+$end))
       {
        print LOC "$chr\t$c1\t$c2\t$c3\t$c4\tOK";
       }
    else
       {
        print LOC "$chr\t$c1\t$c2\t$c1\t$c2\tOverlap";
       }
    if (exists($signs{'+'}) && exists($signs{'-'}))
       {
        print LOC "\t+\/-";
       }
    elsif (exists($signs{'+'}))
       {
        print LOC "\t+";
       }
    else
       {
        print LOC "\t-";
       }
    print LOC "\t$te\t$count";
    foreach $read(@data)
       {
        $sign = ($direction[$read] > 0) ? '+' : '-';
        print LOC "\t$read\_$sign";
       }
    print LOC "\n";
   }
