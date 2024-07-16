#!/usr/bin/env perl
### create he stacking file ###
### individual reads ready for pairing ###
use strict;
use warnings;

my $path = $ARGV[0];
my $sub = $ARGV[1];
my $ver = $ARGV[2];
my $strand = $ARGV[3];
my $maxreads = $ARGV[4];
my ($ou_file, $pe_file, $sr_file, $mei_file, $all_pe, $all_sr);
my ($mei, $line, $group, $pe_num, $sr_num, $total_num, $read, $tag);
my (@data, %pe, %sr, %pe_filter, %sr_filter, @reads, @support, %ape, %asr);
my ($c1, $c2, $c3, $direction, $seq, $up, $r1, $cigar);

$mei_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/LINE/'.$sub.'.LINE.novel.calls';
$ou_file  = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/LINE/'.$sub.'.stacking.txt';
$all_pe   = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/LINE/'.$sub.'.LINE.novel.sites';
$all_sr   = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/'.$sub.'.sr.tabe.discover';
$pe_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/LINE/'.$sub.'.pe.prob3.txt';
$sr_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/LINE/'.$sub.'.sr.prob3.txt';

(open (MEI , "<$mei_file")) || die "cannot open the mei file $mei_file\n";
(open (OU  , ">$ou_file ")) || die "cannot open the ou file $ou_file\n";
(open (APE , "<$all_pe  ")) || die "cannot open the site file $all_pe\n";
(open (ASR , "<$all_sr  ")) || die "cannot open the site file $all_sr\n";
(open (PE  , "<$pe_file ")) || die "cannot open the pe file $pe_file \n";
(open (SR  , "<$sr_file ")) || die "cannot open the SR file $sr_file \n";
print OU "insertion\tstrand\tPE\tread\tprob1\tprob2\tprob3\tupstream\tcord1\tcord2\tcord3\tcord4\tsequence\n";
<PE>;
while ($line = <PE>)
   {
    chomp($line);
    @data = split("\t", $line);
    $pe{$data[3]} = $data[6]."\t".$data[7]."\t".$data[8];
   }
<SR>; 
while ($line = <SR>)
   {
    chomp($line);
    @data = split("\t", $line);
    $sr{$data[3]} = $data[6]."\t".$data[7]."\t".$data[8];
   }

while ($line = <APE>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $c1 = ($data[5] eq '+') ? 1 : -1;
    $c2 = ($data[6] eq '+') ? 1 : -1;
    $c3 = ($data[15] eq '+') ? 1 : -1;
    $direction = -1 * $c1 * $c2 * $c3;
    $up = $c1 * $direction;
    $up = 0 if ($up < 0);
    $seq = substr($data[12], $data[8], $data[9]-$data[8]);
    if ($data[4] =~ /^(lib\d+)\_(.+)/)
       {
        $read = $1.'_'. $r1.':'.$2;
       }
    else
       {
        $read = $r1.':'.$data[4];
       }    
    $ape{$read} = $up."\t".$data[1]."\t".$data[2]."\t".$data[10]."\t".$data[11]."\t".$seq; ### upstream cord1 cord2 cord3 cord4 seq ###
   }

while ($line = <ASR>)
   {
    chomp($line);
    next if ($line !~ /L1/);
    @data = split("\t", $line);
    $cigar = $data[6];
    if ($cigar =~ /(\d+)M(\d+)S/)
       {
        $up = ($strand) ? 1 : 0;
       }
    else
       {
        $up = ($strand) ? 0 : 1;
       }
    $asr{$data[4]} = $up."\t".$data[1]."\t".$data[2]."\t".$data[13]."\t".$data[14]."\t".$data[16]; # lib1_r2:MG01HX08:291:HCMHJALXX:6:2202:10906:11382 #
   }

while ($line = <MEI>)
   {
    ### read one call ###
    chomp($line);
    @data = split("\t", $line);
    $pe_num = $sr_num = $total_num = 0;
    %pe_filter = %sr_filter = ();
    $mei = $data[0].'_'.$data[1].'_'.$data[2];

    ### find all of its supporting reads ###
    @support = split(';', $data[7]);
    foreach $group(@support)
       {
        ($tag, @reads) = split(',', $group);
        if ($tag eq 'pe')
           {
            foreach $read(@reads)
               {
                $pe_filter{$read} = $ape{$read} if (exists($pe{$read}));
               }
           }
        elsif ($tag eq 'sr')
           {
            foreach $read(@reads)
               {
                $sr_filter{$read} = $asr{$read} if (exists($sr{$read}));
               }
           }
       }   
    ### count the remaining support, skip if < 2 ###
    $pe_num = scalar keys %pe_filter;
    $sr_num = scalar keys %sr_filter;
    $total_num = $pe_num + $sr_num;
    next if ( ($total_num < 2) || ($total_num > $maxreads) );
 
    ### print output ###
    ### insertion pos strand PE upstream read1 prob1 prob2 prob3 cord1 cord2 cord3 cord4 sequence ###
    for $read(keys %pe_filter)
       {
        print OU "$mei\t$strand\t1\t$read\t$pe{$read}\t$pe_filter{$read}\n";
       }
    for $read(keys %sr_filter)
       {
        print OU "$mei\t$strand\t0\t$read\t$sr{$read}\t$sr_filter{$read}\n";
       }
   }
