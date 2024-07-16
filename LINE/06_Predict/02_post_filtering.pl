#!/usr/bin/env perl

use strict;
use warnings;

my $outpath = $ARGV[0];
my $ver = $ARGV[1];
my $sub = $ARGV[2];
my $cont = $ARGV[3];
my $masterpath = $ARGV[4];
my $hg = $ARGV[5];
my $cutoff = $ARGV[6];
if ($cutoff !~ /\d/)
  {
   die "Cutoff $cutoff is not valid\n";
  }
my $in1_file = $outpath.'/'.$sub.'/retro_v'.$ver.'_1/LINE/'.$sub.'.pred.P1.txt';
my $in2_file = $outpath.'/'.$sub.'/retro_v'.$ver.'_1/LINE/'.$sub.'.pred.P2.txt';
my $ou1_file = $outpath.'/'.$sub.'/retro_v'.$ver.'_0/LINE/RS'.$cutoff.'.0strand.bed';
my $ou2_file = $outpath.'/'.$sub.'/retro_v'.$ver.'_1/LINE/RS'.$cutoff.'.1strand.bed';
my ($flag, $count, $mei, $line, $read_pair, %anno, @data, @temp, %neg_strand, %pos_strand);
my ($read, $c1, $c2, $c3, $c4, $i, $read1, $read2, $strand, $up1, $up2, $PE1, $PE2);
my (%cand, %additional, @temp1, @temp2);

(open (IN1, "<$in1_file")) || die "cannot open the in1 file $in1_file\n";
(open (IN2, "<$in2_file")) || die "cannot open the in2 file\n";
(open (OU1, ">$ou1_file")) || die "cannot open the out file\n";
(open (OU2, ">$ou2_file")) || die "cannot open the out file\n";

while ($line = <IN1>)
   {
    chomp($line);
    next if ($line !~ /\w/);
    @data = split("\t", $line);
    if ($data[2] < $cutoff)
       {
        @temp = split("_", $data[0]);
        $anno{$data[0]}{$data[3]."\t".$data[7]} = $data[1]."\t".$data[2]."\t".$data[3]."\t".$data[4]."\t".$data[5]."\t".$data[6]."\t".$data[7]."\t".$data[8]."\t".$data[9]."\t".$data[10]."\t".$data[11]."\t".$data[12];
        if ($data[1] == 0)
           {
            $neg_strand{$data[0]}{$data[3]} = 1;
            $neg_strand{$data[0]}{$data[7]} = 1;
           }
        else
           {
            $pos_strand{$data[0]}{$data[3]} = 1;
            $pos_strand{$data[0]}{$data[7]} = 1;
           }
       }
   }

while ($line = <IN2>)
   {
    chomp($line);
    next if ($line !~ /\w/);
    @data = split("\t", $line);
    if ($data[2] < $cutoff)
       {
        @temp = split("_", $data[0]);
        $anno{$data[0]}{$data[3]."\t".$data[7]} = $data[1]."\t".$data[2]."\t".$data[3]."\t".$data[4]."\t".$data[5]."\t".$data[6]."\t".$data[7]."\t".$data[8]."\t".$data[9]."\t".$data[10]."\t".$data[11]."\t".$data[12];
        if ($data[1] == 0)
           {
            $neg_strand{$data[0]}{$data[3]} = 1;
            $neg_strand{$data[0]}{$data[7]} = 1;
           }
        else
           {
            $pos_strand{$data[0]}{$data[3]} = 1;
            $pos_strand{$data[0]}{$data[7]} = 1;
           }
       }
   }

for $mei(keys %neg_strand)
   {
    @temp = split('_', $mei);
    $count = scalar keys %{$neg_strand{$mei}};
    print OU1 "$temp[0]\t$temp[1]\t$temp[2]\t$count\n";
   }

for $mei(keys %pos_strand)
   {
    @temp = split('_', $mei);
    $count = scalar keys %{$pos_strand{$mei}};
    print OU2 "$temp[0]\t$temp[1]\t$temp[2]\t$count\n";
   }

system("$masterpath/LINE/06_Predict/03_filter_control.sh $outpath $ver $sub $cont $cutoff 1 $masterpath $hg");
system("$masterpath/LINE/06_Predict/03_filter_control.sh $outpath $ver $sub $cont $cutoff 0 $masterpath $hg");

my $in_file1 = $outpath.'/'.$sub.'/retro_v'.$ver.'_0/LINE/noCONT.cutoff'.$cutoff.'.calls';
my $in_file2 = $outpath.'/'.$sub.'/retro_v'.$ver.'_1/LINE/noCONT.cutoff'.$cutoff.'.calls';
my $out_file = $outpath.'/'.$sub.'/retro_v'.$ver.'_1/LINE/cutoff'.$cutoff.'.candidate.calls';
(open (IN1, "<$in_file1")) || die "cannot open the in file 1\n";
(open (IN2, "<$in_file2")) || die "cannot open the in file 2\n";
(open (OUT, ">$out_file")) || die "cannot open the out file\n";
my ($pe_file, $sr_file);

    for ($strand=0; $strand <2; $strand ++)
       {
        $pe_file = $outpath.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/'.$sub.'.pe.LINE.matrix';
        $sr_file = $outpath.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/'.$sub.'.sr.LINE.matrix';
        (open (PE, "<$pe_file")) || die "cannot open the pe1 file\n";
        (open (SR, "<$sr_file")) || die "cannot open the pe1 file\n";
        <PE>;
        while ($line = <PE>)
           {
            chomp($line);
            @data = split("\t", $line);
            $read = $data[3];
            $c1 = ($data[1] < $data[2]) ? $data[1] :$data[2];
            $c2 = ($data[1] < $data[2]) ? $data[2] :$data[1];
            $c3 = ($data[4] < $data[5]) ? $data[4] :$data[5];
            $c4 = ($data[4] < $data[5]) ? $data[5] :$data[4];
            $additional{$read} = $data[12]."\t".$data[25]."\t".$data[26]."\t".$c1."\t".$c2."\t".$c3."\t".$c4;
          }
        <SR>;
        while ($line = <SR>)
          {
           chomp($line);
           @data = split("\t", $line);
           $read = $data[3];
           $c1 = ($data[1] < $data[2]) ? $data[1] :$data[2];
           $c2 = ($data[1] < $data[2]) ? $data[2] :$data[1];
           $c3 = ($data[4] < $data[5]) ? $data[4] :$data[5];
           $c4 = ($data[4] < $data[5]) ? $data[5] :$data[4];
           $additional{$read} = $data[12]."\t".$data[14]."\t".$data[15]."\t".$c1."\t".$c2."\t".$c3."\t".$c4;
          }
      }

print OUT "chr_cord1_cord2\tcount\tstrand\tscore\tread1\tup1\tRF1\tmap1\tACAG1\tTAG1\tc1\tc2\tc3\tc4\tread2\tup2\tRF2\tmap2\tACAG2\tTAG2\td1\td2\td3\td4\tdist5\tID\tnotes\n";
### add a few arbitrary filter on dist5 ###
### if up1==up2, warn false positive if abs(dist5) > 1 ###
### if up1==up2 & PE1 == PE2 ==0, warn false positive if abs(dist5)> 0.02 ###
### check similarity of supporting reads ###
while ($line = <IN1>)
   {
    chomp($line);
    @data = split("\t", $line);
    $mei = $data[0].'_'.$data[1].'_'.$data[2];
    $cand{$mei} = $data[3];
   }

while ($line = <IN2>)
   {
    chomp($line);
    @data = split("\t", $line);
    $mei = $data[0].'_'.$data[1].'_'.$data[2];
    $cand{$mei} = $data[3];
   }

for $mei (keys %cand)
   {
    for $read_pair(keys %{$anno{$mei}})
       {
        $flag = '';
        @temp = split("\t",$anno{$mei}{$read_pair});
        if ($temp[4] == $temp[8])
           {
            $flag = 'Dist5_is_too_large_' if (abs($temp[10]) > 1);
            if ($temp[3] == 0 && $temp[7] == 0)
               {
                $flag = 'split_read_is_too_far_away_' if (abs($temp[10]) > 0.02);
               }
           }
        if (($temp[11] =~ /\d/) && ($temp[11] < 95))
          {
           $flag .= 'ID_too_low';
          }
        $PE1  =  ($temp[3] == 1) ? 'PE' : 'SR';
        $PE2  =  ($temp[7] == 1) ? 'PE' : 'SR';
        if ($temp[0] ==1)
           {
            $strand = '+strand';
            $up1  =  ($temp[4] == 1) ? 'upstream' : 'downstream';
            $up2  =  ($temp[8] == 1) ? 'upstream' : 'downstream'; 
           }
        else
           {
            $strand = '-strand';
            $up1  =  ($temp[4] == 1) ? 'downstream' : 'upstream';
            $up2  =  ($temp[8] == 1) ? 'downstream' : 'upstream';
           }
        $read1 = $temp[2]."\t".$PE1."\t".$up1."\t".$temp[5]."\t".$additional{$temp[2]};
        $read2 = $temp[6]."\t".$PE2."\t".$up2."\t".$temp[9]."\t".$additional{$temp[6]};
        @temp1 = split("\t", $additional{$temp[2]});
        @temp2 = split("\t", $additional{$temp[6]});
        if ( ($temp1[0] < 95) || ($temp2[0] < 95) )
           {
            $flag .= 'lowmap';
           }
        elsif ( ($temp1[1] eq '0') || ($temp2[1] eq '0') )
           {
            $flag .= 'noACAG';
           }
        elsif ( ($temp1[2] eq '0') || ($temp2[2] eq '0') )
           {
            $flag .= 'noTAG';
           }

        if (($up1 eq 'upstream') && ($up2 eq 'downstream'))
          {
           $flag .= 'longTSD' if ($temp1[4] - $temp2[3] > 25);
          }
        elsif (($up1 eq 'downstream') && ($up2 eq 'upstream'))
          {
           $flag .= 'longTSD' if ($temp2[4] - $temp1[3] > 25);
          } 
        if ( ($strand eq '+strand') && ((($up1 eq 'downstream') && (6100-$temp1[6]>400)) || (($up2 eq 'downstream') && (6100-$temp2[6]>400))) )
          {
           $flag .= 'LargeployAGap';
          }
        elsif ( ($strand eq '-strand') && ((($up1 eq 'upstream') && (6100-$temp1[6]>400)) || (($up2 eq 'upstream') && (6100-$temp2[6]>400))) )
          {
           $flag .= 'LargeployAGap';
          }
        $flag = 'ok' if ($flag !~ /\w/);
   
        if ($temp1[3] < $temp2[3])
           {
            print OUT "$mei\t$cand{$mei}\t$strand\t$temp[1]\t$read1\t$read2\t$temp[10]\t$temp[11]\t$flag\n";
           }
        else
           {
            print OUT "$mei\t$cand{$mei}\t$strand\t$temp[1]\t$read2\t$read1\t$temp[10]\t$temp[11]\t$flag\n";
           }
       }
    } 

