#!/usr/bin/env perl

use strict;
use warnings;

my $path   = $ARGV[0];
my $sub    = $ARGV[1];
my $ver    = $ARGV[2];
my $strand = $ARGV[3];
my $in_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/LINE/'.$sub.'.stacking.txt';
my $ou_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/LINE/'.$sub.'.pairs';
my $ins_file = $path.'/'.$sub.'/insert.'.$sub.'.txt';

(open (IN, "<$in_file")) || die "cannot open the in file\n";
(open (OU, ">$ou_file")) || die "cannot open the ou file\n";
(open (INS, "<$ins_file")) || die "cannot open the insert file\n";

my ($line, @data, $ins, $mad, %count, %support, @read1, @read2);
my ($mei, $line1, $line2, $sign, $c1, $c2, $c3, $c4, $d1, $d2, $d3, $d4, $over, $dist1, $dist2, $dist3, $dist4, $prob1, $prob2, $i, $j, $pA1, $pA2, $up1, $up2, $dist5);
my ($lib1, $lib2, @temp, %ins);

chomp ($line = <INS>);
@data = split("\t", $line);
if ($line !~ /lib/)
   {  
    $ins = $data[0];
    $mad = $data[1];
   }
else
   {
    $ins{$data[2]} = $data[0];
    while ($line = <INS>)
       {
        chomp($line);
        @data = split("\t",  $line);
        $ins{$data[2]} = $data[0];
       }
   }

<IN>;
while ($line = <IN>)
   {
    chomp($line);
    @data = split("\t", $line);
    $mei  = $data[0];
    $count{$mei} = 0 if (!exists($count{$mei}));
    $support{$mei}[$count{$mei}++] = $line;
   }

print OU "insertion\tstrand\tread1\tPE1\tup1\tread1_RF\tread1_LR\tread1_NB\tread2\tPE2\tup2\tread2_RF\tread2_LR\tread2_NB\tdist1\tdist2\tdist3\tdist4\tid\tpA1\tpA2\tc1\tc2\tc3\tc4\td1\td2\td3\td4\tdist5\n";

for $mei(keys %support)
   {
    for ($i=0; $i<$#{$support{$mei}}; $i++)
       {
        @read1 = split("\t", $support{$mei}[$i]);
        $prob1 = $read1[4];
        $c1 = ($read1[8] < $read1[9]) ? $read1[8] : $read1[9];
        $c2 = ($read1[8] > $read1[9]) ? $read1[8] : $read1[9];
        $c3 = ($read1[10]<$read1[11]) ? $read1[10] : $read1[11];
        $c4 = ($read1[10]>$read1[11]) ? $read1[10] : $read1[11];
        $pA1 = ($c4 > 6021) ? (($c4-6021)/($c4-$c3)) : 0;
        if ($read1[3] =~ /^(lib\d+)\_/)
           {
            $lib1 = $1; 
           }
        for ($j=$i+1; $j<=$#{$support{$mei}}; $j++)
           {
            @read2 = split("\t", $support{$mei}[$j]);
            $prob2 = $read2[5]; 
            $d1 = ($read2[8] < $read2[9]) ? $read2[8] : $read2[9];
            $d2 = ($read2[8] > $read2[9]) ? $read2[8] : $read2[9];
            $d3 = ($read2[10]<$read2[11]) ? $read2[10] : $read2[11];
            $d4 = ($read2[10]>$read2[11]) ? $read2[10] : $read2[11];
            $pA2 = ($d4 > 6021) ? (($d4-6021)/($d4-$d3)) : 0;

            if ($read2[3] =~ /^(lib\d+)\_/)
               {
                $lib2 = $1;
               
                if ($lib1 =~ /lib/)
                  {
                   $ins = ($ins{$lib1} + $ins{$lib2}) / 2;
                  }
               }
            $dist1 = abs($c1-$d1)/ $ins;
            $dist2 = abs($c2-$d2)/ $ins;
            $dist3 = abs($c3-$d3)/ $ins;
            $dist4 = abs($c4-$d4)/ $ins;
            $over  = &over_id($c3, $c4, $read1[12], $d3, $d4, $read2[12]);

            $up1 = $read1[7];
            $up2 = $read2[7];
            if ( ($strand == 1) && ($up1 == 1) && ($up2 == 1) )
              {
               $sign = ($c2<$d2) ? 1 : -1;
               $dist5 = (($d3 - $c3) - ($d2 - $c2)) * $sign / $ins;
              }
            elsif (($strand == 1) && ($up1 == 0) && ($up2 == 0))
              {
               $sign = ($c1<$d1) ? 1 : -1;
               $dist5 = (($d4 - $c4) - ($d1 - $c1)) * $sign / $ins;
              }
            elsif (($strand == 1) && ($up1 != $up2))
               {
                $sign = ($c2<$d2) ? 1 : -1;
                $dist5 = ($d3 - $c3) * $sign;
               }
            elsif ( ($strand == 0) && ($up1 == 0) && ($up2 == 0) )
               {
                $sign = ($c2<$d2) ? 1 : -1;
                $dist5 = (-1*($d4 - $c4) - ($d2 - $c2)) * $sign / $ins;
               }
            elsif (($strand == 0) && ($up1 == 1) && ($up2 == 1))
               {
                $sign = ($c1<$d1) ? 1 : -1;
                $dist5 = (-1*($d3 - $c3) - ($d1 - $c1)) * $sign / $ins;
               }
            elsif (($strand == 0) && ($up1 != $up2))
               {
                $sign = ($c2<$d2) ? 1 : -1;
                $dist5 = -1*($d3 - $c3) * $sign;
               }

           if ($prob1 > $prob2)
              {
               $line1 = $read1[3]."\t".$read1[2]."\t".$read1[7]."\t".$read1[4]."\t".$read1[5]."\t".$read1[6];
               $line2 = $read2[3]."\t".$read2[2]."\t".$read2[7]."\t".$read2[4]."\t".$read2[5]."\t".$read2[6];
              }
           else
             {
              $line2 = $read1[3]."\t".$read1[2]."\t".$read1[7]."\t".$read1[4]."\t".$read1[5]."\t".$read1[6];
              $line1 = $read2[3]."\t".$read2[2]."\t".$read2[7]."\t".$read2[4]."\t".$read2[5]."\t".$read2[6];
             }
           print OU "$read1[0]\t$read1[1]\t$line1\t$line2\t$dist1\t$dist2\t$dist3\t$dist4\t$over\t$pA1\t$pA2\t$c1\t$c2\t$c3\t$c4\t$d1\t$d2\t$d3\t$d4\t$dist5\n";
         }
      }
   }

sub over_id
  {
   my $x1 = shift;
   my $x2 = shift;
   my $seq1 = shift;
   my $x3 = shift;
   my $x4 = shift;
   my $seq2 = shift;
   my ($temp_file1, $temp_file2, $algn_file);
   my ($y, $line, $align, $total, $match, $over);
   my $tmpfolder=$path.'/'.$sub.'/temp';
   ### need to consider polyA in overlapping sequences ###
   $y = 0;
   if ( ($x1 <= $x3) && ($x2 >= $x3) && ($x4 >= $x2) )
      {
       $y = $x2 - $x3;
      }
   elsif ( ($x1 >= $x3) && ($x1 <= $x4) && ($x2 >= $x4) )
      {
       $y = $x4 - $x1;
      }
   elsif ( ($x3 >= $x1) && ($x4 <= $x2) )
      {
       $y = $x4 - $x3;
      }
   elsif ( ($x1 >= $x3) && ($x4 >= $x2) )
      {
       $y = $x2 - $x1;
      }
   $over = '';
   if ($y > 19)
      {
       $temp_file1 = $path.'/'.$sub.'/temp/'.$sub.'.'.$strand.'.read1.temp.fa';
       $temp_file2 = $path.'/'.$sub.'/temp/'.$sub.'.'.$strand.'.read2.temp.fa';
       $algn_file  = $path.'/'.$sub.'/temp/'.$sub.'.'.$strand.'.reads.temp.alignment';
       (open (TEMP1, ">$temp_file1")) || die "cannot open the temp file\n";
       print TEMP1 ">seq1\n$seq1\n";
       (open (TEMP2, ">$temp_file2")) || die "cannot open the temp file\n";
       print TEMP2 ">seq2\n$seq2\n";
       close(TEMP1);
       close(TEMP2);
       #system("exonerate --bestn 1 --model affine:local --ryo \"INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\\n\" $temp_file1 $temp_file2 | grep INFO > $algn_file");
       system("TMPDIR=$tmpfolder exonerate --bestn 1 --model affine:local --ryo \"INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\\n\" $temp_file1 $temp_file2 | grep INFO > $algn_file");
       (open (ALGN, "<$algn_file")) || die "cannot open the align file\n";
       <ALGN>;
       $line = <ALGN>;
       if ($line)
          {
           @data = split(" ", $line);
           $over = $data[3];
          }
       else
          {
           $over = 0;
          }
      }
   return ($over);
  }

