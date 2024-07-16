#!/usr/bin/env perl

use strict;
use warnings;
### add polyA ###
### add dist5 ###
my $path   = $ARGV[0];
my $sub    = $ARGV[1];
my $ver    = $ARGV[2];
my $strand = $ARGV[3];
my $in_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/ALU/'.$sub.'.stacking.txt';
my $read1_file = '';
my $read2_file = '';
my $align_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/ALU/'.$sub.'.alignment.txt';
system("rm $align_file");
my $ou_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/ALU/'.$sub.'.pairs';
my $ins_file = $path.'/'.$sub.'/insert.'.$sub.'.txt';

(open (IN, "<$in_file")) || die "cannot open the in file\n";
(open (OU, ">$ou_file")) || die "cannot open the ou file\n";
(open (INS, "<$ins_file")) || die "cannot open the insert file\n";

my ($sign, $line, @data, $ins, $mad, %count, %support, @read1, @read2, %identity, %scores, %ins);
my ($mei, $line1, $line2, $c1, $c2, $c3, $c4, $d1, $d2, $d3, $d4, $over, $dist1, $dist2, $dist3, $dist4, $prob1, $prob2, $i, $j, $pA1, $pA2, $dist5, $up1, $up2);
my ($lib1, $lib2, @temp);

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

my $count =1;
for $mei(keys %support)
   {
    print "$count\n" if (!($count %1000));
    $read1_file = $path.'/'.$sub.'/temp/'.$sub.'.'.$strand.'.read1.mei.fa';
    $read2_file = $path.'/'.$sub.'/temp/'.$sub.'.'.$strand.'.read2.mei.fa';
    (open (READ1, ">$read1_file")) || die "cannot open the READ1 file\n";
    (open (READ2, ">$read2_file")) || die "cannot open the READ2 file\n";
    for ($i=0; $i<=$#{$support{$mei}}; $i++)
       {
        @read1 = split("\t", $support{$mei}[$i]);
        print READ1 ">$read1[3]\n$read1[12]\n";
        print READ2 ">$read1[3]\n$read1[12]\n";
       }
    close (READ1);
    close (READ2);
    system("env TMPDIR=$path exonerate --model affine:local --ryo \"INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %s\\n\" $read1_file $read2_file | grep INFO >> $align_file");
    $count++;
   }
print "realignment is finished\n";

(open (ALIGN, "<$align_file")) || die "cannot open the align file\n";
while ($line = <ALIGN>)
   {
    next if  ($line =~ /Command/);
    chomp($line);
    @data  = split(" ", $line);
    if (exists($scores{$data[1]}{$data[5]}) && ($scores{$data[1]}{$data[5]} > $data[10]))
       { ### old alignment is better alignment ###
        next;
       }
    elsif (exists($scores{$data[1]}{$data[5]}) && ($scores{$data[1]}{$data[5]} == $data[10]))
       { ### old alignment is of same quality ###
        $identity{$data[1]}{$data[5]} = $data[3] if ($data[3] > $identity{$data[1]}{$data[5]});
       }
    else
       { ### no old alignment, or olf alignment is worse alignment ###
        $identity{$data[1]}{$data[5]} = $data[3];
        $scores{$data[1]}{$data[5]} = $data[10];
       }
   }

print OU "insertion\tstrand\tread1\tPE1\tup1\tread1_RF\tread1_LR\tread1_NB\tread2\tPE2\tup2\tread2_RF\tread2_LR\tread2_NB\tdist1\tdist2\tdist3\tdist4\tid\tc1\tc2\tc3\tc4\td1\td2\td3\td4\tdist5\n";

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
        if ($read1[3] =~ /^(lib\d+)\_/)
           {
            $lib1 = $1;
           }
        for ($j=$i+1; $j<=$#{$support{$mei}}; $j++)
           {
            @read2 = split("\t", $support{$mei}[$j]);
            $prob2 = $read2[4]; 
            $d1 = ($read2[8] < $read2[9]) ? $read2[8] : $read2[9];
            $d2 = ($read2[8] > $read2[9]) ? $read2[8] : $read2[9];
            $d3 = ($read2[10]<$read2[11]) ? $read2[10] : $read2[11];
            $d4 = ($read2[10]>$read2[11]) ? $read2[10] : $read2[11];
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
            $over  = &over_id($c3, $c4, $read1[3], $d3, $d4, $read2[3]);

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
            print OU "$read1[0]\t$read1[1]\t$line1\t$line2\t$dist1\t$dist2\t$dist3\t$dist4\t$over\t$c1\t$c2\t$c3\t$c4\t$d1\t$d2\t$d3\t$d4\t$dist5\n";
          }
      }
   }

sub over_id
  {
   my $x1 = shift;
   my $x2 = shift;
   my $read1 = shift;
   my $x3 = shift;
   my $x4 = shift;
   my $read2 = shift;
   my ($y, $over);
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
       $over = (exists($identity{$read1}{$read2})) ? $identity{$read1}{$read2} : 0;
      }
   return ($over);
  }      

