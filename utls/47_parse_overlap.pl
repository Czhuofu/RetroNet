#!/usr/bin/env perl

use warnings;
use strict;
use Math::CDF;

my $retro = $ARGV[0];
my $TE    = $ARGV[1];
my $path  = $ARGV[2];
my $sub   = $ARGV[3].'_Brain';
my $id_CUT= 0.95;
my $masterpath= '/home/xwzhu/levinson/RetroSom';
my $nov_file  = $path.'/'.$sub.'/'.$retro.'/'.$TE.'/'.$sub.'.'.$TE.'.noref.calls';
my $over_file = $path.'/'.$sub.'/'.$retro.'/'.$TE.'/'.$sub.'.'.$TE.'.over.calls';
my $out_file  = $path.'/'.$sub.'/'.$retro.'/'.$TE.'/'.$sub.'.'.$TE.'.overlap.txt';
my $pe_file   = $path.'/'.$sub.'/'.$retro.'/'.$TE.'/'.$sub.'.'.$TE.'.novel.sites';
my $sr_file   = $path.'/'.$sub.'/'.$retro.'/'.$sub.'.'.$TE.'.sr.tabe.discover';
my $r1_file   = $masterpath.'/temp/'.$sub.'.retro_v'.$ARGV[0].'.read1.fa';
my $r2_file   = $masterpath.'/temp/'.$sub.'.retro_v'.$ARGV[0].'.read2.fa';
my $align_file= $masterpath.'/temp/'.$sub.'.retro_v'.$ARGV[0].'.align.temp';
my ($line, @data, %all, %calls, %count, %PE_all, %sr_uniq);
my (%novel, @sr, @temp, @nums, @dat2, @read1, @read2, %iden);
my ($start, $size, $left1, $left2, $left3, $right1, $right2, $right3, $seq, $seq1, $seq2, $prob, %prob);
my ($i, $j, $site, $GAP, $MAP, $MM, $identity, $substr1, $substr2, $sr);
my $over_count  = 0;
(open (OUT , ">$out_file")) || die "cannot open the out file $out_file\n";
(open (OVER, ">$over_file"))|| die "cannot open the over file\n";
(open (PE  , "<$pe_file ")) || die "cannot open the PE file\n";
(open (SR  , "<$sr_file ")) || die "cannot open the SR file\n";
(open (NOV , "<$nov_file")) || die "cannot open the novel file\n";
my %tagTE = (
   'ALU' => 'Alu',
   'LINE' => 'L1HS',
  );
my $tag = $tagTE{$TE};
$i = 0;
while ($line = <PE>)
   {
    chomp($line);
    $i ++;
    @data = split("\t", $line);
    $PE_all{$i} = $data[8]."\t".$data[9]."\t".$data[10]."\t".$data[11]."\t".$data[12]."\t".$data[7]; 
   }

while ($line = <SR>)
   {
    chomp($line);
    @data = split("\t", $line);
    $sr_uniq{$data[4]} = $data[11]."\t".$data[12]."\t".$data[13]."\t".$data[14]."\t".$data[17]."\t".$data[8] if ($data[3] =~ /$tag/);
   }

while ($line = <NOV>)
   {
    chomp($line);
    @data = split("\t", $line);
    $site = $data[0]."\t".$data[1]."\t".$data[2]."\t".$data[4]."\t".$data[5]."\t".$data[6];
    # next if ( ($data[4] < 2) || ($data[4] > 3) );
    $calls{$site} = $line;
    $count{$site} = 0;
    if ( ($data[7] =~ /pe\,/) && ($data[7] =~ /sr\,/) )
       {
        @temp = split(';', $data[7]);
        for ($j=0; $j<=$#temp; $j++)
           {
            if ($temp[$j] =~ /pe/)
              {
               @nums = split(',', $temp[$j]);
               for ($i=4; $i<=$#nums; $i++)
                  {
                   if ($nums[$i])
                     {
                      @dat2 = split('_', $nums[$i]);
                      if (exists($PE_all{$dat2[0]}))
                         {
                          $all{$site}[$count{$site}][0] = $PE_all{$dat2[0]}; 
                          $all{$site}[$count{$site}][1] = $dat2[0];
                          $count{$site} ++;
                         }
                     }
                  }
              }
            elsif ($temp[$j] =~ /sr/)
              {
               @nums = split(',', $temp[$j]);
               @sr   = split(/\|/, $nums[2]);
               foreach $sr(@sr)
                  {
                   if (exists($sr_uniq{$sr}))
                      {
                       $all{$site}[$count{$site}][0] = $sr_uniq{$sr};
                       $all{$site}[$count{$site}][1] = $sr;
                       $count{$site} ++;
                      }
                  }
              }
           }
       }
    elsif ($data[7] =~ /pe\,/)
       {
        @nums = split(',', $data[7]);
        for ($i=4; $i<=$#nums; $i++)
           {
            if ($nums[$i] && ($nums[$i] =~ /\_/))
              {
               @dat2 = split('_', $nums[$i]);
               if (exists($PE_all{$dat2[0]}))
                  {
                   $all{$site}[$count{$site}][0] = $PE_all{$dat2[0]};  
                   $all{$site}[$count{$site}][1] = $dat2[0]; 
                   $count{$site} ++;
                  }
              }
           }
       }
    else
       {
        @sr = split(/[\,\|\;]/, $data[7]);
        foreach $sr(@sr)
          {
           next if (($sr eq 'sr') || ($sr =~ /^\d+$/));
           if (exists($sr_uniq{$sr}))
              {
               $all{$site}[$count{$site}][0] = $sr_uniq{$sr};
               $all{$site}[$count{$site}][1] = $sr;
               $count{$site} ++;
              }
          }
       }
   }

print OUT "chr\tcord1\tcord2\tsupport\tSR\tPE\tread1\tread2\tsize\tmap\tgap\tmismatch\tidentity\tmap1\tmap2\tpvalue\n";

my $stop=0;
for $site(keys %calls)
   {
    $stop ++;
    @data = split("\t", $site);
    if ( ($data[3] != 2) && ($data[3] != 3) )
       {
        print OVER "$calls{$site}\n";
        next;
       }
    for ($i=0; $i<($count{$site} - 1); $i++)
       {
        @read1 = split("\t", $all{$site}[$i][0]);
        $seq   = substr($read1[4], $read1[0], $read1[1]-$read1[0]);
        $seq1  = ($read1[2] < $read1[3]) ? $seq : &rc($seq);
        $left1 = ($read1[2] < $read1[3]) ? $read1[2] : $read1[3];
        $right1= ($read1[2] < $read1[3]) ? $read1[3] : $read1[2];
        for ($j=$i+1; $j<$count{$site}; $j++)
           {
            @read2 = split("\t", $all{$site}[$j][0]);
            $seq   = substr($read2[4], $read2[0], $read2[1]-$read2[0]);
            $seq2  = ($read2[2] < $read2[3]) ? $seq : &rc($seq);
            $left2 = ($read2[2] < $read2[3]) ? $read2[2] : $read2[3];
            $right2= ($read2[2] < $read2[3]) ? $read2[3] : $read2[2];
            $left3 = ($left1  < $left2 ) ? $left2  : $left1;
            $right3= ($right1 < $right2) ? $right1 : $right2;
            if (($right3 - $left3) > 20)
               {
                (open (R1, ">$r1_file")) || die "cannot open read1 file\n";
                (open (R2, ">$r2_file")) || die "cannot open read2 file\n";
                $start = $left3 - $left1;
                $size  = (($right3 - $left3) > (length($seq1) - $start)) ? (length($seq1) - $start) : ($right3 - $left3);
                $substr1 = substr($seq1, $start, $size);
                $start = $left3 - $left2;
                $size  = (($right3 - $left3) > (length($seq2) - $start)) ? (length($seq2) - $start) : ($right3 - $left3);
                $substr2 = substr($seq2, $start, $size);
                print R1 ">read1\n$substr1\n";
                print R2 ">read2\n$substr2\n";
                system("$masterpath/utls/18_align_overlap.sh $r1_file $r2_file $align_file");
                (open (ALIGN, "<$align_file")) || die "cannot open the align file\n";
                $identity = -1;
                $GAP = $MAP = $MM = 0;
                while ($line = <ALIGN>)
                   {
                    if ($line =~ /vulgar/)
                       {
                        while ($line =~ /G (\d+) (\d+)/g)
                           {
                            $GAP += $1;
                            $GAP += $2;
                           }
                        while ($line =~ /M (\d+) (\d+)/g)
                           {
                            $MAP += $1;
                           }
                        chomp($line = <ALIGN>);
                        @data = split(" ", $line);
                        $MM = $MAP * (1 - $data[3] / 100);
                        $MM = int($MM + 0.5);
                        $identity = 1 - ($MM + $GAP) / ($GAP + $MAP);
                        $prob = 1 - &Math::CDF::pbinom($MM+$GAP,$GAP+$MAP,0.01);
                        $identity = substr($identity, 0, 6);
                        last;
                       }
                   }
                $prob{$site} = $prob;
                $iden{$site} = $identity;
                $over_count ++; 
                print OUT "$site\t$all{$site}[$i][1]\t$all{$site}[$j][1]\t$size\t$MAP\t$GAP\t$MM\t$identity\t$read1[5]\t$read2[5]\t$prob\n";
               }
           }
       }    
   }

for $site(keys %calls)
   {
    @data = split("\t", $site);
    if ( (($data[3] == 2) || ($data[3] == 3)) && (!exists($prob{$site})) )
       {
        print OVER "$calls{$site}\n";
       }
    #elsif ( exists($prob{site}) && ($prob{$site} > (0.05/$over_count)) && ($iden{$site} > $id_CUT) )
    elsif ( exists($prob{$site}) && ($prob{$site} > 0.05) && ($iden{$site} > $id_CUT) )
       {
        print OVER "$calls{$site}\n";
       }
   }

sub rc 
  { ### reverse complement ###
   my $seq = shift;
   my %RCbase = (
      'A' => 'T',
      'T' => 'A',
      'G' => 'C',
      'C' => 'G',
      'a' => 'T',
      't' => 'A',
      'g' => 'C',
      'c' => 'G',
      'n' => 'N',
      'N' => 'N',
     );
   my @bases = split("", $seq);
   my @output;
   my ($seqout, $i);

   for ($i=0; $i<=$#bases; $i++)
     {
      $output[$i] = $RCbase{$bases[length($seq)-$i-1]};
     }
   $seqout=join("", @output);
   return($seqout);
  }
