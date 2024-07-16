#!/usr/bin/env perl

use warnings;
use strict;
### PE reads for LINE1 insertions ###

### ARGV[0]: output folder to all the subjects ###
my $outpath = $ARGV[0];
### ARGV[1]: subject ID ###
my $sub = $ARGV[1];
### ARGV[2]: RetroSom version control ###
### ARGV[3]: 0, -strand; 1, +strand   ###
my $strand = $ARGV[3];
my $retro = 'retro_v'.$ARGV[2].'_'.$strand;
### ARGV[4]: hg38 or hg19 ###
my $hg = $ARGV[4];
### ARGV[5]: masterpath ###
my $masterpath = ($ARGV[5] =~ /\w/) ? $ARGV[5] : '~/masterpath';

### look for PE supporting reads with split-read anchors ###
my $site_file  = $outpath.'/'.$sub.'/'.$retro.'/LINE/'.$sub.'.LINE.novel.sites';
my $fasta_temp = $outpath.'/'.$sub.'/'.$retro.'/LINE/temp.fa';
my $exo_temp   = $outpath.'/'.$sub.'/'.$retro.'/LINE/temp.exo.txt';
my $seg_temp   = $outpath.'/'.$sub.'/'.$retro.'/LINE/temp.seg.txt';
my $sanch_file = $outpath.'/'.$sub.'/'.$retro.'/LINE/split-reads-anchor.txt';
my $ref_file   = $masterpath.'/refTE/sequence/L1HS.fa';

my $LEN = 20;
(open (IN, "<$site_file")) || die "cannot open the in file\n";
(open (FASTA, ">$fasta_temp")) || die "cannot open the fasta file\n";
(open (OUT, ">$sanch_file")) || die "cannot open the out file\n";

print OUT "read\tdirect\tinsert\tseg\talign\n";
my ($line, @data);
my ($split, $insert, $r1, $direction, $read, $seq);
my (%seg, %exo, $align);

while ($line = <IN>)
   {
    chomp($line);
    @data = split("\t", $line);
    next if ( ($data[0] =~ /^GL/) || ($data[0] =~ /\_/) || ($data[0] =~ /\-/) || ($data[20] eq 'NA') );
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";

    $split = 0;
    if (($data[16] =~ /^(\d+)S\d+M$/) && ($1 > $LEN))
       {
        #print "$data[20]\n";
        $insert = $1;
        $split  = 1;
        $direction = ($data[19] & 0x10) ? 1 : 0;
        $read = $r1.':'.$data[4].'###'.$direction.'###'.$insert;
        $seq = substr($data[20], 0, $insert);
        print FASTA ">$read\n$seq\n";
       }
    elsif (($data[16] =~ /^\d+M(\d+)S$/) && ($1 > $LEN))
       {
        $insert = $1;
        $split  = 1;
        $direction = ($data[19] & 0x10) ? 0 : 1;
        $read = $r1.':'.$data[4].'###'.$direction.'###'.$insert;
        $seq = substr($data[20], length($data[20])-$insert, $insert);
        print FASTA ">$read\n$seq\n";
       }
    }
    
close(FASTA);
system("$masterpath/SEG/seg $fasta_temp 350 -h | awk '{if (/>/) print \$1\"\t\"\$2}' |sed s/complexity=// > $seg_temp");
system("env TMPDIR=$outpath exonerate --model affine:local --bestn 1 --ryo \"INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %qs\\n\" $fasta_temp $ref_file > $exo_temp");

(open (SEG, "<$seg_temp")) || die "cannot open the seg file\n";
(open (EXO, "<$exo_temp")) || die "cannot open the exo file\n";

while ($line = <SEG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $data[0] =~ /^\>(.+)\(/;
    $seg{$1} = $data[1];
   }

while ($line = <EXO>)
   {
    chomp($line);
    if ($line =~ /^INFO/)
       {
        @data = split(" ", $line);
        $exo{$data[1]} = $data[3];
       }
    }

for $read(keys %seg)
   {
    $align = (exists($exo{$read})) ? $exo{$read} : 0;
    @data = split('###', $read);
    print OUT "$data[0]\t$data[1]\t$data[2]\t$seg{$read}\t$align\n";
   }

### build temporal files for checking pAT and blastx ###
my $pe_file    = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.filter.discover';
my $temp1_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.pAT.input';
my $temp2_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.blastx.input';
my $cleavage_file = ($strand) ? $masterpath.'/refTE/position/TTAAAA.'.$hg.'.sorted.pA.bed' : $masterpath.'/refTE/position/TTAAAA.'.$hg.'.sorted.pT.bed';
my $pAT_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.pAT.output';
my $blastx_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.blastx.output';
my $db_file = $masterpath.'/refTE/blastx/L1HS.db';

(open (PE, "<$pe_file")) || die "cannot open the PE discover file\n";
(open (TEMP1, ">$temp1_file")) || die "cannot open the temp1 file\n";
(open (TEMP2, ">$temp2_file")) || die "cannot open the temp2 file\n";

my ($pos1, $pos2);
while ($line = <PE>)
   {
    chomp($line);
    #print "read $count\n";
    @data = split("\t", $line);
    next if ( (($data[0] !~ /^[0-9XY]/) && ($data[0] !~ /^chr[0-9XY]/)) || ($data[0] =~ /\_/) || ($data[0] =~ /\-/) || ($data[3] !~ /L1/) );
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $data[4] = $r1.':'.$data[4];
    print TEMP1 "$data[0]\t$data[1]\t$data[2]\t$data[4]\n";

    $pos1 = ($data[10] < $data[11]) ? $data[10] : $data[11];
    $pos2 = ($data[10] > $data[11]) ? $data[10] : $data[11];
    if ( (($pos1 >= 908) && ($pos2 <=1921)) || (($pos1 >= 1988) && ($pos2 <=5812)) )
       {
        print TEMP2 ">$data[4]\n$data[12]\n";
       }
   }
close (TEMP1);
close (TEMP2);
#system ("closestBed -d -t first -a $temp1_file -b $cleavage_file > $pAT_file");
system ("sort -T $outpath -k1,1 -k2,3n $temp1_file | closestBed -d -t first -a stdin -b $cleavage_file > $pAT_file");
system ("blastx -query $temp2_file -db $db_file -outfmt 6 -out $blastx_file");

### L1 PE reads data matricies ###
my $seg_file  = $outpath.'/'.$sub.'/'.$retro.'/seg.pe.scores';
my $ou_file   = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.pe.LINE.matrix';
my $call_file = $outpath.'/'.$sub.'/'.$retro.'/LINE/'.$sub.'.LINE.novel.calls';
my $aca_file  = $outpath.'/'.$sub.'/fix_ACA/'.$sub.'.discover';
my $acg_file  = $outpath.'/'.$sub.'/fix_ACG/'.$sub.'.discover';
my $aca_all   = $outpath.'/'.$sub.'/fix_ACA/'.$sub.'.L1HS.discover';
my $tag_file  = $outpath.'/'.$sub.'/fix_TAG/'.$sub.'.discover';
my $tag_all   = $outpath.'/'.$sub.'/fix_TAG/'.$sub.'.L1HS.discover';
my $dep1_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.bed';
my $dep2_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.depth';

my $gc5389_file = $outpath.'/'.$sub.'/fix_5389/'.$sub.'.discover';
my $g5533_file  = $outpath.'/'.$sub.'/fix_5533G/'.$sub.'.discover';
my $c5533_file  = $outpath.'/'.$sub.'/fix_5533C/'.$sub.'.discover';
my $at5710_file = $outpath.'/'.$sub.'/fix_5710/'.$sub.'.discover';

my $gc5389_all = $outpath.'/'.$sub.'/fix_5389/'.$sub.'.L1HS.discover';
my $g5533_all = $outpath.'/'.$sub.'/fix_5533G/'.$sub.'.L1HS.discover';
my $at5710_all = $outpath.'/'.$sub.'/fix_5710/'.$sub.'.L1HS.discover';
my %size = (
   'AluYa5' => 282,
   'AluYk13' => 282,
   'AluYb8' => 289,
   'AluYb9' => 289,
   'AluYc1' => 282,
   'AluYa5a2' => 282,
   'L1HS' => 6064,
   'HERVK' => 7536,
   'HERVH' => 7713,
   'SVA_E' => 1382,
   'SVA_F' => 1375,
   'L1' => 6064,
   'Alu' => 289,
   'SVA' => 1382,
    );
my $CUT3 = 1500; ### longest fragment in sequencing ###
my $CUT4 = 200;  ### size cutoff for too short L1 insertions ###

(open (PE  , "<$pe_file  ")) || die "cannot open the pe file\n";
(open (SEG , "<$seg_file ")) || die "cannot open the seg file\n";
(open (CALL, "<$call_file")) || die "cannot open the call file\n";
(open (OU  , ">$ou_file  ")) || die "cannot open the ou file\n";
(open (SITE, "<$site_file")) || die "cannot open the site file\n";
(open (ACA , "<$aca_file ")) || die "cannot open the ACA file\n";
(open (ACG , "<$acg_file ")) || die "cannot open the ACG file\n";
(open (TAG , "<$tag_file ")) || die "cannot open the ACG file\n";
(open (GC5389, "<$gc5389_file")) || die "cannot open the gc5389 file\n";
(open (G5533, "<$g5533_file")) || die "cannot open the g5533 file\n";
(open (C5533, "<$c5533_file")) || die "cannot open the c5533 file\n";
(open (AT5710, "<$at5710_file")) || die "cannot open the at5710 file\n";
(open (ALL1, "<$aca_all  ")) || die "cannot open the all1 file\n";
(open (ALL2, "<$tag_all  ")) || die "cannot open the all2 file\n";
(open (ALL3, "<$gc5389_all ")) || die "cannot open the all3 file\n";
(open (ALL4, "<$g5533_all ")) || die "cannot open the all4 file\n";
(open (ALL5, "<$at5710_all ")) || die "cannot open the all5 file\n";
(open (DEP1, "<$dep1_file")) || die "cannot open the depth file\n";
(open (DEP2, "<$dep2_file")) || die "cannot open the depth file\n";
(open (SPLIT, "<$sanch_file")) || die "cannot open the split file\n";
(open (BLAST, "<$blastx_file")) || die "cannot open the blast file\n";
(open (PAT, "<$pAT_file")) || die "cannot open the blast file\n";

my (@reads, %aca, %acg, %tag, %uniq, %score);
my ($i, $j, $k, @temp, @temp1, @temp2, @temp3, @site, %call);
%seg=();
my ($call, $sr, $pe, $ref, $anchor_alignment, $anchor_length, $anchor_mm, $pos, $neg, $comp, $refpos);
my ($acag, $tag, $a_size, $a_map, $map, $count);
my ($depth, $len, $size, %cord, %depth, %all1, %all2, %all3, %all4, %all5);
my ($short, $cigar, $end3, $end5, $c1, $c2, $c3, $min, $max);
my ($med, $mad, @dep_all, $med_x, $mad_x, @Xdep_all, $med_y, $mad_y, @Ydep_all);
my (%gc5389, %g5533, %c5533, %at5710);
my ($up, $gc5389, $g5533, $c5533, $at5710);
my (%split, $anchor_split, $anchor_direct, $anchor_insert, $anchor_seg, $anchor_map, $anchor_XS, $anchor_AS, $anchor_mapQ);
my (%orf, $ORF, $score, %pat);

$i = 0;
while ($line = <SITE>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $site[++$i] = $r1.':'.$data[4];
   }
print "1: read all support PE reads\n";

while ($line = <SEG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $data[0] =~ /^\>(.+)\(/;
    @temp = split('###', $1);
    if ($temp[1] =~ /L1/)
      {
       $temp[1] = 'L1';
      }
    elsif ($temp[1] =~ /Alu/)
      {
       $temp[1] = 'Alu';
      }
    elsif ($temp[1] =~ /SVA/)
      {
       $temp[1] = 'SVA';
      }
    $line = join ('###', @temp);
    $seg{$line} = $data[1];
   }
print "2: read all seg scores\n";

while ($line = <CALL>)
  {
   chomp($line);
   @data = split("\t", $line);
   if ($data[7] =~ /pe/)
      {
       @temp1 = split(';', $data[7]);
       @temp2 = split(',', $temp1[0]);
       for ($i=1; $i<=$#temp2; $i++)
           {
            $call{$temp2[$i]} = $data[4]."\t".$data[5]."\t".$data[6];
           }
      }
  }
$count = scalar keys %call;
print "3: read all calls $count\n";

while ($line = <ACA>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $aca{$r1.':'.$data[4]} = 1;
   }

while ($line = <ACG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $acg{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL1>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all1{$r1.':'.$data[4]} = 1;
   }

while ($line = <TAG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $tag{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL2>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all2{$r1.':'.$data[4]} = 1;
   }

while ($line = <GC5389>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $gc5389{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL3>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all3{$r1.':'.$data[4]} = 1;
   }

while ($line = <G5533>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $g5533{$r1.':'.$data[4]} = 1;
   }

while ($line = <C5533>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $c5533{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL4>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all4{$r1.':'.$data[4]} = 1;
   }

while ($line = <AT5710>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $at5710{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL5>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all5{$r1.':'.$data[4]} = 1;
   }
print "6: read ACA ACG TAG calls\n";

while ($line = <DEP1>)
   {
    chomp($line);
    next if ($line !~ /PE/);
    @data = split("\t", $line);
    @temp1 = split('__', $data[3]);
    $cord{$temp1[1]} = $data[0]."\t".$data[1]."\t".$data[2];
   }

$i = $j = $k = 0;
while ($line = <DEP2>)
   {
    chomp($line);
    @data = split("\t", $line);
    if ($data[0] =~ /X/)
       {
        $depth{$data[0]."\t".$data[1]."\t".$data[2]} = $data[3] / ($data[2]-$data[1]);
        $Xdep_all[$i++] = $data[3] / ($data[2]-$data[1]);
       }
    elsif ($data[0] =~ /Y/)
       {
        $depth{$data[0]."\t".$data[1]."\t".$data[2]} = $data[3] / ($data[2]-$data[1]);
        $Ydep_all[$j++] = $data[3] / ($data[2]-$data[1]);
       }
    else
       {
        $depth{$data[0]."\t".$data[1]."\t".$data[2]} = $data[3] / ($data[2]-$data[1]);
        $dep_all[$k++] = $data[3] / ($data[2]-$data[1]);
       }
   }
$med_x = &median(@Xdep_all);
$mad_x = &mad(@Xdep_all);
$med_y = &median(@Ydep_all);
$mad_y = &mad(@Ydep_all);
$med = &median(@dep_all);
$mad = &mad(@dep_all);
print "7: read depth file median=$med MAD=$mad \n";
print "7: read depth file Xmedian=$med_x XMAD=$mad_x \n";
print "7: read depth file Ymedian=$med_y YMAD=$mad_y \n";

<SPLIT>;
while ($line = <SPLIT>)
   {
    chomp($line);
    @data = split("\t", $line);
    $split{$data[0]} = $data[1]."\t".$data[2]."\t".$data[3]."\t".$data[4];
   }

while ($line = <PAT>)
   {
    chomp($line);
    @data = split("\t", $line);
    $pat{$data[3]} = $data[7]."\t".$data[8];
   }

while ($line = <BLAST>)
   {
    chomp($line);
    @data = split("\t", $line);
    if (exists($orf{$data[0]}))
       {
        $orf{$data[0]} = $data[11] if ($data[11] > $orf{$data[0]});
       }
    else
       {
        $orf{$data[0]} = $data[11];
       }
   }
print "8: read the splitanchor, pAT and blastx file\n";

$count =1;
print OU "chr\tcord1\tcord2\tread\trefpos1\trefpos2\tpos\tneg\tseg\tsupport\tsr\tpe\tmap\tanchor_align\tanchor_len\tanchor_mm\tanchor_split\tanchor_direct\tanchor_insert\tanchor_seg\tanchor_map\tanchor_XS\tanchor_AS\tanchor_mapQ\tref\tACAG\tTAG\tdepth\tmap_len\tshort\tend3\tend5\tc1\tc2\tc3\tdirection\tcenter\tGC5389\tG5533\tC5533\tAT5710\tupstream\tgap\tpAT\tdist\tORF\tXscore\trefpos\n";
while ($line = <PE>)
   {
    chomp($line);
    #print "read $count\n";
    $count ++;
    @data = split("\t", $line);
    next if ( (($data[0] !~ /^[0-9XY]/) && ($data[0] !~ /^chr[0-9XY]/)) || ($data[0] =~ /\_/) || ($data[0] =~ /\-/) );
    $read = $data[4].'###'.$data[3].'###'.$data[7].'###'.$data[8].'###'.$data[9];
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $data[4] = $r1.':'.$data[4];

    next if ($data[3] !~ /L1/);
    if (exists($call{$data[4]}))
      {### number of support ###
       $ref = 0;
       ($call, $sr, $pe) = split("\t", $call{$data[4]});
      }
    else
      { ### reference insertions reads ###
       $call = $sr = $pe = 0;
       $ref = 1;
      }

    ### short and 3end ###
    $c1 = ($data[5] eq '+') ? 1 : -1;
    $c2 = ($data[6] eq '+') ? 1 : -1;
    $c3 = ($data[15] eq '+') ? 1 : -1;
    $direction = -1 * $c1 * $c2 * $c3;
    $up = $c1 * $direction;
    $short = $end3 = 0;
    $refpos= abs(($data[10] + $data[11])/2 - 3000);
    $min = ($data[11] < $data[10]) ? $data[11] : $data[10];
    $max = ($data[11] > $data[10]) ? $data[11] : $data[10];
    $anchor_XS = $data[22];
    $anchor_AS = $data[23];
    $anchor_mapQ = $data[24];

    if ($data[5] eq '+')
       { ### upstream read ###
        if ($direction > 0)
           { ### positive strand 5->3 ###
            if (($data[3] !~ /Alu/) && ($min > ($size{$data[3]}-$CUT4) ))
               {
                $short = 1;
               }
           }
        else
           { ### negative strand 3->5 ###
            if ($max  < ($size{$data[3]} - $CUT3))
               {
                $end3 = 1;
               }
           }
       }
    else
       { ### downstream read ###
        if ($direction > 0)
           { ### positive strand 5->3 ###
            if ($max  < ($size{$data[3]} - $CUT3))
               {
                $end3 = 1;
               }
           }
        else
           { ### negative strand 3->5 ###
            if (($data[3] !~ /Alu/) && ($min > ($size{$data[3]}-$CUT4) ))
               {
                $short = 1;
               }
           }
       }

    ### 5' transduction if the 5' is truncated ###
    $end5 = 0;
    $cigar = $data[13];
    if ( $data[10] > $data[11] )
        { ### L1 5' is at the end of read ###
         if ( ($data[11] > 10) && ( (length($data[12]) - $data[9]) > 10 ))
            { ### L1 is missing 5'end and there is extra sequence on 5'end ###
             if ($cigar =~ /(\d+)S$/)
                {
                 $end5 = 1 if ( abs((length($data[12]) - $data[9]) - $1) > 4 );
                }
             else
                {
                 $end5 = 1;

                }
            }
        }
     elsif ($data[10] < $data[11])
        {
         if ( ($data[10] > 10) && ($data[8] > 10) )
            { ### L1 is missing 5'end and there is extra sequence on 5 end ###
             if ($cigar =~ /^(\d+)S/)
                {
                 $end5 = 1 if ( abs($data[8] - $1) > 4 );
                }
             else
                {
                 $end5 = 1;

                }
            }
        }

    $pos = 1;
    $neg = 1;
    if (exists($all1{$data[4]}))
       {
        $acag = (exists($aca{$data[4]})|| exists($acg{$data[4]})) ? 1 : 0;
       }
    else 
       {
        $acag = 'NA';
       }
    if (exists($all2{$data[4]}))
       {
        $tag = exists($tag{$data[4]}) ? 1 : 0;
       }
    else
       {
        $tag = 'NA';
       }
    if (exists($all3{$data[4]}))
       {
        $gc5389 = exists($gc5389{$data[4]}) ? 1 : 0;
       }
    else
       {
        $gc5389 = 'NA';
       }
    if (exists($all4{$data[4]}))
       {
        $g5533 = exists($g5533{$data[4]}) ? 1 : 0;
       }
    else
       {
        $g5533 = 'NA';
       }
    if (exists($all4{$data[4]}))
       {
        $c5533 = exists($c5533{$data[4]}) ? 1 : 0;
       }
    else
       {
        $c5533 = 'NA';
       }
    if (exists($all5{$data[4]}))
       {
        $at5710 = exists($at5710{$data[4]}) ? 1 : 0;
       }
    else
       {
        $at5710 = 'NA';
       }
    if ($data[0] =~ /X/)
       {
        $depth = ( $depth{$cord{$data[4]}} - $med_x ) / $mad_x;
       }
    elsif ($data[0] =~ /Y/)
       {
        $depth = ( $depth{$cord{$data[4]}} - $med_y ) / $mad_y;
       }
    else
       {
        $depth = ( $depth{$cord{$data[4]}} - $med ) / $mad;
       }

    #print "$data[4]\t$depth\n" if (!$depth);

    $comp = $seg{$read}; ### seq complexity ###
    $map  = $data[7];
    $len  = abs($data[9] - $data[8]);  ### frag length mappable to TE ###
    $size = length($data[12]);      ### sequence length ###
    $len /= $size;
    $anchor_alignment = &check_XA($data[17]); ### ignoring alt contigs in XA tags ###
    
    $a_size = $a_map = 0;
    while ($data[16] =~ /(\d+)/g)
       {
        $a_size += $1; # anchor size, all parts #
       }
    while ($data[16] =~ /(\d+)M/g)
       {
        $a_map += $1;  # anchor, all mappable size #
       }
    $anchor_length = $a_map / $a_size; ### mapped portion on anchor reads ###
    $anchor_mm = $data[18]/$a_size;    ### number of anchor mismatches ###

    if (exists($split{$data[4]}))
       {
        $anchor_split = 1;
        ($anchor_direct, $anchor_insert, $anchor_seg, $anchor_map) = split("\t", $split{$data[4]});
        $anchor_insert /= length($data[12]);
       }
    else
       {
        $anchor_split = 0;
        $anchor_direct = $anchor_insert = $anchor_seg = $anchor_map = 'NA';
       }

    if ( exists($orf{$data[4]})) 
       {
        $ORF = 1;
        $score = $orf{$data[4]};
       }
    else
       {
        $ORF = $score = 0;
       }
    $pos1 = ( $data[10] + $data[11] ) / 2;
    
    if (!exists($score{$data[4]}))
       {
        $uniq{$data[4]} = "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[10]\t$data[11]\t$pos\t$neg\t$comp\t$call\t$sr\t$pe\t$map\t$anchor_alignment\t$anchor_length\t$anchor_mm\t$anchor_split\t$anchor_direct\t$anchor_insert\t$anchor_seg\t$anchor_map\t$anchor_XS\t$anchor_AS\t$anchor_mapQ\t$ref\t$acag\t$tag\t$depth\t$len\t$short\t$end3\t$end5\t$c1\t$c2\t$c3\t$direction\t$refpos\t$gc5389\t$g5533\t$c5533\t$at5710\t$up\t$data[25]\t$pat{$data[4]}\t$ORF\t$score\t$pos1";
        $score{$data[4]} = $map;
       }
    elsif ($map > $score{$data[4]})
       {
        $uniq{$data[4]} = "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[10]\t$data[11]\t$pos\t$neg\t$comp\t$call\t$sr\t$pe\t$map\t$anchor_alignment\t$anchor_length\t$anchor_mm\t$anchor_split\t$anchor_direct\t$anchor_insert\t$anchor_seg\t$anchor_map\t$anchor_XS\t$anchor_AS\t$anchor_mapQ\t$ref\t$acag\t$tag\t$depth\t$len\t$short\t$end3\t$end5\t$c1\t$c2\t$c3\t$direction\t$refpos\t$gc5389\t$g5533\t$c5533\t$at5710\t$up\t$data[25]\t$pat{$data[4]}\t$ORF\t$score\t$pos1";
        $score{$data[4]} = $map;
       }
   }

@reads = ();
$i = 0;
for $read(keys %uniq)
   {
    @data = split("\t", $uniq{$read});
    $reads[$i][0] = $data[0];
    $reads[$i][1] = $data[1];
    $reads[$i][2] = $data[2];
    $reads[$i][3] = $uniq{$read};
    $i++;
   }

@reads = sort{ $a->[0] cmp $b ->[0] || $a->[1] <=> $b ->[1] || $a->[2] <=> $b ->[2]} @reads;

for ($i=0; $i<=$#reads; $i++)
   {
    print OU "$reads[$i][3]\n";
   }
close (OU);

sub median
   {
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
       {
        return $vals[int($len/2)];
       }
    else #even
       {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
       }
   }

sub mad
   {
    my @vals = @_;
    my $med = &median(@vals);
    my ($i, $mad, @diff);
    my $constant = 1.4826;
    for ($i=0; $i<=$#vals; $i++)
       {
        $diff[$i] = abs($vals[$i] - $med);
       }
    $mad = &median(@diff);
    return($mad*$constant);
   }

sub check_XA
   { ### no alternative alignment, return 1 ###
    my $tag = shift;
    return (1) if ($tag =~ /NULL/);  ### good unique alignment for anchors ###
    $tag =~ /XA:Z:(.+)$/;
    my $str = $1;
    my ($xa, @temp);
  
    @temp = split(/\;/, $str);
    foreach $xa(@temp)
        {
         if ($xa && ($xa !~ /alt/))
            {
             return (0);
            }
        }
    return(1);
   }
