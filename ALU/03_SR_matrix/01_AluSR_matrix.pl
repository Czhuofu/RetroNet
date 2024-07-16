#!/usr/bin/env perl

use warnings;
use strict; 

### Feature matrix for Alu SR supporting reads ###
### ARGV[0]: output outpath to all the subjects ###
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

### build temporal files for checking pAT ###
my $sr_file    = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.tabe.discover';
my $temp1_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.ALU.pAT.input';
my $cleavage_file = ($strand) ? $masterpath.'/refTE/position/TTAAAA.'.$hg.'.sorted.pA.bed' : $masterpath.'/refTE/position/TTAAAA.'.$hg.'.sorted.pT.bed';
my $pAT_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.ALU.pAT.output';

(open (SR, "<$sr_file")) || die "cannot open the PE discover file\n";
(open (TEMP1, ">$temp1_file")) || die "cannot open the temp1 file\n";
my ($line, @data, $r1, @temp, $read);

while ($line = <SR>)
   {
    chomp($line);
    @data = split("\t", $line);
    next if ( (($data[0] !~ /^[0-9XY]/) && ($data[0] !~ /^chr[0-9XY]/)) || ($data[0] =~ /\_/) || ($data[0] =~ /\-/) || ($data[3] !~ /Alu/) );
    ($r1, @temp) = split(/\:/, $data[4]);
    $data[4] = join(':', @temp);
    $read = $data[4]."\t".$data[3]."\t".$data[6]."\t".$data[11]."\t".$data[12]."\t".$data[13]."\t".$data[14];
    $data[4] = $r1.':'.$data[4];
    print TEMP1 "$data[0]\t$data[1]\t$data[2]\t$data[4]\n";
   }
close (TEMP1);
#system ("closestBed -d -t first -a $temp1_file -b $cleavage_file > $pAT_file");
system ("sort -k1,1 -k2,3n $temp1_file | closestBed -d -t first -a stdin -b $cleavage_file > $pAT_file");

### Feature matrix for Alu SR supporting reads ###
my $seg_file  = $outpath.'/'.$sub.'/'.$retro.'/seg.sr.scores';
my $sr_old    = $outpath.'/'.$sub.'/retro_v'.$ARGV[2].'/'.$sub.'.sr2.discover';
my $ou_file   = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.ALU.matrix';
my $call_file = $outpath.'/'.$sub.'/'.$retro.'/ALU/'.$sub.'.ALU.novel.calls';

my $ggcg_file = $outpath.'/'.$sub.'/fix_GGCG/'.$sub.'.sr.discover';
my $ggcg_all = $outpath.'/'.$sub.'/fix_GGCG/'.$sub.'.sr.Alu.discover';
my $gagc_file = $outpath.'/'.$sub.'/fix_GAGC/'.$sub.'.sr.discover';
my $gagc_all = $outpath.'/'.$sub.'/fix_GAGC/'.$sub.'.sr.Alu.discover';
my $cc_file = $outpath.'/'.$sub.'/fix_CC/'.$sub.'.sr.discover';
my $cc_all = $outpath.'/'.$sub.'/fix_CC/'.$sub.'.sr.Alu.discover';
my $tg_file = $outpath.'/'.$sub.'/fix_TG/'.$sub.'.sr.discover';
my $tg_all = $outpath.'/'.$sub.'/fix_TG/'.$sub.'.sr.Alu.discover';
my $cg_file = $outpath.'/'.$sub.'/fix_CG/'.$sub.'.sr.discover';
my $cg_all = $outpath.'/'.$sub.'/fix_CG/'.$sub.'.sr.Alu.discover';
my $at_file = $outpath.'/'.$sub.'/fix_AT/'.$sub.'.sr.discover';
my $at_all = $outpath.'/'.$sub.'/fix_AT/'.$sub.'.sr.Alu.discover';
my $Yb_file = $outpath.'/'.$sub.'/fix_AluYb/'.$sub.'.sr.discover';
my $Yb_all = $outpath.'/'.$sub.'/fix_AluYb/'.$sub.'.sr.Alu.discover';

my $dep1_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.bed';
my $dep2_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.depth';

my $CUT = 200;
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

(open (SR  , "<$sr_file  ")) || die "cannot open the sr file\n";
(open (SEG , "<$seg_file ")) || die "cannot open the seg file\n";
(open (CALL, "<$call_file")) || die "cannot open the call file\n";
(open (OU  , ">$ou_file  ")) || die "cannot open the ou file\n";

(open (GGCG, "<$ggcg_file")) || die "cannot open the GGCG file\n";
(open (GAGC, "<$gagc_file")) || die "cannot open the GAGC file\n";
(open (CC  , "<$cc_file  ")) || die "cannot open the CC file\n";
(open (TG  , "<$tg_file  ")) || die "cannot open the TG file\n";
(open (CG  , "<$cg_file  ")) || die "cannot open the CG file\n";
(open (AT  , "<$at_file  ")) || die "cannot open the AT file\n";
(open (Yb  , "<$Yb_file  ")) || die "cannot open the AluYb file\n";

(open (ALL1, "<$ggcg_all ")) || die "cannot open the all1 file\n";
(open (ALL2, "<$gagc_all ")) || die "cannot open the all2 file\n";
(open (ALL3, "<$cc_all   ")) || die "cannot open the all3 file\n";
(open (ALL6, "<$at_all   ")) || die "cannot open the all4 file\n";
(open (ALL7, "<$Yb_all   ")) || die "cannot open the all5 file\n";

(open (DEP1, "<$dep1_file")) || die "cannot open the depth file\n";
(open (DEP2, "<$dep2_file")) || die "cannot open the depth file\n";
(open (SOLD, "<$sr_old   ")) || die "cannot open the sr old file\n";
(open (PAT , "<$pAT_file ")) || die "cannot open the sr old file\n";

my ($temp2_file);
my (@reads, %aca, %acg, %tag, %uniq, %score);
my ($i, $j, $k, @temp1, @temp2, @temp3, @site, %seg, %call);
my ($call, $sr, $pe, $ref, $anchor_alignment, $anchor_length, $anchor_mm, $pos, $neg, $comp);
my ($tag, $map, $count);
my ($depth, $len, $size, %cord, %depth);
my ($te, $short, $end3, $end5, $direction, $min, $max);
my ($med, $mad, @dep_all, $med_x, $mad_x, @Xdep_all, $med_y, $mad_y, @Ydep_all);
my ($start, $end, $cigar, $ratio, $refpos, $anchor, $insert);
my ($A_pair, $A_insert, $a_size, $a_map,  $A_mm, $A_MapQ, $A_AS, $A_XS, %old_sr, $old);
my (%pat, %ggcg, %gagc, %cc, %tg, %cg, %at, %Yb);
my ($up, $ggcg, $gagc, $cc, $tg, $cg, $at, $yb);
my (%all1, %all2, %all3, %all4, %all5, %all6, %all7);

while ($line = <SOLD>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $old_sr{$data[3]."\t".$r1.':'.$data[4]} = 1;
   }

while ($line = <SEG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $data[0] =~ /^\>(.+)\>(\w+)\>(\w+)\>(\w+)\>(\w+)\>(\w+)\>(\w+)\(/;
    $tag = $1;
    $i   = $2;
    $cigar = $3;
    $seg{$tag."\t".$i."\t".$cigar."\t".$4."\t".$5."\t".$6."\t".$7} = $data[1];
   }
print "2: read all seg scores\n";

while ($line = <CALL>)
  {
   chomp($line);
   @data = split("\t", $line);
   if ($data[7] =~ /sr/)
      {
       if ($data[7] =~ /pe/)
          {
           @temp1 = split(';', $data[7]);
           @temp2 = split(',', $temp1[1]);
           for ($i=1; $i<=$#temp2; $i++)
              {
               $call{$temp2[$i]} = $data[4]."\t".$data[5]."\t".$data[6];
              }
          }
       else
          {
           @temp2 = split(',', $data[7]);
           for ($i=1; $i<=$#temp2; $i++)
              {
               $call{$temp2[$i]} = $data[4]."\t".$data[5]."\t".$data[6];
              }
          }
      }
   }
$count = scalar keys %call;
print "3: read all calls $count\n";

while ($line = <GGCG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $ggcg{$r1.':'.$data[4]} = 1;
   }

while ($line = <GAGC>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $gagc{$r1.':'.$data[4]} = 1;
   }

while ($line = <CC>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $cc{$r1.':'.$data[4]} = 1;
   }

while ($line = <TG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $tg{$r1.':'.$data[4]} = 1;
   }

while ($line = <CG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $cg{$r1.':'.$data[4]} = 1;
   }

while ($line = <AT>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $at{$r1.':'.$data[4]} = 1;
   }

while ($line = <Yb>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $Yb{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL1>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $all1{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL2>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $all2{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL3>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $all3{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL6>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $all6{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL7>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $all7{$r1.':'.$data[4]} = 1;
   }
print "6: read Allele calls\n";

while ($line = <DEP1>)
   {
    chomp($line);
    next if ($line !~ /SR/);
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
$med = &median(@dep_all);
$mad = &mad(@dep_all);
$med_x = &median(@Xdep_all);
$mad_x = &mad(@Xdep_all);
$med_y = &median(@Ydep_all);
$mad_y = &mad(@Ydep_all);
print "7: read depth file median=$med MAD=$mad medianX=$med_x MADX=$mad_x medianY=$med_y MADY=$mad_y\n";

while ($line = <PAT>)
   {
    chomp($line);
    @data = split("\t", $line);
    $pat{$data[3]} = $data[7]."\t".$data[8];
   }

$count = 0;
print OU "chr\tcord1\tcord2\tread\trefpos1\trefpos2\tpos\tneg\tseg\tsupport\tsr\tpe\tmap\tref\tGGCG\tGAGC\tCC\tTG\tCG\tAT\tYB\tdepth\tmap_len\tmap_size\tmap_ratio\tshort\tend3\tend5\tdirection\trefpos\tanchor\tinsert\tA_pair\tA_insert\tA_mm\tA_MapQ\tA_AS\tA_XS\toldSR\tgap\tpAT\tdist\n";
while ($line = <SR>)
   {
    chomp($line);
    #print "read $count\n";
    $count ++;
    @data = split("\t", $line);
    next if ( (($data[0] !~ /^[0-9XY]/) && ($data[0] !~ /^chr[0-9XY]/)) || ($data[0] =~ /\_/) || ($data[0] =~ /\-/) );
    ($r1, @temp) = split(/\:/, $data[4]);
    $data[4] = join(':', @temp);
    $read = $data[4]."\t".$data[3]."\t".$data[6]."\t".$data[11]."\t".$data[12]."\t".$data[13]."\t".$data[14];
    $data[4] = $r1.':'.$data[4];

    next if ($data[3] !~ /Alu/);
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
    $direction = ($data[13] < $data[14]) ? 1 : -1;
    $short = $end3 = 0;
    if ($data[9] > 0)
      {
       $anchor = 1;
      }
    else
      {
       $anchor = -1;
      }

    ($end3, $end5, $insert) = &parse_end($data[11], $data[12], $data[13], $data[14], $data[17], $data[6]);

    $pos = 1;
    $neg = 1;
    if (exists($all1{$data[4]}))
       {
        $ggcg = exists($ggcg{$data[4]}) ? 1 : 0;
       }
    else
       {
        $ggcg = 'NA';
       }
    if (exists($all2{$data[4]}))
       {
        $gagc = exists($gagc{$data[4]}) ? 1 : 0;
       }
    else
       {
        $gagc = 'NA';
       }
    if (exists($all3{$data[4]}))
       {
        $cc = exists($cc{$data[4]}) ? 1 : 0;
        $tg = exists($tg{$data[4]}) ? 1 : 0;
        $cg = exists($cg{$data[4]}) ? 1 : 0;
       }
    else
       {
        $cc = $tg = $cg = 'NA';
       }
    if (exists($all6{$data[4]}))
       {
        $at = exists($at{$data[4]}) ? 1 : 0;
       }
    else
       {
        $at = 'NA';
       }
    if (exists($all7{$data[4]}))
       {
        $yb = exists($Yb{$data[4]}) ? 1 : 0;
       }
    else
       {
        $yb = 'NA';
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
    $map  = $data[8];
    $len  = abs($data[12] - $data[11]);  ### frag length mappable to TE ###
    $size = split("", $data[17]);        ### sequence length ###
    $ratio= $len / $size;
    #$center = abs(($data[13] + $data[14])/2 - 3000);
    $refpos = ($data[13] + $data[14])/2;

    ###A_pair\tA_insert\tA_mm\tA_MapQ\tA_AS\tA_XS
    $A_pair = ($data[18] & 0x2) ? 1 : 0;
    $old = 1;
    if (exists($old_sr{$data[3]."\t".$data[4]}))
        {
         $A_pair= 1;
         $old   = 0;
        }

    $a_size = $a_map = 0;
    while ($data[19] =~ /(\d+)/g)
       {
        $a_size += $1;
       }
    while ($data[19] =~ /(\d+)M/g)
       {
        $a_map += $1;
       }

    $A_insert = 1 - $a_map/$a_size;
    $A_mm     = $data[20]/$a_size;
    $A_MapQ   = $data[21];
    $A_AS     = $data[22];
    $A_XS     = $data[23]; 
    
    ### unique entries ###
    if (!exists($score{$data[4]}))
       {
        $uniq{$data[4]} = "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[13]\t$data[14]\t$pos\t$neg\t$comp\t$call\t$sr\t$pe\t$map\t$ref\t$ggcg\t$gagc\t$cc\t$tg\t$cg\t$at\t$yb\t$depth\t$len\t$size\t$ratio\t$short\t$end3\t$end5\t$direction\t$refpos\t$anchor\t$insert\t$A_pair\t$A_insert\t$A_mm\t$A_MapQ\t$A_AS\t$A_XS\t$old\t$data[24]\t$pat{$data[4]}";
        $score{$data[4]} = $map * $len;
       }
    elsif (($map * $len) > $score{$data[4]})
       {
        $uniq{$data[4]} = "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[13]\t$data[14]\t$pos\t$neg\t$comp\t$call\t$sr\t$pe\t$map\t$ref\t$ggcg\t$gagc\t$cc\t$tg\t$cg\t$at\t$yb\t$depth\t$len\t$size\t$ratio\t$short\t$end3\t$end5\t$direction\t$refpos\t$anchor\t$insert\t$A_pair\t$A_insert\t$A_mm\t$A_MapQ\t$A_AS\t$A_XS\t$old\t$data[24]\t$pat{$data[4]}";
        $score{$data[4]} = $map * $len;
       }   
   }
print "$count SR reads \n";

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

sub parse_end
   {
    my $p1 = shift;
    my $p2 = shift;
    my $p3 = shift;
    my $p4 = shift;
    my $seq = shift;
    my $cigar = shift;
    my ($end3, $end5, $insert);

    $end3 = $end5 = 0;
    if ($cigar =~ /\d+M\d+S$/)
       {
        $insert = 1;
        if ( ($p3 > $p4) && ($p3 < 275) )
           {
            $end3 = 1;   ### no 3' ###
           }
        elsif (($p4 < 275) && ($p3 < $p4) && ((length($seq)-$p2) > 5) )
           {
            $end3 = 1; ### something inserted before the full length 3'end ###
           }
        }
    elsif ($cigar =~ /^\d+S\d+M/)
        {
         $insert = -1;
         if ( ($p4 < 275) && ($p3 < $p4) )
            {
             $end3 = 1;
            }
         elsif ( ($p3 > $p4) && ($p3 < 275) && ($p1 > 5) )
            {
             $end3 = 1;
            }
        }
    if ( ($p3 < $p4) && ($p3 > 4) && ($p1>4) )
       {
        $end5 = 1;
       }
    elsif ( ($p3 > $p4) && ((length($seq)-$p2)>4) && ($p4 > 4) )
       {
        $end5 = 1;
       }
     return($end3, $end5, $insert);
    }
