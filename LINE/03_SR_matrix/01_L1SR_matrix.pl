#!/usr/bin/env perl

use warnings;
use strict; 
### SR reads for LINE1 insertions ###
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

### build temporal files for checking pAT and blastx ###
my $sr_file    = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.tabe.discover';
my $temp1_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.pAT.input';
my $temp2_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.blastx.input';
my $cleavage_file = ($strand) ? $masterpath.'/refTE/position/TTAAAA.'.$hg.'.sorted.pA.bed' : $masterpath.'/refTE/position/TTAAAA.'.$hg.'.sorted.pT.bed';
my $pAT_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.pAT.output';
my $blastx_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.blastx.output';
my $db_file = $masterpath.'/refTE/blastx/L1HS.db';

(open (SR, "<$sr_file")) || die "cannot open the PE discover file\n";
(open (TEMP1, ">$temp1_file")) || die "cannot open the temp1 file\n";
(open (TEMP2, ">$temp2_file")) || die "cannot open the temp2 file\n";

my ($line, @data, $pos1, $pos2, $r1, @temp);
while ($line = <SR>)
   {
    chomp($line);
    @data = split("\t", $line);
    next if ( (($data[0] !~ /^[0-9XY]/) && ($data[0] !~ /^chr[0-9XY]/)) || ($data[0] =~ /\_/) || ($data[0] =~ /\-/) || ($data[3] !~ /L1/) );
    ($r1, @temp) = split(/\:/, $data[4]);
    $data[4] = join(':', @temp);
    $data[4] = $r1.':'.$data[4];
    print TEMP1 "$data[0]\t$data[1]\t$data[2]\t$data[4]\n";

    $pos1 = ($data[13] < $data[14]) ? $data[13] : $data[14];
    $pos2 = ($data[13] > $data[14]) ? $data[13] : $data[14];
    if ( (($pos1 >= 908) && ($pos2 <=1921)) || (($pos1 >= 1988) && ($pos2 <=5812)) )
       {
        print TEMP2 ">$data[4]\n$data[16]\n";
       }
   }
close (TEMP1);
close (TEMP2);
#system ("closestBed -d -t first -a $temp1_file -b $cleavage_file > $pAT_file");
system ("sort -T $outpath -k1,1 -k2,3n $temp1_file | closestBed -d -t first -a stdin -b $cleavage_file > $pAT_file");
system ("blastx -query $temp2_file -db $db_file -outfmt 6 -out $blastx_file");

### L1 SR reads data matricies ###
my $seg_file  = $outpath.'/'.$sub.'/'.$retro.'/seg.sr.scores';
my $sr_old    = $outpath.'/'.$sub.'/retro_v'.$ARGV[2].'/'.$sub.'.sr2.discover';
my $ou_file   = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.sr.LINE.matrix';
my $call_file = $outpath.'/'.$sub.'/'.$retro.'/LINE/'.$sub.'.LINE.novel.calls';
my $aca_file  = $outpath.'/'.$sub.'/fix_ACA/'.$sub.'.sr.discover';
my $acg_file  = $outpath.'/'.$sub.'/fix_ACG/'.$sub.'.sr.discover';
my $tag_file  = $outpath.'/'.$sub.'/fix_TAG/'.$sub.'.sr.discover';
my $aca_all   = $outpath.'/'.$sub.'/fix_ACA/'.$sub.'.sr.L1HS.discover';
my $tag_all   = $outpath.'/'.$sub.'/fix_TAG/'.$sub.'.sr.L1HS.discover';
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
(open (ACA , "<$aca_file ")) || die "cannot open the ACA file\n";
(open (ACG , "<$acg_file ")) || die "cannot open the ACG file\n";
(open (TAG , "<$tag_file ")) || die "cannot open the ACG file\n";
(open (ALL1, "<$aca_all  ")) || die "cannot open the aca all file\n";
(open (ALL2, "<$tag_all  ")) || die "cannot open the tag all file\n";
(open (DEP1, "<$dep1_file")) || die "cannot open the depth file\n";
(open (DEP2, "<$dep2_file")) || die "cannot open the depth file\n";
(open (SOLD, "<$sr_old   ")) || die "cannot open the sr old file $sr_old\n";
(open (PAT , "<$pAT_file ")) || die "cannot open the pat file\n";
(open (BLAST, "<$blastx_file")) || die "cannot open the blast file\n";

my (@reads, %aca, %acg, %tag, %uniq, %score);
my ($i, $j, $k, @temp1, @temp2, @temp3, @site, %seg, %call, %all1, %all2);
my ($call, $sr, $pe, $ref, $anchor_alignment, $anchor_length, $anchor_mm, $pos, $neg, $comp);
my ($acag, $tag, $read, $map, $count);
my ($depth, $len, $size, %cord, %depth);
my ($te, $short, $end3, $end5, $direction, $min, $max);
my ($med, $mad, @dep_all, $med_x, $mad_x, @Xdep_all, $med_y, $mad_y, @Ydep_all);
my ($start, $end, $cigar, $ratio, $center, $anchor, $insert);
my ($A_pair, $A_insert, $a_size, $a_map,  $A_mm, $A_MapQ, $A_AS, $A_XS, %old_sr, $old);
my (%orf, $ORF, $score, %pat);

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

while ($line = <ACA>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $aca{$r1.':'.$data[4]} = 1;
   }

while ($line = <ACG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $acg{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL1>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $all1{$r1.':'.$data[4]} = 1;
   }

while ($line = <TAG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $tag{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL2>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[5] & 0x40) ? "r1" : "r2";
    $all2{$r1.':'.$data[4]} = 1;
   }

print "6: read ACA ACG TAG calls\n";

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
$med_x= &median(@Xdep_all);
$mad_x = &mad(@Xdep_all);
$med_y = &median(@Ydep_all);
$mad_y = &mad(@Ydep_all);
print "7: read depth file median=$med\nMAD=$mad\nXmedian=$med_x\nXmad = $mad_x\nYmedian=$med_y\nYmad = $mad_y\n";

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

$count = 0;
print OU "chr\tcord1\tcord2\tread\trefpos1\trefpos2\tpos\tneg\tseg\tsupport\tsr\tpe\tmap\tref\tACAG\tTAG\tdepth\tmap_len\tmap_size\tmap_ratio\tshort\tend3\tend5\tdirection\tcenter\tanchor\tinsert\tA_pair\tA_insert\tA_mm\tA_MapQ\tA_AS\tA_XS\toldSR\tgap\tpAT\tdist\tORF\tXscore\trefpos\n";
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

    if ($data[6] =~ /(\d+)M\d+S$/)
      {
       $insert = 1;
       $start = $data[13];
       $end   = $data[14];
       $te    = $data[3];
       if ( ( ($te =~ 'Alu' && $start < 270) || ($te =~ 'L1' && $start < 5900) ) && ($start > $end) )
          {
           $end3 = 1;   ### no 3' ###
          }
       elsif ( ( ($te =~ 'Alu' && $end < 270) || ($te =~ 'L1' && $end < 5900) ) && ($start < $end) && ($data[12] < ($data[10] - 10) ) )
          {
           $end3 = 1; ### something inserted before the full length 3'end ###
          }
       elsif ( ($te =~ /L1/) && ( abs(($start + $end)/2-6064) < $CUT))
          {
           $short = 1; ### too short ###
          }
      }
    elsif ($data[6] =~ /^\d+S\d+M/)
      {
       $insert= -1;
       $start = $data[14];
       $end   = $data[13];
       $te    = $data[3];
       if ( ( ($te =~ 'Alu' && $start < 270) || ($te =~ 'L1' && $start < 5900) ) && ($start > $end) )
          {
           $end3 = 1;
          }
       elsif ( ( ($te =~ 'Alu' && $end < 270) || ($te =~ 'L1' && $end < 5900) ) && ($start < $end) && ($data[11] > 10 ) )
          {
           $end3 = 1;
          }
       elsif ( ($te =~ /L1/) && ( abs(($start + $end)/2-6064) < $CUT))
          {
           $short = 1; ### too short ###
          }
       }

    $end5 = 0;
    if ( ($data[13] > $data[14]) && ($data[14]>4) && ( (length($data[17]) - $data[12]) > 4) )
       {
        $end5 = 1;
       }
    elsif ( ($data[13] < $data[14]) && ($data[13]>4) && ($data[11] > 4) )
       {
        $end5 = 1;
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
    $size = length($data[17]);        ### sequence length ###
    $ratio= $len / $size;
    $center = abs(($data[13] + $data[14])/2 - 3000);
    #$refpos = ($data[13] + $data[14])/2;

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

    ### unique entries ###
    if (!exists($score{$data[4]}))
       {
        $uniq{$data[4]} = "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[13]\t$data[14]\t$pos\t$neg\t$comp\t$call\t$sr\t$pe\t$map\t$ref\t$acag\t$tag\t$depth\t$len\t$size\t$ratio\t$short\t$end3\t$end5\t$direction\t$center\t$anchor\t$insert\t$A_pair\t$A_insert\t$A_mm\t$A_MapQ\t$A_AS\t$A_XS\t$old\t$data[24]\t$pat{$data[4]}\t$ORF\t$score\t$pos1";
        $score{$data[4]} = $map * $len;
       }
    elsif (($map * $len) > $score{$data[4]})
       {
        $uniq{$data[4]} = "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[13]\t$data[14]\t$pos\t$neg\t$comp\t$call\t$sr\t$pe\t$map\t$ref\t$acag\t$tag\t$depth\t$len\t$size\t$ratio\t$short\t$end3\t$end5\t$direction\t$center\t$anchor\t$insert\t$A_pair\t$A_insert\t$A_mm\t$A_MapQ\t$A_AS\t$A_XS\t$old\t$data[24]\t$pat{$data[4]}\t$ORF\t$score\t$pos1";
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

