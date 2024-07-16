#!/usr/bin/env perl

use warnings;
use strict;
### Feature matrix for Alu PE reads ###

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

### look for Alu PE supporting reads with split-read anchors ###
my $in_file    = $outpath.'/'.$sub.'/'.$retro.'/ALU/'.$sub.'.ALU.novel.sites';
my $fasta_temp = $outpath.'/'.$sub.'/'.$retro.'/ALU/temp.fa';
my $exo_temp   = $outpath.'/'.$sub.'/'.$retro.'/ALU/temp.exo.txt';
my $seg_temp   = $outpath.'/'.$sub.'/'.$retro.'/ALU/temp.seg.txt';
my $sanch_file = $outpath.'/'.$sub.'/'.$retro.'/ALU/split-reads-anchor.txt';
my $ref_file   = $masterpath.'/refTE/sequence/ALU.fa';

my $LEN = 20;
(open (IN   , "<$in_file   ")) || die "cannot open the in file\n";
(open (FASTA, ">$fasta_temp")) || die "cannot open the fasta file\n";
(open (OUT  , ">$sanch_file")) || die "cannot open the out file\n";

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

### build temporal files for checking pAT ###
my $pe_file    = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.filter.discover';
my $temp1_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.pe.ALU.pAT.input';
my $cleavage_file = ($strand) ? $masterpath.'/refTE/position/TTAAAA.'.$hg.'.sorted.pA.bed' : $masterpath.'/refTE/position/TTAAAA.'.$hg.'.sorted.pT.bed';
my $pAT_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.pe.ALU.pAT.output';

(open (PE, "<$pe_file")) || die "cannot open the PE discover file\n";
(open (TEMP1, ">$temp1_file")) || die "cannot open the temp1 file\n";

my ($pos1, $pos2);
while ($line = <PE>)
   {
    chomp($line);
    @data = split("\t", $line);
    next if ( (($data[0] !~ /^[0-9XY]/) && ($data[0] !~ /^chr[0-9XY]/)) || ($data[0] =~ /\_/) || ($data[0] =~ /\-/) || ($data[3] !~ /Alu/) );
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $data[4] = $r1.':'.$data[4];
    print TEMP1 "$data[0]\t$data[1]\t$data[2]\t$data[4]\n";
   }
close (TEMP1);
#system ("closestBed -d -t first -a $temp1_file -b $cleavage_file > $pAT_file");
system ("sort -T $outpath -k1,1 -k2,3n $temp1_file | closestBed -d -t first -a stdin -b $cleavage_file > $pAT_file");

### Feature matrix for Alu PE supporting reads ###
my $seg_file  = $outpath.'/'.$sub.'/'.$retro.'/seg.pe.scores';
my $ou_file   = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.pe.ALU.matrix';
my $call_file = $outpath.'/'.$sub.'/'.$retro.'/ALU/'.$sub.'.ALU.novel.calls';
my $site_file = $outpath.'/'.$sub.'/'.$retro.'/ALU/'.$sub.'.ALU.novel.sites';
my $dep1_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.bed';
my $dep2_file = $outpath.'/'.$sub.'/'.$retro.'/'.$sub.'.depth';

my $ggcg_file = $outpath.'/'.$sub.'/fix_GGCG/'.$sub.'.discover';
my $ggcg_all = $outpath.'/'.$sub.'/fix_GGCG/'.$sub.'.Alu.discover';
my $gagc_file = $outpath.'/'.$sub.'/fix_GAGC/'.$sub.'.discover';
my $gagc_all = $outpath.'/'.$sub.'/fix_GAGC/'.$sub.'.Alu.discover';
my $cc_file = $outpath.'/'.$sub.'/fix_CC/'.$sub.'.discover';
my $cc_all = $outpath.'/'.$sub.'/fix_CC/'.$sub.'.Alu.discover';
my $tg_file = $outpath.'/'.$sub.'/fix_TG/'.$sub.'.discover';
my $tg_all = $outpath.'/'.$sub.'/fix_TG/'.$sub.'.Alu.discover';
my $cg_file = $outpath.'/'.$sub.'/fix_CG/'.$sub.'.discover';
my $cg_all = $outpath.'/'.$sub.'/fix_CG/'.$sub.'.Alu.discover';
my $at_file = $outpath.'/'.$sub.'/fix_AT/'.$sub.'.discover';
my $at_all = $outpath.'/'.$sub.'/fix_AT/'.$sub.'.Alu.discover';
my $Yb_file = $outpath.'/'.$sub.'/fix_AluYb/'.$sub.'.discover';
my $Yb_all = $outpath.'/'.$sub.'/fix_AluYb/'.$sub.'.Alu.discover';

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
(open (SPLIT, "<$sanch_file")) || die "cannot open the split file\n";
(open (PAT, "<$pAT_file")) || die "cannot open the split file\n";

my ($temp2_file);
my (@reads, %aca, %acg, %tag, %uniq, %score);
my ($i, $j, $k, @temp, @temp1, @temp2, @temp3, @site, %call);
my ($call, $sr, $pe, $ref, $anchor_alignment, $anchor_length, $anchor_mm, $pos, $neg, $comp, $refpos);
my ($acag, $tag, $a_size, $a_map, $map, $count);
my ($depth, $len, $size, %cord, %depth);
my ($short, $end3, $end5, $cigar, $trans2, $transduction, $pA_perc, $pA_length, $c1, $c2, $c3, $min, $max);
my ($med, $mad, @dep_all, $med_X, $mad_X, @Xdep_all, $med_Y, $mad_Y, @Ydep_all);
my (%ggcg, %gagc, %cc, %tg, %cg, %at, %Yb);
my ($up, $ggcg, $gagc, $cc, $tg, $cg, $at, $yb);
my (%all1, %all2, %all3, %all4, %all5, %all6, %all7);
my (%pat, %split, $anchor_split, $anchor_direct, $anchor_insert, $anchor_seg, $anchor_map, $anchor_XS, $anchor_AS, $anchor_mapQ);
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

while ($line = <GGCG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $ggcg{$r1.':'.$data[4]} = 1;
   }

while ($line = <GAGC>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $gagc{$r1.':'.$data[4]} = 1;
   }

while ($line = <CC>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $cc{$r1.':'.$data[4]} = 1;
   }

while ($line = <TG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $tg{$r1.':'.$data[4]} = 1;
   }

while ($line = <CG>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $cg{$r1.':'.$data[4]} = 1;
   }

while ($line = <AT>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $at{$r1.':'.$data[4]} = 1;
   }

while ($line = <Yb>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $Yb{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL1>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all1{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL2>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all2{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL3>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all3{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL6>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all6{$r1.':'.$data[4]} = 1;
   }

while ($line = <ALL7>)
   {
    chomp($line);
    @data = split("\t", $line);
    $r1 = ($data[14] & 0x40) ? "r1" : "r2";
    $all7{$r1.':'.$data[4]} = 1;
   }
print "6: read alleles calls\n";

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
$med = &median(@dep_all);
$mad = &mad(@dep_all);
$med_X = &median(@Xdep_all);
$mad_X = &mad(@Xdep_all);
$med_Y = &median(@Ydep_all);
$mad_Y = &mad(@Ydep_all);
print "7: read depth file median=$med\nMAD=$mad\nXmedian=$med_X\nXmad = $mad_X\nYmedian=$med_Y\nYmad = $mad_Y\n";

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
print "8: read the splitanchor and pAT file\n";

$count =1;
print OU "chr\tcord1\tcord2\tread\trefpos1\trefpos2\tpos\tneg\tseg\tsupport\tsr\tpe\tmap\tanchor_align\tanchor_len\tanchor_mm\tanchor_split\tanchor_direct\tanchor_insert\tanchor_seg\tanchor_map\tanchor_XS\tanchor_AS\tanchor_mapQ\tref\tGGCG\tGAGC\tCC\tTG\tCG\tAT\tYB\tdepth\tmap_len\tshort\tend3\tend5\ttransduction\ttrans2\tpA_perc\tpA_length\tc1\tc2\tc3\tdirection\trefpos\tupstream\tgap\tpAT\tdist\n";
while ($line = <PE>)
   {
    chomp($line);
    #print "read $count\n";
    $count ++;
    @data = split("\t", $line);
    next if ( (($data[0] !~ /^[0-9XY]/) && ($data[0] !~ /^chr[0-9XY]/)) || ($data[0] =~ /\_/) || ($data[0] =~ /\-/) || ($data[3] !~ /Alu/) );
    $read = $data[4].'###'.$data[3].'###'.$data[7].'###'.$data[8].'###'.$data[9];
    $r1   = ($data[14] & 0x40) ? "r1" : "r2";
    $data[4] = $r1.':'.$data[4];

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

    $c1 = ($data[5] eq '+') ? 1 : -1;
    $c2 = ($data[6] eq '+') ? 1 : -1;
    $c3 = ($data[15] eq '+') ? 1 : -1;
    $direction = -1 * $c1 * $c2 * $c3;
    $up = $c1 * $direction;
    $short = 0;
    $refpos= ($data[10] + $data[11])/2;
    $min = ($data[11] < $data[10]) ? $data[11] : $data[10];
    $max = ($data[11] > $data[10]) ? $data[11] : $data[10];
    $anchor_XS = $data[22];
    $anchor_AS = $data[23];
    $anchor_mapQ = $data[24];

    $end3 = &parse_end3($data[8], $data[9], $data[10], $data[11], $data[12]);
    $end5 = &parse_end5($data[8], $data[9], $data[10], $data[11], $data[12], $data[13]);
    $transduction = &parse_trans($data[8], $data[9], $data[10], $data[11], $data[12], $data[13]);
    $trans2 = &parse_trans2($data[8], $data[9], $data[10], $data[11], $data[12], $data[13]);
    ($pA_perc, $pA_length) = &parse_polyA($data[8], $data[9], $data[10], $data[11], $data[12]);

    $pos = 1; #exists($pos{$data[4]}) ? 1 : 0;
    $neg = 1; #exists($neg{$data[4]}) ? 1 : 0;
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
        $depth = ( $depth{$cord{$data[4]}} - $med_X ) / $mad_X;
       }
    elsif ($data[0] =~ /Y/)
       {
        $depth = ( $depth{$cord{$data[4]}} - $med_Y ) / $mad_Y;
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
    if (!exists($score{$data[4]}))
       {
        $uniq{$data[4]} = "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[10]\t$data[11]\t$pos\t$neg\t$comp\t$call\t$sr\t$pe\t$map\t$anchor_alignment\t$anchor_length\t$anchor_mm\t$anchor_split\t$anchor_direct\t$anchor_insert\t$anchor_seg\t$anchor_map\t$anchor_XS\t$anchor_AS\t$anchor_mapQ\t$ref\t$ggcg\t$gagc\t$cc\t$tg\t$cg\t$at\t$yb\t$depth\t$len\t$short\t$end3\t$end5\t$transduction\t$trans2\t$pA_perc\t$pA_length\t$c1\t$c2\t$c3\t$direction\t$refpos\t$up\t$data[25]\t$pat{$data[4]}";
        $score{$data[4]} = $map;
       }
    elsif ($map > $score{$data[4]})
       {
        $uniq{$data[4]} = "$data[0]\t$data[1]\t$data[2]\t$data[4]\t$data[10]\t$data[11]\t$pos\t$neg\t$comp\t$call\t$sr\t$pe\t$map\t$anchor_alignment\t$anchor_length\t$anchor_mm\t$anchor_split\t$anchor_direct\t$anchor_insert\t$anchor_seg\t$anchor_map\t$anchor_XS\t$anchor_AS\t$anchor_mapQ\t$ref\t$ggcg\t$gagc\t$cc\t$tg\t$cg\t$at\t$yb\t$depth\t$len\t$short\t$end3\t$end5\t$transduction\t$trans2\t$pA_perc\t$pA_length\t$c1\t$c2\t$c3\t$direction\t$refpos\t$up\t$data[25]\t$pat{$data[4]}";
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

sub parse_end3
   {
    my $p1 = shift;
    my $p2 = shift;
    my $p3 = shift;
    my $p4 = shift;
    my $seq = shift;
    my $end3 = 0;

    if ( ($p3 < $p4) && ($p4 < 280) && ((length($seq) - $p2) > 5 ))
       {
        $end3 = 1;
       }
    elsif ( ($p3 > $p4) && ($p3 < 280) && ($p1 > 5) )
       {
        $end3 = 1;
       }
    return($end3);
   }

sub parse_end5
   {
    my $p1 = shift;
    my $p2 = shift;
    my $p3 = shift;
    my $p4 = shift;
    my $seq = shift;
    my $cigar=shift;
    my $end5 = 0;

    if ( ($p3 > $p4) && ($p4 > 4) && ( (length($seq) - $p2) > 4 ))
       { ### Alu is missing 5'end and there is extra sequence on 5'end ###
        if ($cigar =~ /(\d+)S$/)
           {
            $end5 = 1 if ( abs((length($seq) - $p2) - $1) > 4 );
           }
        else
           {
            $end5 = 1;
           }
       }
    elsif ( ($p3 < $p4) && ($p3 > 10) && ($p1 > 10) )
       { ### Alu is missing 5'end and there is extra sequence on 5 end ###
        if ($cigar =~ /^(\d+)S/)
           {
            $end5 = 1 if ( abs($p1 - $1) > 4 );
           }
        else
           {
            $end5 = 1;
           }
       }
    return($end5);
   }
 
sub parse_trans
   {
    my $p1 = shift;
    my $p2 = shift;
    my $p3 = shift;
    my $p4 = shift;
    my $seq= shift;
    my $cigar = shift;
    my $trans = 0;
    my ($max, $start, $end, $substr);

    if ( ($cigar !~ /\d+S$/) && ($cigar !~ /^\d+S/) )
       {
        if ( $p1 > 5 )
           {
            $substr = substr($seq, 0, $p1);
           }
        elsif ((length($seq) - $p2) > 5)
           {
            $substr = substr($seq, $p2, (length($seq) - $p2));
           }
        if ( ( ($p1 > 5) || ((length($seq) - $p2) > 5) ) && ($p3 < 280) && ($p4 < 280) )
           {
            $trans = 1;
           }
        elsif ( ($p1 > 5) || ((length($seq) - $p2) > 5) )
           {
            $max = 0;
            if ($p3 < $p4)
               {
                while ($substr =~ /(A+)/g)
                   {
                    $start = $-[0];
                    $end   = $+[0];
                    $max  += $end - $start;
                   }
               }
            else
               {
                while ($substr =~ /(T+)/g)
                   {
                    $start = $-[0];
                    $end   = $+[0];
                    $max  += $end - $start;
                   }
               }
            $trans = 1 if ( (length($substr) - $max) > 4);
           }
        }
    return($trans);
   }

sub parse_polyA
   {
    my $p1 = shift;
    my $p2 = shift;
    my $p3 = shift;
    my $p4 = shift;
    my $seq = shift;
    my ($substr, $start, $end);
    my $polyA = 'NA';
    my $max_pA = 'NA';

    if ( $p1 > 5 )
       {
        $substr = substr($seq, 0, $p1);
       }
    elsif ((length($seq) - $p2) > 5)
       {
        $substr = substr($seq, $p2, (length($seq) - $p2));
       }
    if ( (($p1 > 5) || ((length($seq) - $p2) > 5)) && ($p3 > 280 || $p4 > 280) )
       {
        $max_pA = 0;
        if ($p3 < $p4)
           {
            while ($substr =~ /AAAA(A+)/g)
               {
                $start = $-[0];
                $end   = $+[0];
                $max_pA= (($end - $start) > $max_pA) ? ($end - $start) : $max_pA;
               }
           }
        else
           {
            while ($substr =~ /TTTT(T+)/g)
               {
                $start = $-[0];
                $end   = $+[0];
                $max_pA= (($end - $start) > $max_pA) ? ($end - $start) : $max_pA;
               }
           }
         $polyA = $max_pA / length($substr);
         $polyA = substr($polyA, 0, 5);
        }
     return ($polyA, $max_pA);
    }

sub parse_trans2
   {
    my $p1 = shift;
    my $p2 = shift;
    my $p3 = shift;
    my $p4 = shift;
    my $seq= shift;
    my $cigar = shift;
    my ($max, $start, $end, $substr, $substr1, $substr2);

    if ( ($cigar !~ /\d+S$/) && ($cigar !~ /^\d+S/) )
       {
        if ( $p1 > 5 )
           {
            $substr1 = substr($seq, 0, $p1);
           }
        if ((length($seq) - $p2) > 5)
           {
            $substr2 = substr($seq, $p2, (length($seq) - $p2));
           }
        if ( ( ($p1 > 5) || ((length($seq) - $p2) > 5) ) && ($p3 < 280) && ($p4 < 280) )
           {
            return(1);
           }
        if ($p1 > 5)
           {
            $max = 0;
            $substr= $substr1;
            if ($p3 < $p4)
               {
                while ($substr =~ /(A+)/g)
                   {
                    $start = $-[0];
                    $end   = $+[0];
                    $max  += $end - $start;
                   }
               }
            else
               {
                while ($substr =~ /(T+)/g)
                   {
                    $start = $-[0];
                    $end   = $+[0];
                    $max  += $end - $start;
                   }
               }
            return(1) if ( (length($substr) - $max) > 4);
           }
        if ((length($seq) - $p2) > 5)
           {
            $substr=$substr2;
            $max = 0;
            if ($p3 < $p4)
               {
                while ($substr =~ /(A+)/g)
                   {
                    $start = $-[0];
                    $end   = $+[0];
                    $max  += $end - $start;
                   }
               }
            else
               {
                while ($substr =~ /(T+)/g)
                   {
                    $start = $-[0];
                    $end   = $+[0];
                    $max  += $end - $start;
                   }
               }
            return(1) if ( (length($substr) - $max) > 4);
           }
        }
    return(0);
   }
