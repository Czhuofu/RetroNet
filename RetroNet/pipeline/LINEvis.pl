#!/usr/bin/env perl

use warnings;
use strict;

use GD;
use GD::Arrow;

my $scale = 1;
my $edge  = 8;
my $masterpath= $ARGV[0];
my $in_file = $ARGV[9].'/04_PE.alignment';
my $sr_file = $ARGV[9].'/05_SR.alignment';
my $chr     = $ARGV[2];
my $te      = $ARGV[3];
my $hg      = $ARGV[5];
my $strand  = $ARGV[6];
my $map_file = $ARGV[10];
my $mat_file = $ARGV[7];
my $CORD1    = $ARGV[8]; ## coordiantes of the insertion ###
(open (MAP, "<$map_file")) || die "cannot open the map file\n $!";

my $refTE = $masterpath.'/refTE/sequence/'.$te.'.fa';
my $min_hg = 1000000000;
my $max_hg = 0;
my (@data, @map, @name, @reads, @redcolors, %transduct);
my ($im, $ou_file, $white, $blue, $black, $red, $green, $magenta);
my ($name, $i, $j, $k, $x, $y, $line, $width);
my ($x1, $x2, $y1, $y2, $col);
my ($seq1, $seq2, $seq3, @colon, %anchor, %anchor_split, %split_read, %unique);
my $TE_length;
my (@data1, @data2, @name1, @name2);

if (($te =~ /^L1/) || ($te eq 'LINE'))
   {
    $TE_length = 6300;
   }
elsif ($te =~ /Alu/)
   {
    $TE_length = 400;
    (open (MAT, "<$mat_file")) || die "cannot opne the mat file\n";
    <MAT>;
    while ($line = <MAT>)
       {
        chomp($line);
        @data = split("\t", $line);
        if ($data[38])
           {
            $data[3] =~ s/^r[12]\://;
            $transduct{$data[3]} = 1;
           }
       }
   }

$i = 0;
if ( (-e $in_file) && (! (-z $in_file)) )
   {
    (open (ALN, "<$in_file")) || die "cannot open the in file\n";
    
    while ($line = <ALN>)
       {
        chomp($line);
        if ($line =~ /Target range/)
           {
            <ALN>;
            chomp($line = <ALN>);
            $seq1 = $seq2 = $seq3 = '';
            while ( $line !~ /vulgar/ )
               {
                @colon = ();
                while ($line  =~ /\:/g)
                   {
                    push @colon, $-[0];
                   }
                $seq1 .= substr($line, $colon[0]+2, $colon[1]-$colon[0]-3);
                chomp($line = <ALN>);
                $seq2 .= substr($line, $colon[0]+2, $colon[1]-$colon[0]-3);
                chomp($line = <ALN>);
                $seq3 .= substr($line, $colon[0]+2, $colon[1]-$colon[0]-3);
                <ALN>;
                chomp($line = <ALN>);
               }
            if ($line =~ /vulgar/)
               {
                chomp($line = <ALN>);
                @data = split(" ", $line);
                next if ($data[5] ne $te);
                @name = split(';', $data[1]);
                $min_hg = ($min_hg < $name[2]) ? $min_hg : $name[2];
                $max_hg = ($max_hg > $name[3]) ? $max_hg : $name[3];
                next if (exists($unique{$name[0]}));
                $unique{$name[0]} = 1;
                if ( ($name[0] !~ /anchor/ ) && ($data[8] < $data[9]) )
                   {
                    $reads[$i][0] = $name[0];
                    $reads[$i][1] = $line;
                    $reads[$i][2] = $seq1;
                    $reads[$i][3] = $seq2;
                    $reads[$i][4] = $seq3;
                    $reads[$i][5] = $data[8];
                    $reads[$i][6] = $name[2];
                    $i ++;
                   }
                elsif ( ($name[0] !~ /anchor/ ) && ($data[8] > $data[9]) )
                   {
                    $reads[$i][0] = $name[0];
                    $reads[$i][1] = $line;
                    $reads[$i][2] = &rc($seq1);
                    $reads[$i][3] = reverse($seq2);
                    $reads[$i][4] = &rc($seq3);
                    $reads[$i][5] = $data[9];
                    $reads[$i][6] = $name[2];
                    $i++;
                   }
                if ( ($name[0] =~ /anchor/ ) && ($data[8] < $data[9]) )
                   { ### supporting read anchor end is a split read ###
                    $anchor_split{$name[0]} = $name[4];  ### anchor_direct ###
                    $anchor{$name[0]}[0] = $line;
                    $anchor{$name[0]}[1] = $seq1;
                    $anchor{$name[0]}[2] = $seq2;
                    $anchor{$name[0]}[3] = $seq3;
                    $anchor{$name[0]}[4] = $data[8];  ### ME reference coord ###
                    $anchor{$name[0]}[5] = $name[5]."\t".$name[6];  ### anchor_length and cigar ###
                   }
                elsif ( ($name[0] =~ /anchor/ ) && ($data[8] > $data[9]) )
                   {
                    $anchor_split{$name[0]} = $name[4];
                    $anchor{$name[0]}[0] = $line;
                    $anchor{$name[0]}[1] = &rc($seq1);
                    $anchor{$name[0]}[2] = reverse($seq2);
                    $anchor{$name[0]}[3] = &rc($seq3);
                    $anchor{$name[0]}[4] = $data[9];
                    $anchor{$name[0]}[5] = $name[5]."\t".$name[6];
                   }
               }
           }
       }
   close (ALN);
   }

if ( (-e $sr_file) && (! (-z $sr_file)) )
   {
    (open (ALN, "<$sr_file")) || die "cannot open the sr file $sr_file\n";
    while ($line = <ALN>)
       {
        chomp($line);
        if ($line =~ /Target range/)
          {
           <ALN>;
           chomp($line = <ALN>);
           $seq1 = $seq2 = $seq3 = '';
           while ( $line !~ /vulgar/ )
             {
              @colon = ();
              while ($line  =~ /\:/g)
                     {
                    push @colon, $-[0];
                   }
              $seq1 .= substr($line, $colon[0]+2, $colon[1]-$colon[0]-3);
              chomp($line = <ALN>);
              $seq2 .= substr($line, $colon[0]+2, $colon[1]-$colon[0]-3);
              chomp($line = <ALN>);
              $seq3 .= substr($line, $colon[0]+2, $colon[1]-$colon[0]-3);
              <ALN>;
              chomp($line = <ALN>);
             }
           if ($line =~ /vulgar/)
            {
             chomp($line = <ALN>);
             @data = split(" ", $line);
             @name = split(';', $data[1]);
             next if (exists($split_read{$name[0]}));
             next if ( ($name[1] ne $chr) || (abs($name[2] - $CORD1) > 5000) );
             $min_hg = ($min_hg < $name[2]-150) ? $min_hg : ($name[2]-150);
             $max_hg = ($max_hg > $name[3]+150) ? $max_hg : ($name[3]+150);
             $split_read{$name[0]} = 1;
             if ($data[8] < $data[9])
                {
                 $reads[$i][0] = $name[0];
                 $reads[$i][1] = $line;
                 $reads[$i][2] = $seq1;
                 $reads[$i][3] = $seq2;
                 $reads[$i][4] = $seq3;
                 $reads[$i][5] = $data[8];
                 $reads[$i][6] = $name[2];
                }
             else
                {
                 $reads[$i][0] = $name[0];
                 $reads[$i][1] = $line;
                 $reads[$i][2] = &rc($seq1);
                 $reads[$i][3] = reverse($seq2);
                 $reads[$i][4] = &rc($seq3);
                 $reads[$i][5] = $data[9];
                 $reads[$i][6] = $name[2];
                }
             $i++;
           }
         }
     }
   close (ALN);
   }
my $flank = 450;
$min_hg = ($min_hg > $flank ) ? $min_hg-$flank : 0;
$max_hg += $flank;
#print "min=$min_hg\tmax=$max_hg\n";

$i = 0;
while ($line = <MAP>)
    {
     chomp($line);
     @data = split("\t", $line);
     next if ($data[0] ne $chr);
     $map[$i][0] = $data[1];
     $map[$i][1] = $data[2];
     $map[$i][2] = $data[3];
     $i ++;
     last if ($data[2] > $max_hg);
    }

my %ks;
my $count = 10;
my $count_i = 0;
my $size = $#reads + 1;
# print "$size\n";
$j = $k = -1;

while ( ($count_i <= $count) && ($count_i < ($size * ($size-1))/2) )
   {
    while (1)
        {
         $j = int(rand $size);
         $k = int(rand $size); 
#next if (!( (($j == 8) || ($j == 22)) && (($k == 8) || ($k == 22)) ));
#next if (($j != 8) && ($j != 22));
#print "$j\t$reads[$j][0]\n";
         if (($j != $k) && (!exists($ks{$j."\t".$k})))
            {
             $ks{$j."\t".$k} = 1;
             $ks{$k."\t".$j} = 1;
             last;
            }
        }
    $count_i++;

        $ou_file = $ARGV[4].'/j'.$j.'_k'.$k.'_'.$ARGV[1];
        (open (OU , ">$ou_file ")) || die "cannot open the ou file\n";
        $width = ($max_hg - $min_hg +300 > $TE_length) ? $max_hg - $min_hg + 300 : $TE_length;
        $x = $edge*20 + $width*$scale + $edge*20;
        $y = $edge*2 + 5*$scale + $edge + 5*$scale + $edge + 5*$scale + $edge + 5*$scale;
        $im = GD::Image->new($x,$y);

        @data1 = split(" ", $reads[$j][1]);
        @name1 = split(';', $data1[1]);
        @data2 = split(" ", $reads[$k][1]);
        @name2 = split(';', $data2[1]);
        $min_hg = ($name1[2] < $name2[2]) ? $name1[2] : $name2[2];
        $max_hg = ($name1[3] > $name2[3]) ? $name1[3] : $name2[3];   
        $min_hg = ($min_hg > $flank ) ? $min_hg-$flank : 0;
        $max_hg += $flank;

        $white = $im->colorAllocate(255,255,255);
        $black = $im->colorAllocate(0,0,0);
        $red   = $im->colorAllocate(255,0,0);
        $blue  = $im->colorAllocate(0,0, 255);
        $green = $im->colorAllocate(0,255,0);
        $magenta= $im->colorAllocate(255,0,255);
        if ( ($strand && ($reads[$j][6] < $reads[$k][6])) || (!$strand && ($reads[$j][6] > $reads[$k][6])) )
           {
            &plot_str($reads[$j][0], $reads[$j][1], $reads[$j][2], $reads[$j][3], $reads[$j][4], $reads[$j][5], $reads[$k][0],$reads[$k][1], $reads[$k][2], $reads[$k][3], $reads[$k][4], $reads[$k][5]);
           }
        else
           {
            &plot_str($reads[$k][0], $reads[$k][1], $reads[$k][2], $reads[$k][3], $reads[$k][4], $reads[$k][5], $reads[$j][0], $reads[$j][1], $reads[$j][2], $reads[$j][3], $reads[$j][4], $reads[$j][5]);
           }
        binmode OU;
        print OU $im->png;
        close OU;
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
      '-' => '-',
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

### plot one pair ###
sub plot_str
  {
   my $n1 = shift;
   my $l1 = shift;
   my $s1 = shift;
   my $s2 = shift;
   my $s3 = shift;
   my $c1 = shift;

   my $n2 = shift;
   my $l2 = shift;
   my $s4 = shift;
   my $s5 = shift;
   my $s6 = shift;
   my $c2 = shift;

   my ($x1, $y1, $x2, $y2);
   my ($i, $line, $seq, @data);
   my ($row_read1, $row_read2, $row_anchor1, $row_anchor2);
  
   (open (REF, "<$refTE")) || die "cannot open the ref file\n";
   <REF>;
   $seq = '';
   while ($line = <REF>)
      {
       chomp($line);
       $seq .= $line;
      }

    @data = split("", $seq);
    for ($i=0; $i<=$#data; $i++)
       {
        if (($data[$i] eq 'A') || ($data[$i] eq 'a'))
          {
           $y1 = $edge*2 + $scale;
          }
        elsif ($data[$i] eq 'C' || ($data[$i] eq 'c'))
          {
           $y1 = $edge*2 + $scale * 2;
          }
        elsif ($data[$i] eq 'T' || ($data[$i] eq 't'))
          {
           $y1 = $edge*2 + $scale * 3;
          }
        elsif ($data[$i] eq 'G' || ($data[$i] eq 'g'))
          {
           $y1 = $edge*2 + $scale * 4;
          }
        else
          {
           next;
          }
        $x1 = $edge*20 + $i * $scale;
        $im->setPixel($x1,$y1,$black);
        $y1 += $edge + 5 * $scale;
        $y2 = $y1 + $scale;
        $im->setPixel($x1,$y1,$black);
        $y1 += $edge + 5 * $scale;
        $y2 = $y1 + $scale;
        $im->setPixel($x1,$y1,$black);
        $y1 += $edge + 5 * $scale;
        $y2 = $y1 + $scale;
        $im->setPixel($x1,$y1,$black);
       }
   ### plot mappability ### 
   ### adjust to the coordinates of read_i and read_j ###
   for ($i=0; $i<=$#map; $i++)
      {
       if ( ($map[$i][0] < $min_hg) && ($map[$i][1] >= $min_hg) )
         {
         if ($max_hg > $map[$i][1])
            {
             $x1 = ($strand) ? $edge : ($edge + ($max_hg - $map[$i][1]) * $scale);
             $x2 = ($strand) ? ($x1 + ($map[$i][1] - $min_hg)*$scale) : ($edge + ($max_hg-$min_hg) * $scale);
            }
         else
            {
             $x1 = $edge;
             $x2 = $edge + ($max_hg - $min_hg)*$scale;
            }

         $y1 = (1 - $map[$i][2]) * ($edge-1);
         $y2 = $edge-1;
         $im -> filledRectangle($x1,$y1,$x2,$y2,$black); 
         last if ($max_hg < $map[$i][1]);
         while ( $i++ && ($i <= $#map) )
            {
              if ($map[$i][1] <  $max_hg)
                 {
                  $x1 = ($strand) ? $edge + ($map[$i][0] - $min_hg)*$scale : $edge + ($max_hg - $map[$i][1])*$scale;
                  $x2 = ($strand) ? $edge + ($map[$i][1] - $min_hg)*$scale : $edge + ($max_hg - $map[$i][0])*$scale;
                 }
              else
                 {
                  $x1 = ($strand) ? $edge + ($map[$i][0] - $min_hg)*$scale : $edge;
                  $x2 = ($strand) ? ($edge + ($max_hg - $min_hg)*$scale) : ($edge + ($max_hg - $map[$i][0]) * $scale); 
                 }
              # print "x1=$x1 x2=$x2 y1=$y1 y2=$y2 col=$col\n";
              $y1 = (1 - $map[$i][2]) * ($edge-1);
              $y2 = $edge-1;
              $im -> filledRectangle($x1,$y1,$x2,$y2,$black); 
              last if ($map[$i][1] > $max_hg);
            }
          last;
         }
      }
    ### order of the 4 rows ###
    ### read pair 1 ###
    @data = split(" ", $l1);
    @name = split(';', $data[1]);
    if ( ($strand  && ($name[4] eq '+')) || ((!$strand) && ($name[4] eq '-')) )
       { ### anchor row1 : read row 2 ###
        $row_read1   = $edge * 3 + $scale * 5; ### row 2 ###
        $row_anchor1 = $edge; ### row 1 ###
       }
    else
       { ### anchor row2: read row 1 ###
        $row_read1 = $edge*2;  ### row 1 ###
        $row_anchor1 = (exists($split_read{$name[0]})) ? $edge : $edge*2 + $scale * 5; ### row 2 ###
       }

    ### read pair 2 ###
    @data = split(" ", $l2);
    @name = split(';', $data[1]);
    if ( ($strand && ($name[4] eq '+')) || ((!$strand) && ($name[4] eq '-')) )
       { ### anchor row3 : read row 4 ###
        $row_read2   = $edge * 5 + $scale * 15;
        $row_anchor2 = $edge * 3 + $scale * 10;
       }
    else
       { ### anchor row 4: read row 3 ###
        $row_read2 = $edge * 4 + $scale * 10;
        $row_anchor2 = (exists($split_read{$name[0]})) ?  $edge * 3 + $scale * 10 : $edge * 4 + $scale * 15;
       }

    ### draw read 1 ###
    # print "load MEend sub for $n1\n";
    &draw_MEend($n1, $l1, $s1, $s2, $s3, $c1, $row_read1);
    &draw_MEend($n2, $l2, $s4, $s5, $s6, $c2, $row_read2);

    ### draw anchor 1 ###
    &draw_anchor($l1, $row_anchor1);
    &draw_anchor($l2, $row_anchor2);
   }

sub draw_anchor
   {
    my $l1 = shift;
    my $row_anchor1 = shift;
     
    my @data = split(" ", $l1);
    my @name = split(';', $data[1]);
    my ($cigar, $direction, $anchor_len, $arrow);
    my ($ins, $min, $max, $x1, $x2, $x3, $x4, $y1, $y2, $y3, $y4, $gap);
    my (@data1, @data2, @data3);
    my $n1 = $name[0].'_anchor';

    if (!exists($split_read{$name[0]}) && ($name[7] =~ /S/) )
       {
        if ($name[7] =~ /^(\d+)S/)
           {
            if ($strand)
               {
                $x1 = $edge + ($name[2] - $min_hg - $1) * $scale;
                $x2 = $edge + ($name[2] - $min_hg - 1) * $scale;
               }
            else
               {
                $x1 = $edge + ($max_hg - $name[2]) * $scale;
                $x2 = $edge + ($max_hg - $name[2] + $1 - 1) * $scale;
               }
            $y1 = $row_anchor1 + $edge / 2;
            $y2 = $y1;
            $im -> line($x1, $y1, $x2, $y2, $black);
           }
       if ($name[7] =~ /(\d+)S$/)
           {
            if ($strand)
               {
                $x1 = $edge + ($name[3] - $min_hg) * $scale;
                $x2 = $edge + ($name[3] - $min_hg + $1 - 1) * $scale;
               }
            else
               {
                $x1 = $edge + ($max_hg - $name[3] - $1) * $scale;
                $x2 = $edge + ($max_hg - $name[3] - 1) * $scale;
               }
            $y1 = $row_anchor1 + $edge / 2;
            $y2 = $y1;
            $im -> line($x1, $y1, $x2, $y2, $black);
           }
       }

    if (exists($split_read{$name[0]}))
       { ### plot the split read suppot read ###
        if ( $name[5] =~ /^(\d+)S(\d+)M/ )
           {
            $name[3] = $name[2] + $2;
            $name[4] = '-';
            $ins = $1;
            if ($strand)
              {
               $x1 = $edge + ($name[2] - $min_hg - $ins + $data[6]) * $scale;
               $x2 = $edge + ($name[2] - $min_hg - $ins + $data[7] - 1) * $scale;
               $x3 = $edge + ($name[2] - $min_hg - $ins) * $scale;
               $x4 = $edge + ($name[2] - $min_hg - 1) * $scale;
              }
            else
              {
               $x1 = $edge + ($max_hg - $name[2] + $ins - $data[7]) * $scale;
               $x2 = $edge + ($max_hg - $name[2] + $ins - $data[6] - 1) * $scale;
               $x3 = $edge + ($max_hg - $name[2]) * $scale;
               $x4 = $edge + ($max_hg - $name[2] + $ins - 1) * $scale;
              }
           }
        if ( $name[5] =~ /(\d+)M(\d+)S$/ )
           {
            $name[2] = $name[2] - $1;
            $name[4] = '+';
            $ins = $2;
            if ($strand)
               {
                $x1 = $edge + ($name[3] - $min_hg + $data[6]) * $scale;
                $x2 = $edge + ($name[3] - $min_hg + $data[7] - 1) * $scale;
                $x3 = $edge + ($name[3] - $min_hg) * $scale;
                $x4 = $edge + ($name[3] - $min_hg + $ins - 1) * $scale;
               }
            else
               {
                $x1 = $edge + ($max_hg - $name[3] - $data[7]) * $scale;
                $x2 = $edge + ($max_hg - $name[3] - $data[6] - 1) * $scale;
                $x3 = $edge + ($max_hg - $name[3] - $ins) * $scale;
                $x4 = $edge + ($max_hg - $name[3] - 1) * $scale;
               }
           }
        $y3 = $row_anchor1 + $edge/2;
        $y4 = $y3;
        $im -> line($x1,$y3+1,$x2,$y4+1,$red); 
        $im -> line($x3,$y3,$x4,$y4,$black); 
       }
    
    ### plot the PE anchor end ###
    if ( $strand && ($name[4] eq '+')) 
       {
        ### upstream flank for a +strand insertion ###
        $x1 = $edge + ($name[3] - $min_hg - 1)*$scale;
        $y1 = $row_anchor1 + $edge/2;
        $x2 = $edge + ($name[2] - $min_hg)*$scale;
        $y2 = $y1;
       } 
    elsif ( $strand && ($name[4] eq '-') ) 
       {
        ### downstream read ###
        $x1 = $edge + ($name[2] - $min_hg)*$scale;
        $y1 = $row_anchor1 + $edge/2;
        $x2 = $edge + ($name[3] - $min_hg - 1)*$scale;
        $y2 = $y1;
       }
    elsif ( (!$strand) && ($name[4] eq '+') )
       {
        ### upstream flank for a +strand insertion ###
        $x2 = $edge + ($max_hg - $name[2] - 1)*$scale;
        $y1 = $row_anchor1 + $edge/2;
        $x1 = $edge + ($max_hg - $name[3])*$scale;
        $y2 = $y1;
       }
    elsif ( (!$strand) && ($name[4] eq '-') ) 
       {
        ### downstream read ###
        $x2 = $edge + ($max_hg - $name[3])*$scale;
        $y1 = $row_anchor1 + $edge/2;
        $x1 = $edge + ($max_hg - $name[2] - 1)*$scale;
        $y2 = $y1;
       }
    # print "$x1\t$x2\t$y1\t$y2\n";
    $arrow = GD::Arrow::Full->new( 
                  -X1    => $x1, 
                  -Y1    => $y1, 
                  -X2    => $x2, 
                  -Y2    => $y2, 
                  -WIDTH => $edge/8,
              );
    $im->filledPolygon($arrow,$blue);
    if (exists($anchor_split{$n1}))
       {
        ### PE anchor_end is also a split read ###
        ($anchor_len, $cigar) = split("\t", $anchor{$n1}[5]);
        @data = split(" ", $anchor{$n1}[0]);
        if ($data[8] < $data[9])
           {
            $direction = ($anchor_split{$n1}) ? 1 : 0;
           }
        else
           {
            $direction = ($anchor_split{$n1}) ? 0 : 1;
           }

        $y3 = $row_anchor1 + $edge/2;
        $y4 = $y3;
        $min = ($x1 < $x2) ? $x1 : $x2;
        $max = ($x1 > $x2) ? $x1 : $x2;
        if ( $strand && ($cigar =~ /^(\d+)S/) )
           {
            $x3 = $min - $anchor_len*$scale;
            $x4 = $x3 + ($anchor_len-1)*$scale;
            $im -> line($x3,$y3,$x4,$y4,$black);
            $x3 = $min - ($anchor_len - $data[6]) * $scale;
            $x4 = $min - ($anchor_len - $data[7] + 1) * $scale;
            $im -> line($x3,$y3+1,$x4,$y4+1,$red);    
           }
        elsif ( !$strand && ($cigar =~ /^(\d+)S/) )
           {
            $x3 = $max+1;
            $x4 = $x3 + ($anchor_len-1) * $scale;
            $im -> line($x3,$y3,$x4,$y4,$black);
            $x3 = $max + 1 + ($anchor_len - $data[7]) * $scale;
            $x4 = $max + 1 + ($anchor_len - $data[6] - 1) * $scale;
            $im -> line($x3,$y3+1,$x4,$y4+1,$red);    
           }

        if ( $strand && ($cigar =~ /(\d+)S$/) )
           {
            $x3 = $max+1;
            $x4 = $x3  + ($anchor_len-1) * $scale;
            $im -> line($x3,$y3,$x4,$y4,$black);
            $x3 = $max + 1 + $data[6] * $scale;
            $x4 = $max + 1 + ($data[7]-1) * $scale;
            $im -> line($x3,$y3+1,$x4,$y4+1,$red);
           }
        elsif ( !$strand && ($cigar =~ /(\d+)S$/) )
           {
            $x3 = $min - $anchor_len * $scale;
            $x4 = $x3 + ($anchor_len-1) * $scale;
            $im -> line($x3,$y3,$x4,$y4,$black);
            $x3 = $min - $data[7] * $scale;
            $x4 = $min - ($data[6]+1) * $scale;
            $im -> line($x3,$y3+1,$x4,$y4+1,$red);
           }

        ### PE ME end ###
        ### draw read 1 ###
        @data1 = split("", $anchor{$n1}[1]);
        @data2 = split("", $anchor{$n1}[2]);
        @data3 = split("", $anchor{$n1}[3]);
        $gap = 0;
        for ($i=0; $i<=$#data1; $i++)
          {
           if (($data1[$i] eq 'A') || ($data1[$i] eq 'a'))
             {
              $y1 = $row_anchor1 + $edge + $scale;
             }
           elsif ($data1[$i] eq 'C' || ($data1[$i] eq 'c'))
             {
              $y1 = $row_anchor1 + $edge + $scale * 2;
             }
           elsif ($data1[$i] eq 'T' || ($data1[$i] eq 't'))
             {
              $y1 = $row_anchor1 + $edge + $scale * 3;
             }
           elsif ($data1[$i] eq 'G' || ($data1[$i] eq 'g'))
             {
              $y1 = $row_anchor1 + $edge + $scale * 4;
             }
           else
             { ### insertion on the reference read ###
              next;
             }

           $x1 = $edge*20 + ($anchor{$n1}[4] + $i - $gap) * $scale;
           $x2 = $x1 + $scale;
           $y2 = $y1 + $scale;
           if ($data3[$i] ne '-')
             {
              if ($direction eq $strand)
                 { ### anchor split support the same strand MEI ###
                  $im->setPixel($x1,$y1,$red);
                 }
             else
                {
                 $im->setPixel($x1,$y1,$magenta);
                }
             }
           else
             { ### insertion on sequencing read ###
              $gap ++;
              $im->setPixel($x1,$row_anchor1,$blue);
             }
          }
        }
    }

sub draw_MEend
   {
    my $n1 = shift;
    my $l1 = shift;
    my $s1 = shift;
    my $s2 = shift;
    my $s3 = shift;
    my $c1 = shift;
    my $row_read1 = shift;
  
    my (@data1, @data2, @data3);
    my ($i, $gap);
    my ($x1, $y1, $x2, $y2);
# print "MEend of $n1\n";

    @data1 = split("", $s1);
    @data2 = split("", $s2);
    @data3 = split("", $s3);
    $gap = 0;
    for ($i=0; $i<=$#data1; $i++)
      {
        if (($data1[$i] eq 'A') || ($data1[$i] eq 'a'))
          {
           $y1 = $row_read1 + $scale;
          }
        elsif ($data1[$i] eq 'C' || ($data1[$i] eq 'c'))
          {
           $y1 = $row_read1 + $scale * 2;
          }
        elsif ($data1[$i] eq 'T' || ($data1[$i] eq 't'))
          {
           $y1 = $row_read1 + $scale * 3;
          }
        elsif ($data1[$i] eq 'G' || ($data1[$i] eq 'g'))
          {
           $y1 = $row_read1 + $scale * 4;
          }
        else
          { ### insertion on the reference read ###
           next;
          }

        $x1 = $edge*20 + ($c1+$i-$gap) * $scale;
        $x2 = $x1 + $scale;
        $y2 = $y1 + $scale;
        if ($data3[$i] ne '-')
           {
            $im->setPixel($x1,$y1,$red);
           }
        else
           { ### insertion on sequencing read ###
            $gap ++;
            $im->setPixel($x1,$row_read1,$blue);
           }
       }

    ### unmapped sequence in the ME end ###
    @data = split(" ", $l1);
    # print "$n1\t$data[6]\t$data[7]\t$data[10]\n";
    if ($data[8] < $data[9]) 
       {
        if (!exists($split_read{$n1}) && ($data[6] > 0))
           {### unmapped seq ###
            $x1 = $edge*20 + ($c1-$data[6]-1) * $scale; 
            $x2 = $edge*20 + ($c1-1) * $scale; 
            $y1 = $row_read1;
            $y2 = $y1;
            #print "$x1\t$x2\t$c1\t$data[6]\t$gap\n" if ($n1 =~ /2988/);
            if (!exists($transduct{$n1}))
              {
               $im -> line($x1, $y1, $x2, $y2, $black);
              }
            else
              {
               $im -> line($x1, $y1, $x2, $y2, $green);
              }
           }
        if (!exists($split_read{$n1}) && ($data[7] < $data[10]))
           {
            $x1 = $edge*20 + ($c1 + $i - $gap) * $scale; 
            $x2 = $x1 + ($data[10] - $data[7]) * $scale; 
            $y1 = $row_read1;
            $y2 = $y1;
            if (!exists($transduct{$n1}))
               {
                $im -> line($x1, $y1, $x2, $y2, $black);
               }
            else
               {
                $im -> line($x1, $y1, $x2, $y2, $green);
               }
           }
       }
    else
       {
        if (!exists($split_read{$n1}) && ($data[6] > 0))
           {### unmapped seq ###
            $x1 = $edge*20 + ($c1 + $i - $gap) * $scale;
            $x2 = $x1 + $data[6] * $scale; 
            $y1 = $row_read1;
            $y2 = $y1;
            #print "$x1\t$x2\t$c1\t$data[6]\t$gap\n" if ($n1 =~ /2988/);
            if (!exists($transduct{$n1}))
              {
               $im -> line($x1, $y1, $x2, $y2, $black);
              }
            else
              {
               $im -> line($x1, $y1, $x2, $y2, $green);
              }
           }
        if (!exists($split_read{$n1}) && ($data[7] < $data[10]))
           {
            $x2 = $edge*20 + ($c1 - 1) * $scale; 
            $x1 = $x2 - ($data[10] - $data[7]) * $scale; 
            $y1 = $row_read1;
            $y2 = $y1;
            if (!exists($transduct{$n1}))
               {
                $im -> line($x1, $y1, $x2, $y2, $black);
               }
            else
               {
                $im -> line($x1, $y1, $x2, $y2, $green);
               }
           }
       }
   }
