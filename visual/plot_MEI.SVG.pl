#!/usr/bin/env perl

use warnings;
use strict;

use GD;
use GD::Arrow;
use GD::SVG;

### adding unit bars ###
### fixing repetitive representation ###
### adding read ID ###
### add SR read plot ###
### flip reference TE for -strand insertions ###

my $masterpath= $ARGV[0];
my $chr     = $ARGV[2];
my $te      = $ARGV[3];
my $hg      = $ARGV[5];
my $strand  = $ARGV[6];
my $in_file = $ARGV[4].'/temp'.$te.'/04_PE.alignment';
my $sr_file = $ARGV[4].'/temp'.$te.'/05_SR.alignment';
my $ou_file = $ARGV[4].'/'.$ARGV[1];

(open (OU, ">$ou_file")) || die "cannot open the ou file$ou_file\n";

my $min_hg = 1000000000;
my $max_hg = 0;
my $min_te = 1000000000;
my $max_te = 0;
my $WIN    = 100;
my $max    = 0;
my $unit_t = 5;
my $unit_s = 200;
my ($arrow, $x0, $y0, $x1, $y1, $x2, $y2, $x3, $y3, $x4, $y4, $im);
my ($read_i, $line, $i, $read, $cord, $string);
my ($num, $hg_start, $te_start, $te_end, $x5, $x6);
my (@lib, @data, @name, @reads);
my $left = 20;
my $up = 20;
my $width = 50;
my ($c1, $c2, $c3, $direction);
my ($cigar, $insert, $rev, $name,  $j, @sr_reads, $map);

$i = 0;
if ( (-e $in_file) && (! (-z $in_file)) )
   {
    (open (IN, "<$in_file")) || die "cannot open the in file\n";
    
    while ($line = <IN>)
       {
        chomp($line);
        @data = split(" ", $line);
        next if ($data[5] ne $te);
        #print "$data[1]\n";
        @name = split(';', $data[1]);
        $read = $name[0];
        #print "$name[2]\n";
        $cord = $name[1].'_'.$name[2].'_'.$name[3];
        $min_hg = ($min_hg < $name[2]) ? $min_hg : $name[2];
        $max_hg = ($max_hg > $name[3]) ? $max_hg : $name[3];
        $min_te = ($min_te < $data[8]) ? $min_te : $data[8];
        $min_te = ($min_te < $data[9]) ? $min_te : $data[9];
        $max_te = ($max_te > $data[8]) ? $max_te : $data[8];
        $max_te = ($max_te > $data[9]) ? $max_te : $data[9];
        $reads[$i][0] = $name[2];
        $reads[$i][1] = $line;
        $i ++ ;
       }
   }

if ( (-e $sr_file) && (! (-z $sr_file)) )
   {
    (open (SR, "<$sr_file")) || die "cannot open the sr file $sr_file\n";
    $j =0;
    while ($line = <SR>)
       {
        chomp($line);
        @data = split(" ", $line);
        @name = split(';', $data[1]);
        $min_hg = ($min_hg < $name[2]-150) ? $min_hg : ($name[2]-150);
        $max_hg = ($max_hg > $name[3]+150) ? $max_hg : ($name[3]+150);
        $min_te = ($min_te < $data[8]) ? $min_te : $data[8];
        $min_te = ($min_te < $data[9]) ? $min_te : $data[9];
        $max_te = ($max_te > $data[8]) ? $max_te : $data[8];
        $max_te = ($max_te > $data[9]) ? $max_te : $data[9];
        $sr_reads[$j][0] = $name[0];
        $sr_reads[$j++][1] = $line;
        $i ++;
       }
   }
   
$read_i = $i; ### total support reads ###
#print "hg\t$min_hg\t$max_hg\n";
#print "te\t$min_te\t$max_te\n";
$min_hg = $min_hg - $WIN if ($min_hg > $WIN);
$min_te = $min_te - $WIN if ($min_te > $WIN);
$max = ($max_hg-$min_hg > $max_te-$min_te) ? ($max_hg-$min_hg) : ($max_te-$min_te);
### start the image ###
$x0 = 500 + $max;
$y0 = 500 + $width * $i * 2.1;
print "$x0\t$y0\n";
$im = new GD::SVG::Image($x0,$y0);

### allocate some colors
my $white = $im->colorAllocate(255,255,255);
my $black = $im->colorAllocate(0,0,0);       
my $red   = $im->colorAllocate(190,30,45);      
my $blue  = $im->colorAllocate(0,173,238);
my $magenta= $im->colorAllocate(255,0,255);

### plot the human reference line ###
$x1 = ($x0 - $max_hg + $min_hg) / 2;
$y1 = 200;
$x2 = $x1 + $max_hg - $min_hg;
$y2 = 200;
$im->setThickness(3);
$im->line($x1, $y1, $x2, $y2, $black);
$hg_start = $x1;

### Add 200bp marks ###
$num = int(($max_hg - $min_hg) / $unit_s);
for ($i=1; $i<=$num; $i++)
   {
    $x1 += $unit_s;
    $x2 = $x1;
    $y1 = 200;
    $y2 = $y1 + $unit_t;
    $im->line($x1, $y1, $x2, $y2, $black);
   }
$im->setThickness(1);

### label reference genome and chromosome ID ###
$string = $hg.'_'.$chr;
$x1 = ($x0 - $max_hg + $min_hg) / 2 - 120;
$y1 = $y1 - 5;
$im->string(gdGiantFont,$x1,$y1,$string,$black);

### label human reference genome coordinates ###
$string = $min_hg;
$x1 = ($x0 - $max_hg + $min_hg) / 2 - $left;
$y1 = 200 - $up;
$im->string(gdGiantFont,$x1,$y1,$string,$black);

$string = $max_hg;
$x1 =  $x0/2 + ($max_hg - $min_hg) / 2 - $left;
$y1 = 200 - $up;
$im->string(gdGiantFont,$x1,$y1,$string,$black);

### plot the TE line ###
$x1 = ($x0 - $max_te + $min_te) / 2;
$y1 = $y0 - 200;
$x2 = $x1 + $max_te - $min_te;
$y2 = $y1;
$im->setThickness(3);
$im->line($x1, $y1, $x2, $y2, $black);
$te_start = $x1;
$te_end = $x2;

### mark very 200 bp ###
$num = int(($max_te - $min_te) / $unit_s);
for ($i=1; $i<=$num; $i++)
   {
    $x1 += $unit_s;
    $x2 = $x1;
    $y1 = $y0 - 200;
    $y2 = $y1 + $unit_t;
    $im->line($x1, $y1, $x2, $y2, $black);
   }
$im->setThickness(1);

### label the type of TE ###
$string = ($strand) ? $te.'+' : $te.'-';
$x1 = ($x0 - $max_te + $min_te) / 2 - 80;
$y1 = $y0 - 200 - 10;
$im->string(gdGiantFont,$x1,$y1,$string,$black);

### label the TE coordinates ###
### flip the line if insertion is on - strand ###
$string = ($strand) ? $min_te : $max_te;
$x1 = ($x0 - $max_te + $min_te) / 2 - $left;
$y1 = $y0 - 200 - $up;
$im->string(gdGiantFont,$x1,$y1,$string,$black);

$string = ($strand) ? $max_te : $min_te;
$x1 = $x0/2 + ($max_te - $min_te) / 2 - $left;
$y1 = $y0 - 200 - $up;
$im->string(gdGiantFont,$x1,$y1,$string,$black);

### plot the PE reads ###
@reads = sort {$a->[0] <=> $b->[0]} @reads;
for ($i = 0; $i <=$#reads;$i++)
   {
    @data = split(" ", $reads[$i][1]);
    @name = split(';', $data[1]);
    $c1 = ($name[4] eq '+') ? 1 : -1;
    $c2 = ($name[5] eq '+') ? 1 : -1;
    $c3 = ($name[6] eq '+') ? 1 : -1;
    $direction = -1 * $c1 * $c2 * $c3; ### equal to $strand ###
    @lib  = split('_', $name[0]);
    $string = $i + 1;
    $string = $lib[0];;
    if ($name[4] eq '+')
       {
        ### upstream read ###
        $x1 = $hg_start + $name[3] - $min_hg;
        $y1 = 250 + $i * 2 * $width + $width/2;
        $x2 = $hg_start + $name[2] - $min_hg;
        $y2 = $y1;

        if ($name[6] eq '+')
           {
            $x3 = $te_start + $data[9] - $min_te;
            $x4 = $te_start + $data[8] - $min_te;
           }
        else
           {
            $x4 = $te_start + $data[9] - $min_te;
            $x3 = $te_start + $data[8] - $min_te;
           }     
        $y3 = $y1 + $width; #250 + $read_i*$width + $i * $width + $width/2;
        $y4 = $y3;
       }
    else
       {
        ### downstream read ###
        $x1 = $hg_start + $name[2] - $min_hg;
        $y1 = 250 + $i * 2 * $width + $width/2;
        $x2 = $hg_start + $name[3] - $min_hg;
        $y2 = $y1;

        ### hg is the second end ###
        if ($name[6] eq '+')
           { 
            $x3 = $te_start + $data[9] - $min_te;
            $x4 = $te_start + $data[8] - $min_te;
           }
        else
           {
            $x4 = $te_start + $data[9] - $min_te;
            $x3 = $te_start + $data[8] - $min_te;
           }
        $y3 = $y1 + $width; #250 + $read_i*$width + $i * $width + $width/2;
        $y4 = $y3;
       }
    $arrow = GD::Arrow::Full->new( 
                  -X1    => $x1, 
                  -Y1    => $y1, 
                  -X2    => $x2, 
                  -Y2    => $y2, 
                  -WIDTH => $width/4,
              );
    $im->filledPolygon($arrow,$blue);

    if ($direction > 0)
       { ### insertion on + strand ###
        $arrow = GD::Arrow::Full->new(
                  -X1    => $x3,
                  -Y1    => $y3,
                  -X2    => $x4,
                  -Y2    => $y4,
                  -WIDTH => $width/4,
              );
        $im->filledPolygon($arrow,$red);
        $im->line($x1,$y1,$x3,$y3,$black);
       }
    else
       { ### insertion on - strand ###
        $x5 = $te_end - ($x3 - $te_start);
        $x6 = $te_end - ($x4 - $te_start);
        $arrow = GD::Arrow::Full->new(
                  -X1    => $x5,
                  -Y1    => $y3,
                  -X2    => $x6,
                  -Y2    => $y4,
                  -WIDTH => $width/4,
              );
        $im->filledPolygon($arrow,$magenta);
        $im->line($x1,$y1,$x5,$y3,$black);
       }
   if ($name[4] eq '+')
       {
        $x1 += 5;
       }
    else
       {
        $x1 -=30;
       }
    $y1 -= 15;
    $im->string(gdGiantFont,$x1,$y1,$string,$black);
   }

### plot the SR reads ###
for ($j=0; $j<=$#sr_reads; $j++)
   {
    @data = split(" ", $sr_reads[$j][1]);
    @name = split(';', $data[1]);
    $cigar = $name[5];
    $rev   = ($name[4] & 0x10) ? 1 : 0;
    if ($cigar =~ /(\d+)M(\d+)S/)
      { ### human => L1 ###
       $map = $1;
       $insert = $2;
       if ($rev == 0)
          { ### FWD strand ###
           $x1 = $hg_start + $name[2] - $min_hg;
           $y1 = 250 + $i * 2 * $width + $width/2;
           $x2 = $hg_start + $name[2] - $min_hg - $map;
           $y2 = $y1;
          }
       else
          {
           $x1 = $hg_start + $name[2] - $min_hg - $map;
           $y1 = 250 + $i * 2 * $width + $width/2;
           $x2 = $hg_start + $name[2] - $min_hg;
           $y2 = $y1;
          }
        $arrow = GD::Arrow::Full->new(
               -X1    => $x1,
               -Y1    => $y1,
               -X2    => $x2,
               -Y2    => $y2,
               -WIDTH => $width/4,
           );
        $im->filledPolygon($arrow,$blue);

        $x3 = $hg_start + $name[2] - $min_hg;
        $y3 = $y1 - $width/4;
        $x4 = $x3 + $insert;
        $y4 = $y3 + $width/2;
        $im->rectangle($x3,$y3,$x4,$y4,$black);

        $x3 = $te_start - $min_te + $data[9];
        $y3 = $y1 + $width;
        $x4 = $te_start - $min_te + $data[8];
        $y4 = $y3;
        if ($data[8] < $data[9])
          {
           $arrow = GD::Arrow::Full->new(
                  -X1    => $x3,
                  -Y1    => $y3,
                  -X2    => $x4,
                  -Y2    => $y4,
                  -WIDTH => $width/4,
              );
           $im->filledPolygon($arrow,$red);
          }
        else
          {
           $x5 = $te_end - ($x3 - $te_start);
           $x6 = $te_end - ($x4 - $te_start);
           $arrow = GD::Arrow::Full->new(
                  -X1    => $x5,
                  -Y1    => $y3,
                  -X2    => $x6,
                  -Y2    => $y4,
                  -WIDTH => $width/4,
              );
           $im->filledPolygon($arrow,$magenta);
          }
       }
    elsif ($cigar =~ /(\d+)S(\d+)M/)
       { ### L1seq + human seq ###
        $insert = $1;
        $map = $2;
        if ($rev == 0)
          { ### FWD strand ###
           $x1 = $hg_start + $name[2] - $min_hg + $map;
           $y1 = 250 + $i * 2 * $width + $width/2;
           $x2 = $hg_start + $name[2] - $min_hg;
           $y2 = $y1;
          }
        else
          {
           # print "reverse strand\n";
           $x1 = $hg_start + $name[2] - $min_hg;
           $y1 = 250 + $i * 2 * $width + $width/2;
           $x2 = $hg_start + $name[2] - $min_hg + $map;
           $y2 = $y1;
          }
        $arrow = GD::Arrow::Full->new(
               -X1    => $x1,
               -Y1    => $y1,
               -X2    => $x2,
               -Y2    => $y2,
               -WIDTH => $width/4,
           );
        $im->filledPolygon($arrow,$blue);

        $x3 = $hg_start + $name[2] - $min_hg - $insert;
        $y3 = $y1 - $width/4;
        $x4 = $x3 + $insert;
        $y4 = $y3 + $width/2;
        $im->rectangle($x3,$y3,$x4,$y4,$black);

        $x3 = $te_start - $min_te + $data[9];
        $y3 = $y1 + $width; 
        $x4 = $te_start - $min_te + $data[8];
        $y4 = $y3;
        if ($data[8] < $data[9])
          {
           $arrow = GD::Arrow::Full->new(
                  -X1    => $x3,
                  -Y1    => $y3,
                  -X2    => $x4,
                  -Y2    => $y4,
                  -WIDTH => $width/4,
              );
           $im->filledPolygon($arrow,$red);
          }
        else
          {
           $x5 = $te_end - ($x3 - $te_start);
           $x6 = $te_end - ($x4 - $te_start);
           $arrow = GD::Arrow::Full->new(
                  -X1    => $x5,
                  -Y1    => $y3,
                  -X2    => $x6,
                  -Y2    => $y4,
                  -WIDTH => $width/4,
              );
           $im->filledPolygon($arrow,$magenta);
          }
       }
    $x1 += 5;
    $y1 -= 7;
    @name = split(/[\;\:\_]/, $sr_reads[$j][0]);
    #$name = $cigar.':'.$name[0].':'.$name[4].':'.$name[5].':'.$name[6].':'.$name[7];   
    $name = $sr_reads[$j][0];
    $im->string(gdGiantFont,$x1,$y1,$name,$black);
    $i ++;
   }
  
# make sure we are writing to a binary stream
binmode OU;

# Convert the image to PNG and print it on standard output
print OU $im->svg;
close OU;

