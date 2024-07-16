#!/usr/bin/env perl

use warnings;
use strict;
 
my $sr_file    = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.sr.tabe.discover';
my $ou_file    = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.sr.dedup.discover';
my $dup_file   = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.sr.dup.txt';
my $clust_file = $ARGV[1].'/'.$ARGV[2].'/'.$ARGV[0].'.sr.cluster.txt';

my ($line, @data);
my ($chr, $te, $read, $anchor_left, $anchor_cigar, $anchor_right, $gap, $tepos1, $tepos2);
my ($i, $j, $k, $read1, $read2, $cigar, $flag, $name, $dupl, $first);
my (%point, %id, %dup, %pair, @all, %assign, %remove);

(open (SR   , "<$sr_file   ")) || die "cannot open the sr file$sr_file\n";
(open (OU   , ">$ou_file   ")) || die "cannot open the ou file\n";
(open (DUP  , ">$dup_file  ")) || die "cannot open the dup file\n";
(open (CLUST, ">$clust_file")) || die "cannot open the clust file\n";

while ($line = <SR>)
   {
    chomp($line);
    @data = split("\t", $line);
    $chr          = $data[0];
    #$te           = $chr.'_'.$data[3];
    $read         = $data[4];
    $flag         = $data[5];
    $cigar        = $data[6];
    $anchor_left  = $data[24];
    $anchor_cigar = $data[19];
    $anchor_right = $anchor_left + _getRightmostCord($anchor_left, $anchor_cigar);
    $gap          = abs($data[9]);
    $tepos1       = $data[13];
    $tepos2       = $data[14];
    $te = $anchor_left.'__'.$data[3];
    $point{$te}   = 0 if (!exists($point{$te}));
    $id{$te}[$point{$te}][0] = $flag."\t".$read;
    $id{$te}[$point{$te}][1] = $chr.'___'.$anchor_left.'___'.$anchor_right.'___'.$gap.'___'.$tepos1.'___'.$tepos2.'___'.$cigar.'___'.$anchor_cigar;
    $id{$te}[$point{$te}][2] = $read;
    $point{$te} ++;
    $te = $anchor_right.'__'.$data[3];
    $point{$te}   = 0 if (!exists($point{$te}));
    $id{$te}[$point{$te}][0] = $flag."\t".$read;
    $id{$te}[$point{$te}][1] = $chr.'___'.$anchor_left.'___'.$anchor_right.'___'.$gap.'___'.$tepos1.'___'.$tepos2.'___'.$cigar.'___'.$anchor_cigar;
    $id{$te}[$point{$te}][2] = $read;
    $point{$te} ++;
   }

for $te(keys %id)
   {
    next if ($point{$te} < 2);
    for ($i=0; $i<($point{$te}-1); $i++)
       {
#next if ( ($id{$te}[$i][0] !~ /$ARGV[3]/) && ($id{$te}[$i][0] !~ /$ARGV[4]/) );
        $read1 = $id{$te}[$i][1];
        for ($j=$i+1; $j<$point{$te}; $j++)
           {
#next if ( ($id{$te}[$j][0] !~ /$ARGV[3]/) && ($id{$te}[$j][0] !~ /$ARGV[4]/) );
            $read2 = $id{$te}[$j][1];
            next if ($id{$te}[$i][0] eq $id{$te}[$j][0]);
            $dupl  = _testduplicate($read1, $read2);   
#print "$dupl\t$id{$te}[$i][0]\t$read1\t$id{$te}[$j][0]\t$read2\n" if ( ($id{$te}[$i][0] =~ /$ARGV[3]/) && ($id{$te}[$j][0] =~ /$ARGV[3]/) );
#print "$dupl\t$id{$te}[$i][0]\t$read1\t$id{$te}[$j][0]\t$read2\n" if ( ($id{$te}[$i][0] =~ /$ARGV[4]/) && ($id{$te}[$j][0] =~ /$ARGV[4]/) );
            if ($dupl)
               {
                $dup{$id{$te}[$i][0]} = 1;
                $dup{$id{$te}[$j][0]} = 1;
                $pair{$id{$te}[$i][0]}{$id{$te}[$j][0]} = 1;
                print DUP "$te\t$id{$te}[$i][0]\t$id{$te}[$i][1]\t$id{$te}[$j][0]\t$id{$te}[$j][1]\n";
               }
           }
       }
   }

$i = scalar keys %dup;
print "duplicate reads $i\n";
$i = 0;
for $name(keys %dup)
   {
    $all[$i++] = $name;
   }

$i = 0;
for ($j=0; $j<$#all; $j++)
   {
    if (!exists($assign{$all[$j]}))
       {
        $assign{$all[$j]} = $i;
        $i ++;
       }
    for ($k=$j+1; $k<=$#all; $k++)
       {
        if (exists($pair{$all[$k]}{$all[$j]}) || exists($pair{$all[$j]}{$all[$k]}) )
           {
            $assign{$all[$k]} = $assign{$all[$j]} if (!exists($assign{$all[$k]}));
           }
       }
   }

print "$i cluster\n";
for ($j=0; $j<$i; $j++)
   {
    print CLUST "$j";
    $first = 1;
    for ($k=0; $k<=$#all; $k++)
       {
        if ($assign{$all[$k]} eq $j)
           {
            print CLUST "\t$all[$k]";
            if ($first == 0)
               {
                $remove{$all[$k]} = 1;
               }
            $first = 0;
           }
       }
    print CLUST "\n";
   }

$j = scalar keys %dup;
$j = $j - $i;
print "remove $j reads\n";

$j = scalar keys %remove;
print "removed $j reads\n";

(open (SR, "<$sr_file")) || die "cannot open the sr file\n";
while ($line = <SR>)
   {### print new sites ###
    chomp($line);
    @data = split("\t", $line);
    if (!exists($remove{$data[5]."\t".$data[4]}))
       {
        $line = join("\t", @data);
        print OU "$line\n"
       }
   }

sub _getRightmostCord
{
my $left  = shift;
my $sam_cigar = shift;
my ($M, $D, $N);
my $right = 0;

$M = $D = $N = 0;
while ($sam_cigar =~ /(\d+)M/g)
   {
    $M += $1;
   }
while ($sam_cigar =~ /(\d+)D/g)
   {
    $D += $1;
   }
while ($sam_cigar =~ /(\d+)N/g)
   {
    $N += $1;
   }
$right = $M + $D + $N;
return ($right);
}

sub _testduplicate
{ ###  $chr.'___'.$anchor_left.'___'.$anchor_right.'___'.$gap.'___'.$tepos1.'___'.$tepos2.'___'.$CIGAR ###
### same chromosome                                                              ###
### one of the two anchor coordiantes                                            ###
### similar CIGAR XMYS or XSYM, && identical X or identical Y                    ###
### one of the two TE coordiantes                                                ###
### if anchor is split:similar CIGAR XMYS or XSYM, && identical X or identical Y ###
####################################################################################
my $sr1 = shift;
my $sr2 = shift;
my @data1 = split('___', $sr1);
my @data2 = split('___', $sr2);
my ($flag1, $flag2, $flag3, $flag4, $flag5, $flag6);
my ($c1, $c2, $c3, $c4);
if ($data1[0] ne $data2[0])
   {
    return(0);
   }
if (($data1[1] == $data2[1]) || ($data1[2] == $data2[2]))
   {
    $flag1 = 1;
   }
else
   {
    $flag1 = 0;
    return(0);
   }
$flag2 = ($data1[3] == $data2[3]) ? 1 : 0;
$flag3 = ($data1[4] == $data2[4]) ? 1 : 0;
$flag4 = ($data1[5] == $data2[5]) ? 1 : 0;
$flag5 = 0;
$flag6 = 1;
return (0) if ($flag3 == 0 && $flag4 == 0);

if ($data1[6] =~ /(\d+)M(\d+)S/)
   {
    $c1 = $1;
    $c2 = $2;
    if ($data2[6] =~ /(\d+)M(\d+)S/)
       {
        $c3 = $1;
        $c4 = $2;
        $flag5 = 1 if ( ($c1 == $c3) || ($c2 == $c4) );
       }
   }

if ($data1[6] =~ /(\d+)S(\d+)M/)
   {
    $c1 = $1;
    $c2 = $2;
    if ($data2[6] =~ /(\d+)S(\d+)M/)
       {
        $c3 = $1;
        $c4 = $2;
        $flag5 = 1 if ( ($c1 == $c3) || ($c2 == $c4) );
       }
   }
return (0) if ($flag5 == 0);
### check anchor read cigar ###
if (($data1[7] !~ /S/) && ($data2[7] !~ /S/))
   {
    $flag6 = 1;
   }
elsif (($data1[7] =~ /S/) && ($data2[7] !~ /S/))
   {
    $flag6 = 0;
   }
elsif (($data1[7] !~ /S/) && ($data2[7] =~ /S/))
   {
    $flag6 = 0;
   }
elsif ($data1[7] =~ /(\d+)M(\d+)S/)
   {
    $c1 = $1;
    $c2 = $2;
    if ($data2[7] =~ /(\d+)M(\d+)S/)
       {
        $c3 = $1;
        $c4 = $2;
        $flag6 = 0 if ( ($c1 != $c3) && ($c2 != $c4) );
       }
    else
       {
        $flag6 = 0;
       }
   }
elsif ($data1[7] =~ /(\d+)S(\d+)M/)
   {
    $c1 = $1;
    $c2 = $2;
    if ($data2[7] =~ /(\d+)S(\d+)M/)
       {
        $c3 = $1;
        $c4 = $2;
        $flag6 = 0 if ( ($c1 != $c3) && ($c2 != $c4) );
       }
    else
       {
        $flag6 = 0; 
       }
   }

if ($flag1 && ($flag3 || $flag4) && $flag5 && $flag6)
   {
    return(1);
   }
else
   {
    return(0);
   }
}

