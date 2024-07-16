#! /usr/bin/perl -w

use strict;

### preparing the windows for checking sequencing depth around each supporting read ###
my $masterpath = $ARGV[0]; ### masterpath ###
my $sub  = $ARGV[1];   ### subjec id ###
my $WIN  = $ARGV[2];   ### add $WIN to the left and the right of the supporting reads, then check the sequencing depth in the range ###
my $fold = $ARGV[3];   ### $outpath/$sub ###
my $retro= $ARGV[4];   ### RetroSom version control ###
my $ref  = $ARGV[5];   ### hg38 or hg19 ###
my $MAP  = 85; ### cutoff for the identity between supporting reads and ME consensus ###
my $in1_file  = $fold.'/'.$retro.'/'.$sub.'.discover'; ### PE supporting reads ###
my $in2_file  = $fold.'/'.$retro.'/'.$sub.'.sr.discover'; ### SR supporting reads ###
my $out_file  = $fold.'/'.$retro.'/'.$sub.'.bed'; ### output BED file, each row represents one supporting read ###
my $size_file = $masterpath.'/refTE/size/'.$ref.'.size'; ### length of the chromosomes ###
my ($i, $read, $line, $left, $right, %size, @bed, @data);
(open (SIZE, "<$size_file")) || die "cannot open the size file\n";
(open (IN1 , "<$in1_file ")) || die "cannot open the in1  file\n";
(open (IN2 , "<$in2_file ")) || die "cannot open the in2  file\n";
(open (OUT , ">$out_file ")) || die "cannot open the out file\n";

while ($line = <SIZE>)
   {
    chomp($line);
    @data = split("\t", $line);
    $size{$data[0]} = $data[1];
   }

$i = 0;
while ($line = <IN1>)
   { ### PE reads ###
    chomp($line);
    @data = split("\t", $line);
    $left = $data[1] - $WIN; ### the boudaries are +/- 300 bp of the PE supporting read ###
    $right= $data[2] + $WIN;
    next if ($data[7] < $MAP); ### delete read if the identity is lower than 85% ###
    next if (!exists($size{$data[0]})); ### delete read if it's on a chromosome with no documented length ###
    $left = 0 if ($left < 0); ### the left most boundary is 0 ###
    $right = $size{$data[0]} if ($right > $size{$data[0]}); ### the right most boundary is the length of chromosome ###
    $read = ($data[14] & 0x40) ? "r1" : "r2"; ### read1 or read2 in paired-end sequencing ###
    $bed[$i][0] = $data[0]; ### chromosome name ###
    $bed[$i][1] = $left; ### left boundary ###
    $bed[$i][2] = $right; ### right boundary ###
    $bed[$i][3] = 'PE__'.$read.':'.$data[4]; ###  unique ID for this supporting read ###
    $i ++;
   }

while ($line = <IN2>)
   {### SR reads ###
    chomp($line);
    @data = split("\t", $line);
    $left = $data[1] - $WIN; ### the boudaries are +/- 300 bp of the SR supporting read ###
    $right= $data[2] + $WIN;
    next if ($data[8] < $MAP); ### delete read if the identity is lower than 85% ###
    next if (!exists($size{$data[0]})); ### delete read if it's on a chromosome with no documented length ###
    $left = 0 if ($left < 0); ### the left most boundary is 0 ###
    $right = $size{$data[0]} if ($right > $size{$data[0]}); ### the right most boundary is the length of chromosome ###
    $read = ($data[5] & 0x40) ? "r1" : "r2";  ### read1 or read2 in paired-end sequencing ###
    $bed[$i][0] = $data[0]; ### chromosome name ###
    $bed[$i][1] = $left; ### left boundary ###
    $bed[$i][2] = $right; ### right boundary ###
    $bed[$i][3] = 'SR__'.$read.':'.$data[4]; ### unique ID for this supporting read ###
    $i ++;
   }

@bed = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]} @bed; ### sort the windows by coordinates ###

for ($i=0; $i<=$#bed; $i++)
   {
    print OUT "$bed[$i][0]\t$bed[$i][1]\t$bed[$i][2]\t$bed[$i][3]\n";
   }


