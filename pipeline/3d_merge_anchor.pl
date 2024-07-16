#!/usr/bin/env perl

use warnings;
use strict;
### PE reads for SVA1 insertions ###

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
my $TE_family = $ARGV[6];

### look for PE supporting reads with split-read anchors ###
my $site_file  = $outpath.'/'.$sub.'/'.$retro.'/'.$TE_family.'/'.$sub.'.'.$TE_family.'.novel.sites';
my $fasta_temp = $outpath.'/'.$sub.'/'.$retro.'/'.$TE_family.'/temp.fa';
my $exo_temp   = $outpath.'/'.$sub.'/'.$retro.'/'.$TE_family.'/temp.exo.txt';
my $seg_temp   = $outpath.'/'.$sub.'/'.$retro.'/'.$TE_family.'/temp.seg.txt';
my $sanch_file = $outpath.'/'.$sub.'/'.$retro.'/'.$TE_family.'/split-reads-anchor.txt';
my $TE_familyname = $TE_family;
if ($TE_family eq 'LINE')
   {
   $TE_familyname = 'L1HS';
   }

my $ref_file = $masterpath.'/refTE/sequence/'.$TE_familyname.'.fa';

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

