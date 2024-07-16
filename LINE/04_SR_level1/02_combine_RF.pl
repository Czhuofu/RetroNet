#! /usr/bin/perl -w

use strict;

my $workfolder = $ARGV[0];
my $retro = 'retro_v'.$ARGV[1].'_'.$ARGV[2];
my $sample = $ARGV[3];
my $m1_file  = $workfolder.'/'.$sample.'/'.$retro.'/'.$sample.'.sr.LINE.matrix';
my $g1_file  = $workfolder.'/'.$sample.'/'.$retro.'/LINE/'.$sample.'.sr.pred.G1.txt';
my $g2_file  = $workfolder.'/'.$sample.'/'.$retro.'/LINE/'.$sample.'.sr.pred.G2.txt';
my $ou1_file = $workfolder.'/'.$sample.'/'.$retro.'/LINE/'.$sample.'.sr2.pred.all.txt';
my $ou2_file = $workfolder.'/'.$sample.'/'.$retro.'/LINE/'.$sample.'.sr2.pred.summary';
my @order = (1,2); 
my ($line, @data);
my (%pred1, %pred2);
my ($i, $pred, $id, $pos, $neg, $g1, $g2);

(open (M1, "<$m1_file")) || die "cannot open the mat file\n";
(open (G1, "<$g1_file")) || die "cannot open the g1 file\n";
while ($line = <G1>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred1{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (G2, "<$g2_file")) || die "cannot open the g2 file\n";
while ($line = <G2>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred2{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (OU1, ">$ou1_file")) || die "cannot open the out file\n";
(open (OU2, ">$ou2_file")) || die "cannot open the out file\n";

print OU1 "chr\tcord1\tcord2\tread\tpos\tneg\tG1\tG2\tG3\n";
<M1>;
while ($line = <M1>)
   {
    chomp($line);
    @data = split("\t", $line);
    if ( (($data[6] == 1) || ($data[7] == 1)) && ($data[13] == 0) )
       {
        $id = $data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3];
        $pos = $data[6];
        $neg = $data[7];
        $g1 = (exists($pred1{$id})) ? $pred1{$id} : 'NA';
        $g2 = (exists($pred2{$id})) ? $pred2{$id} : 'NA';
        print OU1 "$id\t$pos\t$neg\t$g1\t$g2\n";
       }
   }

close(OU1);
(open (OU1, "<$ou1_file")) || die "cannot open the out file\n";

my ($pos_0, $pos_1, $pos_na, $pos_total, $neg_0, $neg_1, $neg_na, $neg_total, $recall, $precision);
my ($basep_0, $basep_1, $basep_na, $basen_0, $basen_1, $basen_na, $base_recall, $base_precision);
my ($rnd_recall, $rnd_precision);
my %usage = (
   '1' => 0,
   '2' => 0,
  ); 

<OU1>;
while ($line = <OU1>)
   {
    chomp($line);
    @data = split("\t", $line);
    $pred = 'NA';
    for ($i=0; $i<=$#order; $i++)
       {
        if ($data[5+$order[$i]] ne 'NA')
           {
            $usage{$order[$i]} ++;
            $pred = ($data[5+$order[$i]] > 0.5) ? 1 : 0;
            last; 
           }
       }
    print OU2 "$line\t$pred\n";
   }

