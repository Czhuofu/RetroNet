#!/usr/bin/env perl 

use strict;
use warnings;

my $rf_file = '';
my $lr_file = '';
my $nb_file = '';
my $ou_file = '';
my $path = $ARGV[0];
my $ver  = $ARGV[1]; 
my $strand = $ARGV[2];
my $sub  = $ARGV[3];
my ($line, @data);
my (%nb, %logr, $rf);

    $rf_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/ALU/'.$sub.'.sr2.pred.summary';
    $lr_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/ALU/'.$sub.'.sr.pred.logR.txt';
    $nb_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/ALU/'.$sub.'.sr.pred.NB.txt';
    $ou_file = $path.'/'.$sub.'/retro_v'.$ver.'_'.$strand.'/ALU/'.$sub.'.sr.prob3.txt';
    %logr = %nb = ();  
   (open (RF, "<$rf_file")) || die "cannot open the RF file $rf_file\n";
   (open (LR, "<$lr_file")) || die "cannot open the RF file $lr_file\n";
   (open (NB, "<$nb_file")) || die "cannot open the RF file $nb_file\n";
   (open (OU, ">$ou_file")) || die "cannot open the RF file\n";
   print OU "chr\tcord1\tcord2\tread\tpos\tneg\tRF\tLR\tNB\n";
   while ($line = <LR>)
      {
       chomp($line);
       @data = split("\t", $line);
       $logr{$data[3]} = $data[4];
     }

   while ($line = <NB>)
      {
       chomp($line);
       @data = split("\t", $line);
       $nb{$data[3]} = $data[4];
     }
    
   while ($line = <RF>)
     {
      chomp($line);
      @data = split("\t", $line);
      $rf =  $data[7];
      print OU "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[5]\t$data[6]\t$rf\t$logr{$data[3]}\t$nb{$data[3]}\n";
     }
  close(RF);
  close(LR);
  close(OU);
  close(NB);
