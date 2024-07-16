#!/usr/bin/perl -w

use strict;

my $seq_file = '/home/txw/hg19G37/human_g1k_v37.fasta';

my ($num, $line, $chr, $seq, %refseq);
my $start = 1; 

(open (SEQ, "<$seq_file")) || die "cannot open the seq file\n";

    while ($line = <SEQ>)
      {
       if ( $line =~ /^\>([^ ]+) / ) # || ($line =~ /^\>([XYMT]+) /) )
          {
           if (!$start)
             {
              $num = length($seq);
              $refseq{$chr} = $seq;
              print "$chr\t$num\n";
             }
           else
             {
              $start = 0;
             }
           $seq = '';
           $chr = $1;

           # print "read seq chr $chr\n";
          }
       else
          {
           chomp($line);
           $line =~ tr/[a-z]/[A-Z]/;
           $seq .= $line;
          }
      }
    $refseq{$chr} = $seq;
    $num = length($seq);
    print "$chr\t$num\n";
