#! /usr/bin/perl -w
#
use strict;

my $seq_file = './L1HS.fa';
my $line;
my $start  = $ARGV[0];
my $end    = $ARGV[1];
my $length = $end- $start + 1;
(open (SEQ, "<$seq_file")) || die "cannot open the seq file\n";

<SEQ>;
my $seq = '';
while ($line = <SEQ>)
    {
     chomp($line);
     $seq .= $line;
    }
$line = substr($seq, $start-1, $length);

print "$line\n";

