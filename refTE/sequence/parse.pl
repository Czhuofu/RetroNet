#! /usr/bin/perl -w
#
use strict;

#my $in_file = './LINE1.fa';
 my $in_file = './ALU.fa';
my $ou_file;
my ($line, @data);
my ($name, $seq);

(open (IN, "<$in_file")) || die "cannot open the in file\n";

while ($line = <IN>)
   {
    chomp($line);
    if ($line =~ /^\>/)
       {
        @data = split("\t", $line);
        $data[0] =~ s/\>//;
        $name = $data[0];
        $seq = '';
        while ( ($line = <IN>) && ($line =~ /\w/) )
            {
             chomp ($line);
             $seq .= $line;
            }
        $ou_file = './'.$name.'.fa';
        (open (OU, ">$ou_file")) || die "cannot open the ou file\n";
        print OU "\>$name\n$seq\n";
        close (OU);
       }
    }

