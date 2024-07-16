#!/usr/bin/env perl

use warnings;
use strict;

my $all_PE    = $ARGV[0];
my $all_SR    = $ARGV[1];
my $masterpath= $ARGV[2];
my $site_file = $ARGV[4].'/01_overlapping_insertions.txt';
my $pe_file   = $ARGV[4].'/02_PEsupport.fa';
my $sr_file   = $ARGV[4].'/03_SRsupport.fa';
my $split_anchor = $ARGV[5];

my %tag = (
   'LINE' => 'L1',
   'ALU'  => 'Alu',
   'SVA'  => 'SVA',
  );
my $tetag = $tag{$ARGV[3]}; 
my ($i, $j, $num, $line, $sr, @sr);
my (@data, %num, @temp, @dat2, @nums, %sr_support, %sr_uniq);
my ($string, %splitanchr, $anchor_direct, $anchor_len, @names, $r1);

(open (SITE, "<$site_file")) || die "cannot open the site file\n";
(open (APE , "<$all_PE   ")) || die "cannot open the all PE file $all_PE\n";
(open (ASR , "<$all_SR   ")) || die "cannot open the all SR file $all_SR\n";
(open (PE  , ">$pe_file  ")) || die "cannot open the out file\n";
(open (SR  , ">$sr_file  ")) || die "cannot open the sr file\n";
(open (ANCHOR, "<$split_anchor")) || die "cannot open the split anchor file\n";

while ($line = <ANCHOR>)
   {
    chomp($line);
    @data = split(" ", $line);
    ($r1, @names) = split(':', $data[0]);
    $line = join(':', @names);
    $splitanchr{$line} = $data[1]."\t".$data[2];
   }

%num = ();
while ($line = <SITE>)
   {
    chomp($line);
    @data = split("\t", $line);
    if ( ($data[3] =~ /pe\,/) && ($data[3] =~ /sr\,/) )
       {
        @temp = split(';', $data[3]);
        for ($j=0; $j<=$#temp; $j++)
            {
             if ($temp[$j] =~ /pe/)
                {
                 @nums = split(',', $temp[$j]);
                 for ($i=4; $i<=$#nums; $i++)
                    {
                     if ($nums[$i])
                       {
                        @dat2 = split('_', $nums[$i]);           
                        $num{$dat2[0]} = 1;
                       }
                    }
                }
             elsif ($temp[$j] =~ /sr/)
                {
                 @nums = split(',', $temp[$j]);
                 @sr   = split(/\|/, $nums[2]);
                 foreach $sr(@sr)
                    {
                     $sr_support{$sr} = 1;
                    }
                }
            }
       }
    elsif ($data[3] =~ /pe\,/)
       {
        @nums = split(',', $data[3]);
        for ($i=4; $i<=$#nums; $i++)
           {
            if ($nums[$i] && ($nums[$i] =~ /\_/))
              {
               @dat2 = split('_', $nums[$i]);        
               $num{$dat2[0]} = 1;
              }
          }
       }
    else
       {
        @sr = split(/[\,\|\;]/, $data[3]);
        foreach $sr(@sr)
          {
           next if (($sr eq 'sr') || ($sr =~ /^\d+$/));
           $sr_support{$sr} = 1 if ($sr =~ /\w/);
          }
       }
   }

$i = 0;
while ($line = <APE>)
   {
    chomp($line);
    $i ++;
    if (exists($num{$i}))
       {
        @data = split("\t", $line);
        print PE ">$data[4]\;$data[0]\;$data[1]\;$data[2]\;$data[5]\;$data[6]\;$data[15]\;$data[16]\n$data[12]\n";
        if (exists($splitanchr{$data[4]}))
           {
            ($anchor_direct, $anchor_len) = split("\t", $splitanchr{$data[4]});
            if ( ($data[16] =~ /^(\d+)S/) && ($1 > 19) )
               {
                $string = substr($data[20], 0, $anchor_len);
                print PE ">$data[4]\_anchor\;$data[0]\;$data[1]\;$data[2]\;$anchor_direct\;$anchor_len\;$data[16]\n$string\n";
               }
            if ( ($data[16] =~ /(\d+)S$/) && ($1 > 19) )
               {
                $string = substr($data[20], length($data[20]) - $anchor_len, $anchor_len);
                print PE ">$data[4]\_anchor\;$data[0]\;$data[1]\;$data[2]\;$anchor_direct\;$anchor_len\;$data[16]\n$string\n";
               }
           }    
       }
   }

while ($line = <ASR>)
   {
    chomp($line);
    @data = split("\t", $line);
    if ( ($data[3] =~ /$tetag/) && (exists($sr_support{$data[4]})) )
       {
        $sr = "$data[4]\;$data[0]\;$data[1]\;$data[2]\;$data[5]\;$data[6]";
        $sr_uniq{$sr} = $data[17];
       }
   }

for $sr(sort keys %sr_uniq)
   {
    print SR ">$sr\n$sr_uniq{$sr}\n";
   }

