#!/usr/bin/perl -w

### direction of SR and PE support reads in novel.calls ###

use strict;

my $tag   = $ARGV[0];
my $te    = $ARGV[1];
my $retro = $ARGV[2];
my $path  = $ARGV[3];

my %TEtag = (
   'ALU' => 'Alu',
   'LINE' => 'L1',
  );

my $tetag = $TEtag{$te};
my $site_file = $path.'/'.$retro.'/'.$te.'/'.$tag.'.'.$te.'.novel.sites';
my $sr_file   = $path.'/'.$retro.'/'.$tag.'.sr.dedup.discover';
my $call_file = $path.'/'.$retro.'/'.$te.'/'.$tag.'.'.$te.'.novel.calls';
my $dire_file = $path.'/'.$retro.'/'.$te.'/'.$tag.'.'.$te.'.direction.both.txt';

(open (SITE, "<$site_file")) || die "cannot open the site file $site_file\n";
(open (SR  , "<$sr_file  ")) || die "cannot open the sr file\n";
(open (CALL, "<$call_file")) || die "cannot open the call file\n";
(open (DIRE, ">$dire_file")) || die "cannot open the direction file\n";

my $count = 0;
my ($line, @data, @data1, @data2, @data3, @direction, %sr);
my ($id, $pe, $count_pos, $count_neg, $i, $j, $c1, $c2, $c3, $direction, $read);

while ($line = <SR>)
   {
    chomp($line);
    next if ($line !~ /$tetag/);
    @data = split("\t", $line);
    $direction = ($data[13] < $data[14]) ? 1 : -1;
    $read = $data[4]."\t".$data[0];
    $sr{$read} = $direction;
   }

while ($line = <SITE>)
   {
    $count ++; 
    chomp($line);
    next if ($line !~ /$tetag/);
    @data = split("\t", $line);
    $c1 = ($data[5] eq '+') ? 1 : -1;
    $c2 = ($data[6] eq '+') ? 1 : -1;
    $c3 = ($data[15] eq '+') ? 1 : -1;
    $direction[$count] = -1 * $c1 * $c2 * $c3;
   }

while ($line = <CALL>)
   {
    chomp($line);
    @data = split("\t", $line);
    @data1 = split(/\;/, $data[7]);
    $count_pos = $count_neg = 0;

    for ($j=0; $j<=$#data1; $j++)
       {
        $pe = $data1[$j];
        if ($pe =~ /^pe/)
          {
           @data2 = split(',', $pe);
           for ($i=4; $i<=$#data2; $i++)
              {
               @data3 = split('_', $data2[$i]);
               $id = $data3[0];
               if ($direction[$id] > 0)
                  {
                   $count_pos ++;
                  }
               else
                  {
                   $count_neg ++;
                  }
              }
          }
        elsif ($pe =~ /^sr/)
          {
           @data2 = split(',', $pe);
           @data3 = split(/\|/, $data2[2]);
           for ($i=0; $i<=$#data3; $i++)
              {
               $id = $data3[$i]."\t".$data[0];
               if ($sr{$id} > 0)
                  {
                   $count_pos ++;
                  }
               else
                  {
                   $count_neg ++;
                  }
              }
          }
       }
    print DIRE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$count_pos\t$count_neg\n";
   }
 
