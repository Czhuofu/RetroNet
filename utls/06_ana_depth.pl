#! /usr/bin/perl -w

use strict;

my $sample = $ARGV[0];
my $path   = $ARGV[1];
my $retro  = $ARGV[2];
my $TE     = $ARGV[3];

my $in_file = $path.'/'.$retro.'/'.$TE.'/'.$sample.'.'.$TE.'.novel.sites';
my $ou_file = $path.'/'.$retro.'/'.$TE.'/'.$sample.'.'.$TE.'.PE.calls';

my $Lrange  = -25;
my $Hrange  = 600;
my ($line, @data);
my $count = 0;
my ($range1, $read, $range2, $chr, $i, $sum, $newmin, $newmax, $min, $max);
my (%site, @chrsite, %cluster);
#(open (SR, "<$sr_file")) || die "cannot open the sr file\n";
(open (IN, "<$in_file")) || die "cannot open the in file\n";
(open (OU, ">$ou_file")) || die "cannot open the ou file\n";

while ($line = <IN>)
   {
    chomp($line);
    $count ++;
    print "read $count lines\n" if (!($count % 100000));
    @data = split("\t", $line);
    ### direction of the reads ###
    ### ++ ###
    ### -- ###
    ### +- ###
    ### -+ ###
    if ($data[5] eq '-')
       {
        $range1 = $data[1]-$Hrange;
        $range2 = $data[1]-$Lrange;
       }
    elsif ($data[5] eq '+')
       {
        $range1 = $data[2] + $Lrange;
        $range2 = $data[2] + $Hrange;
       }
    $chr    = $data[0];
    $site{$chr}{$range1}{$count} = 1;
   }

for $chr(keys %site)
   {
    @chrsite = ();
    $i = 0;
    for $range1(sort {$a<=>$b} keys %{$site{$chr}})
       {
        $chrsite[$i++]  = $range1;
       }
    #print "$chr\t$i\n";
    #next;
   ### merge reads ###
   for ($i=0; $i<=$#chrsite; $i++)
      {
       ### looking for sites to merge ###

       $min = $chrsite[$i];
       $max = $min + $Hrange - $Lrange;
       %cluster = ();
       for $count(keys %{$site{$chr}{$chrsite[$i]}})
          {
           $cluster{$count} = 1;
          }
#print "checking $i\n";
       if ($i<$#chrsite)
          {
           $newmin = $chrsite[$i+1];
           $newmax = $newmin + $Hrange - $Lrange;
           while ($newmin <= $max)
               {
                $i ++;
                for $count(keys %{$site{$chr}{$chrsite[$i]}})
                  {
                   $cluster{$count} = 1;
                  }
                $max = $newmax;
                last if ($i == $#chrsite);
                $newmin = $chrsite[$i+1];
                $newmax = $newmin + $Hrange - $Lrange;      
              } 
          }

       ### finish one cluster ###
       $min = 1 if ($min < 1);
       print OU "$chr\t$min\t$max\t$TE";
       $sum = scalar keys %cluster;
       print OU "\t$sum";
       for $count(keys %cluster)
          {
           print OU "\t$count";
          }
       print OU "\n";
      } ### finish one chr ###
   } ### finish all chrs ### 
          

 
